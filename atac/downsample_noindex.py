# DOWNSAMPLE
import pysam
import argparse
import os
import glob
import logging
import time
import sys
from datetime import timedelta
from functools import wraps
from xopen import xopen


def add_log(func):
    '''
    logging start and done.
    '''
    logFormatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    module = func.__module__
    name = func.__name__
    logger_name = f'{module}.{name}'
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)

    consoleHandler = logging.StreamHandler(sys.stdout)
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)

    @wraps(func)
    def wrapper(*args, **kwargs):

        logger.info('start...')
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        used = timedelta(seconds=end - start)
        logger.info('done. time used: %s', used)
        return result

    wrapper.logger = logger
    return wrapper


class Downsample:
    
    def __init__(self, args):
        self.args = args

        self.fq_list = glob.glob(f"{self.args.fq_path}/*_R*.fastq*")
        self.fq_list = sorted([i for i in self.fq_list if os.path.isfile(i)])
        
        print(f" Downsample List : \n{self.fq_list}")

        assert len(self.fq_list) == 3
        
        self.fq1, self.fq2, self.fq3 = self.fq_list[0], self.fq_list[1], self.fq_list[2]
        self.out_fq1 = f"{self.args.outdir}/downsample_{self.fq1.split('/')[-1]}"
        self.out_fq2 = f"{self.args.outdir}/downsample_{self.fq2.split('/')[-1]}"
        self.out_fq3 = f"{self.args.outdir}/downsample_{self.fq3.split('/')[-1]}"
        

    def __call__(self, *args, **kwargs):
        self.run()

    @add_log
    def run(self):
        if not os.path.exists(self.args.outdir):
            os.system(f"mkdir -p {self.args.outdir}")

        fh_fq1 = xopen(self.out_fq1, 'w')
        fh_fq2 = xopen(self.out_fq2, 'w')
        fh_fq3 = xopen(self.out_fq3, 'w')

        count = 0
        with pysam.FastxFile(self.fq1, persist=False) as fq1, \
                pysam.FastxFile(self.fq2, persist=False) as fq2, \
                    pysam.FastxFile(self.fq3, persist=False) as fq3:

            for entry1, entry2, entry3 in zip(fq1, fq2, fq3):
                count += 1
                if count > self.args.reads_num:
                    break
                fh_fq1.write(f"@{entry1.name}\n{entry1.sequence}\n+\n{entry1.quality}\n")
                fh_fq2.write(f"@{entry2.name}\n{entry2.sequence}\n+\n{entry2.quality}\n")
                fh_fq3.write(f"@{entry3.name}\n{entry3.sequence}\n+\n{entry3.quality}\n")

                if count % 1000000 == 0:
                    self.run.logger.info(f'Downsample {count} reads done.')


            self.run.logger.info(self.fq1 + ' finished.')


        fh_fq1.close()
        fh_fq2.close()
        fh_fq3.close()
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Downsample fastq file")
    parser.add_argument("--fq_path", help="fastq_path", required=True)
    parser.add_argument("--outdir", help="output_dir", default='Downsample')
    parser.add_argument("--reads_num", help="downsample reads number", default=10000000, type=int)
    args = parser.parse_args()
    Downsample(args)()