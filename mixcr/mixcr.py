import time
import logging
import subprocess
import argparse
import os
import sys
from datetime import timedelta
from functools import wraps


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


@add_log
def run_mixcr(args, dir_name):
    cmd = (
        'mixcr align '
        '--force-overwrite '
        f'--species {args.species} '
        f'-t 8 '
        f'--not-aligned-R1 {dir_name}/not_align.fq '
        f'--report {dir_name}/align.txt '
        '-OallowPartialAlignments=true '
        '-OvParameters.geneFeatureToAlign=VTranscriptWithP '
        f'{args.fq} {dir_name}/read2.vdjca '
        f'-Xmx 8g '
        '\n'
        'mixcr exportAlignments '
        f'{dir_name}/read2.vdjca {dir_name}/alignments.txt '
        '-readIds --force-overwrite -vGene -dGene -jGene -cGene '
        '-nFeature CDR3 -aaFeature CDR3 '
        f'-Xmx 8g '
        )

    subprocess.check_call(cmd, shell=True)


if __name__ == '__main__':
    dir_name = "mixcr_res"
    if not os.path.exists(dir_name):
        os.system(f"mkdir -p {dir_name}")
    
    parser = argparse.ArgumentParser(description='mixcr cdr3')
    parser.add_argument('--fq', help='contig file', required=True)
    parser.add_argument('--species', help='human or mouse', required=True)
    args = parser.parse_args()
    
    run_mixcr(args, dir_name)