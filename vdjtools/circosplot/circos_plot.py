import pandas as pd
import argparse
import os
import glob
import logging
import time
import sys
import subprocess
from datetime import timedelta
from functools import wraps

ROOT_DIR = os.path.dirname(__file__)
CHAIN = {
    'TCR': ['TRA', 'TRB'],
    'BCR': ['IGH', 'IGL', 'IGK']
}


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


class Circos:
    """ PlotFancyVJUsage of cellranger output
    """

    def __init__(self, args):
        self.seqtype = args.seqtype
        self.outdir = args.outdir
        self.split_pair = args.split_pair
        self.consensus = glob.glob(f"{args.sample_path}/03.assemble/*/outs/consensus_annotations.csv")[0]
        self.clonotype = glob.glob(f"{args.sample_path}/03.assemble/*/outs/clonotypes.csv")[0]
        self.sample = self.clonotype.split('/')[-3]

    def __call__(self, *args, **kwargs):
        self.convert_format()
        os.chdir(self.outdir)
        self.calc_vj()
        self.make_plot()

    @staticmethod
    def check_mkdir(dir_name):
        if not os.path.exists(dir_name):
            os.system(f"mkdir -p {dir_name}")

    @add_log
    def convert_format(self):
        Circos.check_mkdir(self.outdir)

        df_consensus = pd.read_csv(self.consensus)

        df_clonotype = pd.read_csv(self.clonotype)
        df_clonotype = df_clonotype[["clonotype_id", "frequency", "proportion"]]
        df_clonotype.rename(columns={"frequency": "count", "proportion": "freq"}, inplace=True)

        df_consensus = df_consensus[df_consensus["productive"] == True]
        df_consensus = df_consensus[["clonotype_id", "chain", "v_gene", "d_gene", "j_gene", "cdr3", "cdr3_nt"]]
        df_consensus.loc[df_consensus.d_gene == "None", "d_gene"] = '.'
        df_consensus.rename(
            columns={"v_gene": "v", "d_gene": "d", "j_gene": "j", "cdr3_nt": "cdr3nt", "cdr3": "cdr3aa"}, inplace=True)

        if not self.split_pair:
            df_tmp = pd.merge(df_consensus, df_clonotype)
            df_tmp.sort_values(by="count", ascending=False, inplace=True)
            df_tmp = df_tmp[["count", "freq", "cdr3nt", "cdr3aa", "v", "d", "j"]]
            df_tmp.to_csv(f"{self.outdir}/{self.sample}.txt", sep='\t', index=False)
        else:
            for chain in CHAIN[self.seqtype]:
                df_tmp = df_consensus[df_consensus["chain"] == chain]
                df_tmp = pd.merge(df_tmp, df_clonotype)
                df_tmp.sort_values(by="count", ascending=False, inplace=True)
                df_tmp = df_tmp[["count", "freq", "cdr3nt", "cdr3aa", "v", "d", "j"]]
                df_tmp.to_csv(f"{self.outdir}/{self.sample}_{chain}.txt", sep='\t', index=False)

    @add_log
    def calc_vj(self):
        if not self.split_pair:
            cmd = (
                f"java -jar /SGRNJ06/randd/USER/cjj/soft/vdjtools/vdjtools-1.2.1/vdjtools-1.2.1.jar PlotFancyVJUsage "
                f"{self.sample}.txt "
                f"{self.sample} "
            )
            subprocess.check_call(cmd, shell=True)
        else:
            for chain in CHAIN[self.seqtype]:
                cmd = (
                    f"java -jar /SGRNJ06/randd/USER/cjj/soft/vdjtools/vdjtools-1.2.1/vdjtools-1.2.1.jar PlotFancyVJUsage "
                    f"{self.sample}_{chain}.txt "
                    f"{self.sample}_{chain} "
                )
                subprocess.check_call(cmd, shell=True)

    @add_log
    def make_plot(self):
        if not self.split_pair:
            cmd = (
                f"Rscript {ROOT_DIR}/vj_pairing_plot.r "
                f"{self.sample}.fancyvj.wt.txt "
                f"{self.sample}.fancyvj.wt.png "
            )
            subprocess.check_call(cmd, shell=True)
        else:
            for chain in CHAIN[self.seqtype]:
                cmd = (
                    f"Rscript {ROOT_DIR}/vj_pairing_plot.r "
                    f"{self.sample}_{chain}.fancyvj.wt.txt "
                    f"{self.sample}_{chain}.fancyvj.wt.png "
                )
                subprocess.check_call(cmd, shell=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="PlotFancyVJUsage of flv_CR output")
    parser.add_argument("sample_path", help="absolute sample path of flv_CR")
    parser.add_argument("--seqtype", help="TCR or BCR", default='TCR')
    parser.add_argument("--outdir", help="output_dir", default='VJUsage')
    parser.add_argument('--split_pair', help='split pair', action='store_true')
    args = parser.parse_args()
    Circos(args)()
