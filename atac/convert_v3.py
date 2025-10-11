from xopen import xopen
from datetime import timedelta
from functools import wraps
from collections import defaultdict
import subprocess
import pysam
import random
import os
import operator
import argparse
import logging
import time
import sys
import glob
import json
import pandas as pd


# Const
WHITELIST_ATAC = "/SGRNJ06/randd/USER/cjj/Script/Self-Used-Scripts/atac/737K-arc-v1-atac.txt.gz"
WHITELIST_RNA = "/SGRNJ06/randd/USER/cjj/Script/Self-Used-Scripts/atac/737K-arc-v1-rna.txt.gz"
SGR_ATAC_RNA =  "/SGRNJ06/randd/USER/cjj/Script/Self-Used-Scripts/atac/sgr-atac-rna-V3.txt"
UMI_10X_LEN = 12

# Template Swithcing Oligos Sequence. [16-bp cell barcode][10/12-bp UMI][TSO].
TSO = "TTTCTTATATGGG"


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


def dump_dict_to_json(d, json_file):
    with open(json_file, 'w') as f:
        json.dump(d, f, indent=4)
        

class Convert:
    """Convert barcodes to 10X format."""
    def __init__(self, args):
        self.args = args
        self.atac_fq_list = glob.glob(f"{self.args.atac_path}/*.fastq")
        self.rna_fq_list = glob.glob(f"{self.args.rna_path}/*.fq")
        self.atac_fq_list = sorted([i for i in self.atac_fq_list if os.path.isfile(i)])
        self.rna_fq_list = sorted([i for i in self.rna_fq_list if os.path.isfile(i)])
        assert len(self.atac_fq_list) == 3
        assert len(self.rna_fq_list) == 1

        self.whitelist_10X_rna = xopen(WHITELIST_RNA, 'r')
        self.sgr_tenX = {}

        self.atac_outdir = f"{self.args.outdir}/atac"
        self.rna_outdir = f"{self.args.outdir}/rna"
        for outdir in [self.atac_outdir, self.rna_outdir]:
            if not os.path.exists(outdir):
                os.system(f"mkdir -p {outdir}")
        
        self.atac_fq1, self.atac_fq2, self.atac_fq3 = self.atac_fq_list[0], self.atac_fq_list[1], self.atac_fq_list[2]
        self.rna_fq2 = self.rna_fq_list[0]
        self.sample = self.args.sample

        self.atac_out_fq1 = f'{self.atac_outdir}/{self.sample}_S1_L001_R1_001.fastq.gz'
        self.atac_out_fq2 = f'{self.atac_outdir}/{self.sample}_S1_L001_R2_001.fastq.gz'
        self.atac_out_fq3 = f'{self.atac_outdir}/{self.sample}_S1_L001_R3_001.fastq.gz'
        self.rna_out_fq1 = f'{self.rna_outdir}/{self.sample}_S1_L001_R1_001.fastq.gz'
        self.rna_out_fq2 = f'{self.rna_outdir}/{self.sample}_S1_L001_R2_001.fastq.gz'
        self.barcode_convert_json = f'{self.rna_outdir}/barcode_convert.json'
    
    @add_log
    def gen_sgr_tenX_dict(self):
        
        count_dict = defaultdict(int)

        with pysam.FastxFile(self.rna_fq2) as fq2_fh:
            for entry in fq2_fh:
                sgr_barcode = entry.name.split(':')[0]
                count_dict[sgr_barcode] += 1

        count_dict = dict(sorted(count_dict.items(), key=operator.itemgetter(1), reverse=True))

        for sgr_barcode in count_dict:
            self.sgr_tenX[sgr_barcode] = self.whitelist_10X_rna.readline().strip()

        # Add invalid barcode
        for sgr_barcode, barcode_10X in self.sgr_tenX.items():
            if barcode_10X == '':
                self.sgr_tenX[sgr_barcode] = "AAAA" + ''.join(random.choice("ATCG") for _ in range(12))

    @add_log
    def gen_atac_rna_dict(self):
        df_sgr_atac_rna = pd.read_csv(SGR_ATAC_RNA, sep="\t", header=None, names=["atac","rna"])
        # Celescope V2
        df_sgr_atac_rna['rna'] = df_sgr_atac_rna['rna'].apply(lambda x: x[:9] + '_' + x[9:18] + '_' + x[18:])
        self.sgr_atac_rna = df_sgr_atac_rna.set_index("atac").to_dict(orient="dict")["rna"]
        
        df_tenX_rna = pd.read_csv(WHITELIST_RNA, sep="\t", header=None, names=["rna"])
        df_tenX_atac = pd.read_csv(WHITELIST_ATAC, sep="\t", header=None, names=["atac"])
        df_tenX_rna_atac = pd.concat([df_tenX_rna, df_tenX_atac], axis=1)
        self.tenX_rna_atac = df_tenX_rna_atac.set_index("rna").to_dict(orient="dict")["atac"]
        
    @add_log
    def gzip_fq(self):
        cmd1 = f"gzip -c {self.rna_fq2} > {self.rna_out_fq2}"
        cmd2 = f"gzip -c {self.atac_fq1} > {self.atac_out_fq1}"
        cmd3 = f"gzip -c {self.atac_fq3} > {self.atac_out_fq3}"
        
        for cmd in [cmd1, cmd2, cmd3]:
            subprocess.check_call(cmd, shell=True)

    @add_log
    def write_rna_fq(self):
        out_fq1 = xopen(self.rna_out_fq1, 'w')

        with pysam.FastxFile(self.rna_fq2) as fq2_fh:
            for entry in fq2_fh:
                name = entry.name
                attrs = name.split(':')
                sgr_barcode, sgr_umi = attrs[0], attrs[1]
                new_seq1, new_qual1 = self.convert_seq(sgr_barcode, sgr_umi)
                out_fq1.write(f'@{name}\n{new_seq1}\n+\n{new_qual1}\n')

        out_fq1.close()
    
    @add_log
    def write_atac_fq(self):
        out_fq2 = xopen(self.atac_out_fq2, 'w')
        
        with pysam.FastxFile(self.atac_fq2) as fq2_fh:
            for entry in fq2_fh:
                name = entry.name
                sgr_barcode_atac = entry.sequence
                sgr_barcode_rna = self.sgr_atac_rna[sgr_barcode_atac]
                barcode_10X_rna = self.sgr_tenX.get(sgr_barcode_rna, "A"*16)
                barcode_10X_atac = self.tenX_rna_atac.get(barcode_10X_rna, "A"*16)
                new_qual = 'F' * len(barcode_10X_atac)
                out_fq2.write(f'@{name}\n{barcode_10X_atac}\n+\n{new_qual}\n')

        out_fq2.close()

    def convert_seq(self, barcode_sgr, umi_sgr):
        """
        Convert sgr barcode to 10X barcode; change length of sgr UMI to UMI_10X_LEN

        Args:
            barcode_sgr: str
            umi_sgr: str
        Returns:
            new_seq1: str
            new_qual1: str
        """
        barcode_10X = self.sgr_tenX[barcode_sgr]

        umi_len_sgr = len(umi_sgr)
        if umi_len_sgr > UMI_10X_LEN:
            umi_10X = umi_sgr[:UMI_10X_LEN]
        elif umi_len_sgr < UMI_10X_LEN:
            umi_10X = umi_sgr + 'C' * (UMI_10X_LEN - umi_len_sgr)
        else:
            umi_10X = umi_sgr

        new_seq1 = barcode_10X + umi_10X + TSO
        new_qual1 = 'F' * len(new_seq1)

        return new_seq1, new_qual1

    @add_log
    def dump_tenX_sgr_barcode_json(self):
        tenX_sgr = {}
        for sgr, tenX in self.sgr_tenX.items():
            tenX_sgr[tenX] = sgr

        dump_dict_to_json(tenX_sgr, self.barcode_convert_json)

    def run(self):
        self.gzip_fq()
        self.gen_atac_rna_dict()
        self.gen_sgr_tenX_dict()
        self.write_rna_fq()
        self.write_atac_fq()
        self.dump_tenX_sgr_barcode_json()

    def __call__(self, *args, **kwargs):
        self.run()
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert fastq to 10X")
    parser.add_argument("--atac_path", help="atac_path", required=True)
    parser.add_argument("--rna_path", help="atac_path", required=True)
    parser.add_argument("--sample", help="sample name", required=True)
    parser.add_argument("--outdir", help="output_dir", default='convert_dir')
    args = parser.parse_args()
    Convert(args)()