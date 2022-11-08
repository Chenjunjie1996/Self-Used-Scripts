import celescope
import os
from celescope.tools import utils
from concurrent.futures import ProcessPoolExecutor
import subprocess
import argparse
import pysam
import json


TOOLS_DIR = os.path.dirname(celescope.tools.__file__) + '/trust4'
REF_DIR = f'{TOOLS_DIR}/database'


OUT_NAME = {
    "TCR": ["bcrtcr", "TRA", "TRB"],
    "BCR": ["bcrtcr", "IGH", "IGK", "IGL"],
}


parser = argparse.ArgumentParser(description='get mapping vdj reads number')
parser.add_argument("--fq1", help="raw fastq1 file path", required=True)
parser.add_argument("--fq2", help="raw fastq1 file path", required=True)
parser.add_argument("--ref", help="reference name: hg38 or GRCm38", choices=["hg38", "GRCm38"], required=True)
parser.add_argument("--seqtype", help="TCR or BCR", choices=["TCR", "BCR"], required=True)
args = parser.parse_args()


def get_fastx_read_number(fastx_file):
    """
    get read number using pysam
    """
    n = 0
    with pysam.FastxFile(fastx_file) as f:
        for _ in f:
            n += 1
    return n
    

def dump_dict_to_json(d, json_file):
    with open(json_file, 'w') as f:
        json.dump(d, f, indent=4)
        
        
class Mapping_vdj:
    """
    ##Features

    - Get Mapping to VDJ reads number.
    
    For TCR library:
    Total Reads: 10,000
    Reads Mapped To Any V(D)J Genes: 9,000
    Reads Mapped To TRA: 4,000
    Reads Mapped To TRB: 5,000

    For BCR library:
    Total Reads: 20,000
    Reads Mapped To Any V(D)J Genes: 15,000
    Reads Mapped To IGH: 5,000
    Reads Mapped To IGK: 5,000
    Reads Mapped To IGL: 5,000

    """

    def __init__(self, args, result=[]):
        self.ref = args.ref
        self.fq1 = args.fq1
        self.fq2 = args.fq2
        self.seqtype = args.seqtype
        self.result_list = result
    
    def __call__(self):
#        with ProcessPoolExecutor(max_workers=4) as executor:
#            for result in executor.map(self.get_mapping_reads, OUT_NAME[self.seqtype]):
#                self.result_list.append(result)
        
        self.write_metrics()

    @utils.add_log
    def get_mapping_reads(self, out_name):
        
        cmd = (
            f"fastq-extractor -t 2 "
            f"-f {REF_DIR}/{self.ref}/{out_name}.fa "
            f"-o {out_name} "
            f"-1 {self.fq1} -2 {self.fq2} "
        )
        Mapping_vdj.get_mapping_reads.logger.info(cmd)
        subprocess.check_call(cmd, shell=True)

    @utils.add_log
    def write_metrics(self):
        metrics = []
        with ProcessPoolExecutor(max_workers=4) as executor:
            for result in executor.map(get_fastx_read_number, [self.fq1] + [i + "_1.fq" for i in OUT_NAME[self.seqtype]]):
                metrics.append(result)
        
        if self.seqtype == "TCR":
            metrics = dict(zip(["Total Reads" ,"Reads Mapped To Any V(D)J Genes", "Reads Mapped To TRA", "Reads Mapped To TRB"], metrics))
        else:
            metrics = dict(zip(["Total Reads" ,"Reads Mapped To Any V(D)J Genes", "Reads Mapped To IGH", "Reads Mapped To IGK", "Reads Mapped To IGL"], metrics))

        dump_dict_to_json(metrics, "./Mapping_VDJ.json")
    
    
if __name__ == '__main__':
    runner = Mapping_vdj(args)
    runner()