import pysam
import os
import argparse
from xopen import xopen


class Del30bp:
    
    def __init__(self, args):

        self.fq1 = args.fq1
        self.fq2 = args.fq2
        self.outdir = args.outdir
        self.out_fq2 = f"{self.outdir}/{self.fq2.split('/')[-1]}"
        

    def __call__(self, *args, **kwargs):
        self.run()


    def run(self):
        if not os.path.exists(self.outdir):
            os.system(f"mkdir -p {self.outdir}")

        os.system(f"cp {self.fq1} {self.outdir}")
        
        fh_fq2 = xopen(self.out_fq2, 'w')
        with pysam.FastxFile(self.fq2, persist=False) as fq2:
            for read in fq2:
                name = read.name
                comment = read.comment
                seq = read.sequence
                qual = read.quality

                fh_fq2.write(f"@{name} {comment}\n{seq[30:]}\n+\n{qual[30:]}\n")
        
        fh_fq2.close()
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="del 30bp of R2")
    parser.add_argument("--fq1", help="raw fq1 path", required=True)
    parser.add_argument("--fq2", help="raw fq2 path", required=True)
    parser.add_argument("--outdir", help="output_dir", default='del_30bp_r2')
    args = parser.parse_args()
    Del30bp(args)()