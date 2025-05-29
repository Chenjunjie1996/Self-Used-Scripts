import argparse
import pysam
import pandas as pd

class Filter_bam:
    """
    """
    def __init__(self, args):
        self.args = args
        
        # in
        self.bam = args.bam
        self.bclist = args.bclist
        
        # out
        self.set3 = set()
        self.set5 = set()
        self.total_count = 0
        self.count3 = 0
        self.count5 = 0
        self.intersec_count = 0
        self.count_file = "bc_umi_count.txt"
        self.set3_file = open("set3p.txt", 'w')
        self.set5_file = open("set5p.txt", 'w')
    
    
    def __call__(self):
        df = pd.read_csv(self.bclist,names=['bc'])
        cell_bc = set(df.bc)
        inbam = pysam.AlignmentFile(self.bam, "rb")
        
        for read in inbam:
            cb = read.get_tag("CB")
            umi = read.get_tag("UB")
            judge = read.query_name.split(':')[-1]
            self.total_count += 1
            if judge == '5p':
                self.count5 += 1
                if cb in cell_bc:
                    self.set5.add( (cb, umi) )
            else:
                self.count3 += 1
                if cb in cell_bc:
                    self.set3.add( (cb, umi) )
        self.intersec_count = len(self.set5.intersection(self.set3))
                
        inbam.close()
        
        with open(self.count_file, 'w') as fp:
            fp.write(f"total read count : {self.total_count}\n")
            fp.write(f"3p read count : {self.count3}\n")
            fp.write(f"5p read count : {self.count5}\n")
            fp.write(f"3p cell bc_umi combination : {len(self.set3)}\n")
            fp.write(f"5p cell bc_umi combination : {len(self.set5)}\n")
            fp.write(f"3p5p cell bc_umi intersection : {self.intersec_count}\n")
        
        
        for i in self.set3:
            self.set3_file.write(f"{':'.join(i)}\n")
            
        for i in self.set5:
            self.set5_file.write(f"{':'.join(i)}\n")
                    
        self.set3_file.close()
        self.set5_file.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--bam", help="bam file", required=True)
    parser.add_argument("--bclist", help="cell barcode file", required=True)
    args = parser.parse_args()
    Filter_bam(args)()