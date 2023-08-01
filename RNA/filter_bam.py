import argparse
import pysam
from Bio.Seq import Seq


def reverse_complement(seq):
    """Reverse complementary sequence

    :param original seq
    :return Reverse complementary sequence
    """
    return str(Seq(seq).reverse_complement())


class Filter_bam:
    """
    读入3prime, 5prime aligned_posSorted_addTag.bam，把(barcode, umi, gene_id)不在3prime的reads过滤掉。
    Barcode, UMI 为反向互补，若gene_id不存在填入None
    """
    def __init__(self, args):
        self.args = args
        self.bam3 = args.bam3
        self.bam5 = args.bam5
        self.out_bam = args.out_bam
        self.prime3_set = set()
        self.filter_count = 0
    
    
    def __call__(self):
        inputFile3 = pysam.AlignmentFile(self.bam3, "rb")
        inputFile5 = pysam.AlignmentFile(self.bam5, "rb")
        out_file = pysam.AlignmentFile(self.out_bam, "wb", header=inputFile5.header)
        
        for read in inputFile3:
            cb = read.get_tag("CB")
            umi = read.get_tag("UB")
            try:
                gene_id = read.get_tag("XT")
            except KeyError:
                gene_id = None
            self.prime3_set.add( (cb, umi, gene_id) )
        
        for read in inputFile5:
            cb = read.get_tag("CB")
            umi = read.get_tag("UB")
            try:
                gene_id = read.get_tag("XT")
            except KeyError:
                gene_id = None
            if (reverse_complement(cb), reverse_complement(umi), gene_id) not in self.prime3_set:
                self.filter_count += 1
                out_file.write(read)
        
        with open("count.txt", 'w') as fp:
            fp.write(f"filted read count : {self.filter_count}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="filter 5-prime bam file")
    parser.add_argument("--bam3", help="3prime aligned_posSorted_addTag.bam", required=True)
    parser.add_argument("--bam5", help="5prime aligned_posSorted_addTag.bam", required=True)
    parser.add_argument("--out_bam", help="filtered aligned_posSorted_addTag.bam", default="aligned_posSorted_addTag_filtered.bam")
    args = parser.parse_args()
    Filter_bam(args)()