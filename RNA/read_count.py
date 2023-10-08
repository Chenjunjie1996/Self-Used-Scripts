from Bio.Seq import Seq 
import pysam 
import pandas as pd 
import argparse 

""" Count read number
5 prime featureCounts assigned reads
3 prime featureCounts assigned reads
5 prime reads in 3 prime CB(cell barcode)
3 prime reads in 3 prime CB
5 prime reads in 3 prime CB + UMI
5 prime reads in 3 prime CB + UMI + gene
"""


prime5_assign_reads = 0
prime3_assign_reads = 0
prime5_in_prime3_cb = 0
prime3_in_prime3_cb = 0
prime5_in_prime3_cb_umi = 0
prime5_in_prime3_cb_umi_gene = 0


def reverse_complement(seq):
    """Reverse complementary sequence

    :param original seq
    :return Reverse complementary sequence
    """
    return str(Seq(seq).reverse_complement())


def count_read_num(bam3, bam5, cb):
    
    prime3_cb = pd.read_csv(prime3_cb, header=None)
    prime3_cb = prime3_cb.rename(columns={0: "bc"})
    prime3_cb_set = set(prime3_cb.bc)
    prime3_cb_umi_set = set()
    prime3_cb_umi_gene_set = set()

    inputFile3 = pysam.AlignmentFile(bam3, "rb")
    inputFile5 = pysam.AlignmentFile(bam5, "rb")
        
    for read in inputFile3:
        stat = read.get_tag("XS")
        if stat == "Assigned":
            prime3_assign_reads += 1
            cb = read.get_tag("CB")
            umi = read.get_tag("UB")
            try:
                gene_id = read.get_tag("XT")
            except KeyError:
                gene_id = None
            if cb in prime3_cb_set:
                prime3_in_prime3_cb += 1
                prime3_cb_set.add( reverse_complement(cb) )
                prime3_cb_umi_set.add( (reverse_complement(cb), reverse_complement(umi)) )
                prime3_cb_umi_gene_set.add( (reverse_complement(cb), reverse_complement(umi), gene_id) )
        
    for read in inputFile5:
        stat = read.get_tag("XS")
        if stat == "Assigned":
            prime5_assign_reads += 1
            cb = read.get_tag("CB")
            umi = read.get_tag("UB")
            try:
                gene_id = read.get_tag("XT")
            except KeyError:
                gene_id = None
            if cb in prime3_cb_set:
                prime5_in_prime3_cb += 1
            if (cb, umi) in prime3_cb_umi_set:
                prime5_in_prime3_cb_umi += 1
            if (cb, umi, gene_id) in prime3_cb_umi_gene_set:
                prime5_in_prime3_cb_umi_gene += 1
    
    sample_name = bam5.split('/')[-3]
    out_file = open(f"{sample_name}.txt", 'w')
    out_file.write(f"5 prime featureCounts assigned reads : {format(prime5_assign_reads, ',')}\n")
    out_file.write(f"3 prime featureCounts assigned reads : {format(prime3_assign_reads, ',')}\n")
    out_file.write(f"5 prime reads in 3 prime CB : {format(prime5_in_prime3_cb, ',')}\n")
    out_file.write(f"3 prime reads in 3 prime CB : {format(prime3_in_prime3_cb, ',')}\n")
    out_file.write(f"5 prime reads in 3 prime CB + UMI : {format(prime5_in_prime3_cb_umi, ',')}\n")
    out_file.write(f"5 prime reads in 3 prime CB + UMI + gene : {format(prime5_in_prime3_cb_umi_gene, ',')}\n")
    out_file.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="filter 5-prime bam file")
    parser.add_argument("--bam3", help="3 prime aligned_posSorted_addTag.bam", required=True)
    parser.add_argument("--bam5", help="5 prime aligned_posSorted_addTag.bam", required=True)
    parser.add_argument("--cb", help="3 prime cell barcode", required=True)
    args = parser.parse_args()
    count_read_num(args.bam3, args.bam5, args.cb)

