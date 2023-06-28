import pysam
import glob
import sys


def main(fq_dir):
    # 记录index
    valid_fq = glob.glob("*_2.fq")[0]
    index_set = set()
    with pysam.FastxFile(valid_fq) as f:
        for read in f:
            index_set.add(read.name.split('_')[-1])
    
    invalid_fq1 = open("invalid_1.fq", 'w')
    invalid_fq2 = open("invalid_2.fq", 'w')
    
    raw_fq1 = glob.glob(f'{fq_dir}/*_R1.fastq.gz')[0]
    raw_fq2 = glob.glob(f'{fq_dir}/*_R2.fastq.gz')[0]
    
    total_num = 0
    with pysam.FastxFile(raw_fq1, persist=False) as fq1, \
            pysam.FastxFile(raw_fq2, persist=False) as fq2:
        for entry1, entry2 in zip(fq1, fq2):
            total_num += 1
            if str(total_num) not in index_set:
                invalid_fq1.write(f"@{entry1.name}\n{entry1.sequence}\n+\n{entry1.quality}\n")
                invalid_fq2.write(f"@{entry2.name}\n{entry2.sequence}\n+\n{entry2.quality}\n")
    
    invalid_fq1.close()
    invalid_fq2.close()


if __name__ == '__main__':
    fq_dir = sys.argv[1]
    main(fq_dir)