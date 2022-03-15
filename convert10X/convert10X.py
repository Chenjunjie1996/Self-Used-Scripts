import pandas as pd
import pysam
import os
from xopen import xopen


BARCODES_10X_FILE = "/SGRNJ/Database/script/soft/cellranger/cellranger-3.0.2/cellranger-cs/3.0.2/lib/python/cellranger/barcodes/737K-august-2016.txt"
UMI_10X_LEN = 10
TSO = "TTTCTTATATGGG"
SEQ_LEN = 150


def fastq_line(name, seq, qual):
    return f'@{name}\n{seq}\n+\n{qual}\n'


def convert_seq(sgr_barcode, umi, barcode_dict, barcodes_10X, seq2, qual2):
    '''
    barcode_dict - key:SGR barcode; value:10X barcode
    '''
    if sgr_barcode in barcode_dict:
        barcode_10X = barcode_dict[sgr_barcode]
    else:
        barcode_10X = barcodes_10X.readline().strip()
        barcode_dict[sgr_barcode] = barcode_10X

    if len(umi) > UMI_10X_LEN:
        umi_10X = umi[0:UMI_10X_LEN]
    elif len(umi) < UMI_10X_LEN:
        umi_10X = umi + 'C' * (UMI_10X_LEN - len(umi))
    else:
        umi_10X = umi

    seq2_insert = 90
    seq2_cut = 60

    new_seq2_1 = seq2[0:seq2_insert]
    new_seq2_2 = seq2[seq2_cut:]

    new_seq1 = barcode_10X + umi_10X + TSO
    new_qual1 = 'J' * len(new_seq1)
    new_qual2_1 = qual2[0:len(new_seq2_1)]
    new_qual2_2 = qual2[seq2_cut:]

    return new_seq1, new_qual1, new_seq2_1, new_qual2_1, new_seq2_2, new_qual2_2
    
class Convert(Step):
    """
    Features

    - Format barcodes and UMIs.

    Output
    - `02.convert/barcode_correspond.txt` Recording barcodes correspondence.

    - `02.convert/{sample}_S1_L001_R1_001.fastq.gz` New R1 reads.

    - `02.convert/{sample}_S1_L001_R2_001.fastq.gz` New R2 reads.
    """

    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)

        # common parameter
        self.outdir = args.outdir

        # input
        self.fq2 = args.fq2
        self.not_split_R2 = args.not_split_R2

        # output
        self.out_fq1_file = f'{self.outdir}/{self.sample}_S1_L001_R1_001.fastq.gz'
        self.out_fq2_file = f'{self.outdir}/{self.sample}_S1_L001_R2_001.fastq.gz'
        self.barcode_correspondence_file = f'{self.outdir}/barcode_correspond.txt'
        self.BARCODES_10X_FILE = os.path.dirname(
            soft_dict[args.soft]) + "/lib/python/cellranger/barcodes/737K-august-2016.txt"
        if args.soft_path:
            self.BARCODES_10X_FILE = os.path.dirname(
                args.soft_path) + "/lib/python/cellranger/barcodes/737K-august-2016.txt"

    @utils.add_log
    def run_convert(self):
        # read file
        barcodes_10X = open(self.BARCODES_10X_FILE, 'r')
        fq_file = pysam.FastxFile(self.fq2)

        # open out file
        out_fq1 = xopen(self.out_fq1_file, 'w')
        out_fq2 = xopen(self.out_fq2_file, 'w')

        # define var
        barcode_dict = {}

        # write out put fastq file
        for entry in fq_file:
            name = entry.name
            attrs = name.split('_')
            barcode = attrs[0]
            umi = attrs[1]
            seq = entry.sequence
            qual = entry.quality
            new_seq1, new_qual1, new_seq2_1, new_qual2_1, new_seq2_2, new_qual2_2 = convert_seq(barcode, umi,
                                                                                                barcode_dict,
                                                                                                barcodes_10X, seq, qual)

            if not self.not_split_R2:
                out_fq1.write(fastq_line(f'{name}_1', new_seq1, new_qual1))
                out_fq1.write(fastq_line(f'{name}_2', new_seq1, new_qual1))
                out_fq2.write(fastq_line(f'{name}_1', new_seq2_1, new_qual2_1))
                out_fq2.write(fastq_line(f'{name}_2', new_seq2_2, new_qual2_2))
            else:
                out_fq1.write(fastq_line(name, new_seq1, new_qual1))
                out_fq2.write(fastq_line(name, seq, qual))

        out_fq1.close()
        out_fq2.close()
        barcodes_10X.close()

        # os.system(f'gzip {out_fq1}')
        # os.system(f'gzip {out_fq2}')
        # write barcode correspond file
        barcode_record = pd.DataFrame()
        barcode_record['sgr'] = list(barcode_dict.keys())
        barcode_record['10X'] = [barcode_dict[i] for i in barcode_dict]
        barcode_record.to_csv(self.barcode_correspondence_file, sep='\t', index=False)