from concurrent.futures import ProcessPoolExecutor
import pysam
import os
import glob
import argparse
import pandas as pd


class Add_Bp:
    def __init__(self, args):
        self.mapfile = args.mapfile
    
    def parse_mapfile(self):
        df = pd.read_csv(self.mapfile, sep='\t', header=None)
        df.columns = ['library_id', 'rawdata_path', 'sample_name', 'match_dir']
        df['raw_fq'] = df[['rawdata_path', 'library_id']].apply('/'.join, axis=1)

        fq_list, sample_list = list(df['raw_fq']), list(df['sample_name'])
        
        # 原始R1数据 用于提取末尾33bp
        fq_list = [glob.glob(f'{i}_*R1.fastq.gz')[0] for i in fq_list]
        print(fq_list)
        if not fq_list:
            raise FileNotFoundError("check raw fastq file")
        
        # 150bp的原始convert fq文件
        convert_fq_list = [f'./{sample}/02.convert/{sample}_S1_L001_R1_001.fastq.gz' for sample in sample_list]
        # tmp fastq记录R1序列最后33bp序列
        tmp_fq_list = [f'./{sample}/02.convert/tmp.fq' for sample in sample_list]
        # 拼接R1 33bp序列的 R2序列
        out_convert_fq_list = [f'./{sample}/02.convert/{sample}_S1_L001_R1_001.fastq' for sample in sample_list]
        
        return fq_list, convert_fq_list, tmp_fq_list, out_convert_fq_list
    
    @staticmethod
    def add_33bp_2r2(raw_fq, convert_fq, tmp_fq, out_convert_fq):
        # 找出序列index
        index_set = set()
        with pysam.FastxFile(convert_fq) as f:
            for read in f:
                index_set.add(read.name.split('_')[-1])

        # tmp文件记录目标序列最后33bp序列
        tmp = open(tmp_fq, 'w')
        total_num = 0
        with pysam.FastxFile(raw_fq) as f:
            for read in f:
                total_num += 1
                if str(total_num) in index_set:
                    last33bp = read.sequence[-33:]
                    newqual = len(last33bp) * 'F'
                    tmp.write(f"@{read.name}\n{last33bp}\n+\n{newqual}\n")
        tmp.close()

        # 拼接序列
        out = open(out_convert_fq, 'w')
        with pysam.FastxFile(tmp_fq, persist=False) as fq1, \
                pysam.FastxFile(convert_fq, persist=False) as fq2:
            for entry1, entry2 in zip(fq1, fq2):
                header1, seq1, qual1 = entry1.name, entry1.sequence, entry1.quality
                header2, seq2, qual2 = entry2.name, entry2.sequence, entry2.quality
            
                last33bp = seq1
                newseq = seq2 + last33bp
                newqual = len(newseq) * 'F'
                out.write(f"@{header2}\n{newseq}\n+\n{newqual}\n")
        out.close()

        os.system(f"rm {tmp_fq}")
        os.system(f"rm {convert_fq}")
        os.system(f"gzip -c {out_convert_fq} > {out_convert_fq}.gz")
        os.system(f"rm {out_convert_fq}")

    def __call__(self):
        
        fq_list, convert_fq_list, tmp_fq_list, out_convert_fq_list = self.parse_mapfile()
        
        with ProcessPoolExecutor(max_workers=len(fq_list)) as executor:
            for result in executor.map(
                self.add_33bp_2r2, fq_list, convert_fq_list, tmp_fq_list, out_convert_fq_list
            ):
                print(result, 'done')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Add 33bp R1 reads to R2")
    parser.add_argument("--mapfile", help="mapfile", required=True)
    args = parser.parse_args()
    Add_Bp(args)()
