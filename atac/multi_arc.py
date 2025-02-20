from collections import defaultdict
import glob
import itertools
import os
import argparse


def check_mkdir(dir_name):
    """if dir_name is not exist, make one"""
    if not os.path.exists(dir_name):
        os.system(f"mkdir -p {dir_name}")
        

def get_read(library_id, library_path, read='1'):
    read1_list = [f'_{read}', f'R{read}', f'R{read}_001']
    fq_list = ['fq', 'fastq']
    suffix_list = ["", ".gz"]
    read_pattern_list = [
        f'{library_path}/{library_id}*{read}.{fq_str}{suffix}'
        for read in read1_list
        for fq_str in fq_list
        for suffix in suffix_list
    ]
    fq_list = [glob.glob(read1_pattern) for read1_pattern in read_pattern_list]
    fq_list = (non_empty for non_empty in fq_list if non_empty)
    fq_list = sorted(list(itertools.chain(*fq_list)))
    if len(fq_list) == 0:
        print("Allowed R1 patterns:")
        for pattern in read_pattern_list:
            print(pattern)
        raise Exception(
            '\n'
            f'Invalid Read{read} path! \n'
            f'library_id: {library_id}\n'
            f'library_path: {library_path}\n'
        )
    return fq_list


def get_fq(library_id, library_path, use_R3=False):
    fq1_list = get_read(library_id, library_path, read='1')
    fq2_list = get_read(library_id, library_path, read='2')
    if use_R3:
        fq3_list = get_read(library_id, library_path, read='3')
    if len(fq1_list) != len(fq2_list):
        raise Exception("Read1 and Read2 fastq number do not match!")
    fq1 = ",".join(fq1_list)
    fq2 = ",".join(fq2_list)
    if use_R3:
        fq3 = ",".join(fq3_list)
        return fq1, fq2, fq3
    return fq1, fq2

    
def parse_mapfile(mapfile):
    atac_fq_dict = defaultdict(list)
    rna_fq_dict = defaultdict(list)
    with open(mapfile) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            line_split = line.split()
            assert len(line_split) == 4
            library_id, library_path, sample_name, library_type = line_split
            
            if library_type == "rna":
                fq1, fq2 = get_fq(library_id, library_path)
                if sample_name in rna_fq_dict:
                    rna_fq_dict[sample_name][0].append(fq1)
                    rna_fq_dict[sample_name][1].append(fq2)
                else:
                    rna_fq_dict[sample_name] = [[fq1], [fq2]]
                    
            elif library_type == 'atac':
                fq1, fq2, fq3 = get_fq(library_id, library_path, use_R3=True)
                if sample_name in atac_fq_dict:
                    atac_fq_dict[sample_name][0].append(fq1)
                    atac_fq_dict[sample_name][1].append(fq2)
                    atac_fq_dict[sample_name][2].append(fq3)
                else:
                    atac_fq_dict[sample_name] = [[fq1], [fq2], [fq3]]

    for sample_name in rna_fq_dict:
        rna_fq_dict[sample_name][0] = ",".join(rna_fq_dict[sample_name][0])
        rna_fq_dict[sample_name][1] = ",".join(rna_fq_dict[sample_name][1])
            
    for sample_name in atac_fq_dict:
        atac_fq_dict[sample_name][0] = ",".join(atac_fq_dict[sample_name][0])
        atac_fq_dict[sample_name][1] = ",".join(atac_fq_dict[sample_name][1])
        atac_fq_dict[sample_name][2] = ",".join(atac_fq_dict[sample_name][2])

    if not atac_fq_dict:
        raise Exception('empty atac path!')
    if not rna_fq_dict:
        raise Exception('empty rna path!')
    
    return rna_fq_dict, atac_fq_dict


def prepare_sjm(mapfile, reference):
    """
    parse_mapfile, make log dir, init script variables, init outdir_dic
    make sjm dir, sjm file
    """
    sjm_dir = './sjm'
    check_mkdir(sjm_dir)
    logdir = './log'
    check_mkdir(logdir)

    sjm_file = f'{sjm_dir}/sjm.job'
    sjm_cmd = f'log_dir {logdir}\n'

    # parse_mapfile
    rna_fq_dict, atac_fq_dict = parse_mapfile(mapfile)

    with open(sjm_file, 'w') as fh:
        for sample in rna_fq_dict:
            fh.write(sjm_cmd + '\n')
            fh.write(f"""job_begin
    name rna_sample_{sample}
    sched_options -w n -cwd -V -l vf=1g,p=1
    cmd source activate celescope1.15.0; celescope rna sample --outdir ./{sample}/rna_dir/00.sample --sample {sample} --thread 4 --chemistry customized --fq1 {rna_fq_dict[sample][0]}
job_end

job_begin
    name rna_barcode_{sample}
    sched_options -w n -cwd -V -l vf=5g,p=1
    cmd source activate celescope1.15.0; celescope rna barcode --outdir ./{sample}/rna_dir/01.barcode --sample {sample} --thread 4 --chemistry customized --pattern C9L4C9L4C9U9 --whitelist /SGRNJ06/randd/USER/cjj/TESTDATA/test_rna/20220822_vx/whitelist --linker /SGRNJ06/randd/PROJECT/RD20073101_ScRNA_VDJ/shortlinker --lowNum 2 --allowNoLinker --fq1 {rna_fq_dict[sample][0]} --fq2 {rna_fq_dict[sample][1]}
job_end

""")

            fh.write(f"""job_begin
    name atac_sample_{sample}
    sched_options -w n -cwd -V -l vf=1g,p=1
    cmd source activate celescope_atac1.1.0; celescope atac sample --outdir ./{sample}/atac_dir/00.sample --sample {sample} --thread 8 --chemistry auto --wells 384 --fq2 {atac_fq_dict[sample][1]}
job_end

job_begin
    name atac_barcode_{sample}
    sched_options -w n -cwd -V -l vf=5g,p=1
    cmd source activate celescope_atac1.1.0; celescope atac barcode --outdir ./{sample}/atac_dir/01.barcode --sample {sample} --thread 8 --chemistry auto --lowNum 2 --allowNoLinker --wells 384 --fq1 {atac_fq_dict[sample][0]} --fq2 {atac_fq_dict[sample][1]} --fq3 {atac_fq_dict[sample][2]}
job_end

""")

            # 生成 Convert 和 ARC 作业内容
            fh.write(f"""job_begin
    name convert_{sample}
    sched_options -w n -cwd -V -l vf=5g,p=1
    cmd source activate cjj_atac; python /SGRNJ06/randd/USER/cjj/Script/Self-Used-Scripts/atac/convert.py --rna_path ./{sample}/rna_dir/01.barcode/ --atac_path ./{sample}/atac_dir/01.barcode/ --outdir ./{sample}/convert --sample {sample}
job_end

job_begin
    name arc_{sample}
    sched_options -w n -cwd -V -l vf=20g,p=4
    cmd /SGRNJ06/randd/PROJECT/scATAC_scRNA_multiOmics/scRNA/cell_ranger_test/multi_omics/cellranger-arc-2.0.2/cellranger-arc count --id={sample} --reference={reference} --libraries=./{sample}/cr_arc/libraries.csv  --localcores=4 --localmem=20 --min-atac-count 100 --min-gex-count 100
job_end

order atac_barcode_{sample} after atac_sample_{sample}
order rna_barcode_{sample} after rna_sample_{sample}
order convert_{sample} after rna_barcode_{sample}
order convert_{sample} after atac_barcode_{sample}
order arc_{sample} after convert_{sample}

""")

            # 生成 libraries.csv 文件
            check_mkdir(f"{sample}/cr_arc")
            libraries_csv_path = f"{sample}/cr_arc/libraries.csv"

            current_dir = os.getcwd()
            with open(libraries_csv_path, 'w') as libraries_file:
                libraries_file.write("fastqs,sample,library_type\n")
                libraries_file.write(f"{current_dir}/{sample}/convert/rna/,{sample},Gene Expression\n")
                libraries_file.write(f"{current_dir}/{sample}/convert/atac/,{sample},Chromatin Accessibility\n")

    print(f"SJM 作业文件已生成到: {sjm_dir}")
    print("样本对应的 libraries.csv 文件已生成！")


if __name__ == "__main__":
    # 使用 argparse 解析命令行参数
    parser = argparse.ArgumentParser(description="根据 mapfile 生成 SJM 和 libraries.csv 文件")
    parser.add_argument("--mapfile", required=True, help="mapfile")
    parser.add_argument("--reference", required=True, help="reference path")
    args = parser.parse_args()

    prepare_sjm(args.mapfile, args.reference)