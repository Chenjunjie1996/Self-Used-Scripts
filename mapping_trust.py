import os
import subprocess
import glob
from concurrent.futures import ProcessPoolExecutor

TOOL_PATH = '/SGRNJ03/randd/cjj/Script/mapping.R'


def get_rds():
    rds_list = []
    contig_list = []
    sample_list = []
    assign_list = []

    shell_file = './run.sh'
    shell_file = os.path.abspath(shell_file)
    with open(shell_file, 'r') as f:
        lines = (line.strip() for line in f if line.strip().startswith('--mapfile'))
        for line in lines:
            mapfile_path = line.split()[1]
    if mapfile_path.startswith('..'):
        mapfile_path = '/'.join(shell_file.split('/')[:-2]) + '/' + mapfile_path.split('/')[-1]

    with open(mapfile_path, 'r') as f:
        lines = (line.strip() for line in f)
        for line in lines:
            trans_path = line.strip().split()[3]
            rds_file = glob.glob(f'{trans_path}/06.analysis/*.rds')[0]
            rds_list.append(rds_file)
            assign_file = glob.glob(f'{trans_path}/06.analysis/*_auto_assign/*_auto_cluster_type.tsv')[0]
            assign_list.append(assign_file)

            sample = line.split()[2]
            contig_file = glob.glob(f'./{sample}/04.summarize/{sample}_filtered_contig.csv')[0]
            contig_list.append(contig_file)
            sample_list.append(sample)

    return rds_list, contig_list, sample_list, assign_list


def run_mapping(rds, contig, sample, outdir, assign):
    cmd = (
        f'Rscript /SGRNJ03/randd/cjj/Script/mapping.R '
        f'--rds {rds} '
        f'--VDJ {contig} '
        f'--sample {sample} '
        f'--outdir {outdir} '
        f'--assign_file {assign}'
    )
    subprocess.check_call(cmd, shell=True)


if __name__ == '__main__':
    outdir = 'mapping_res'
    os.system(f'mkdir {outdir}')
    rds_list, contig_list, sample_list, assign_list = get_rds()
    outdir_list = [outdir] * len(rds_list)
    with ProcessPoolExecutor(max_workers=8) as executor:
        for result in executor.map(run_mapping, rds_list, contig_list, sample_list, outdir_list, assign_list):
            print(result, 'done')

