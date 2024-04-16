import glob
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor


def check_mkdir(dir_name):
    """if dir_name is not exist, make one"""
    if not os.path.exists(dir_name):
        os.system(f"mkdir -p {dir_name}")


def run_tss(species, count, meta, fragment, sample, outdir):
    cmd = (
        f'Rscript /SGRNJ06/randd/USER/cjj/Script/Self-Used-Scripts/atac/tssplot.R '
        f'--species {species} '
        f'--count {count} '
        f'--meta {meta} '
        f'--fragment {fragment} '
        f'--sample {sample} '
        f'--outdir {outdir} '
    )
    subprocess.check_call(cmd, shell=True)
        
        
def main():
    species = sys.argv[1]
    
    outdir = 'TssPlot'
    check_mkdir('outdir')
    
    dir_list = glob.glob('*/outs')
    sample_list = [i.split('/')[0] for i in dir_list]
    count_list = [f'{i}/outs/{i}_filtered_peak_count.h5' for i in sample_list]
    meta_list = [f'{i}/outs/cell_qc_metrics.tsv' for i in sample_list]
    fragment_list = [f'{i}/outs/fragments_corrected_dedup_count.tsv.gz' for i in sample_list]
    species = [species] * len(dir_list)
    outdir = [outdir] * len(dir_list)
        
    with ProcessPoolExecutor(max_workers=8) as executor:
        for result in executor.map(run_tss, species, count_list, meta_list, fragment_list, sample_list, outdir):
            print(result, 'done')


if __name__ == '__main__':
    main()