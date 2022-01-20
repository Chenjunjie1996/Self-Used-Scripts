import pandas as pd
import sys
import subprocess
from concurrent.futures import ProcessPoolExecutor


def parse_mapfile(mapfile):
    mapfile = pd.read_csv(mapfile,sep = '\t',header=None)
    mapfile.columns = ['sample','path','merged_sample']
    mapfile['fq1'] = mapfile['path']+ '/'+ mapfile['sample'] + '*_R1.fastq.gz' 
    mapfile['fq2'] = mapfile['path']+ '/'+ mapfile['sample'] + '*_R2.fastq.gz'
    del mapfile['sample']
    del mapfile['path']
    fq1_dict = mapfile.groupby('merged_sample')['fq1'].apply(lambda x:x.tolist()).to_dict()
    fq2_dict = mapfile.groupby('merged_sample')['fq2'].apply(lambda x:x.tolist()).to_dict()
    return fq1_dict,fq2_dict

def get_cmd(fq1_dict,fq2_dict):
    cmd = []
    for key in fq1_dict:
        cmd1 = 'cat ' + ' '.join(fq1_dict[key]) + ' ' + '>' + ' ' + key + '_R1.fq.gz'
        cmd.append(cmd1)
    for key in fq2_dict:
        cmd2 = 'cat ' + ' '.join(fq2_dict[key]) + ' ' + '>' + ' ' + key + '_R2.fq.gz'
        cmd.append(cmd2)
    return cmd


def run_cmd(cmd):
    subprocess.check_call(cmd,shell=True)
    return cmd


def main():
    mapfile = sys.argv[1]
    fq1,fq2=parse_mapfile(mapfile)
    cmds = get_cmd(fq1,fq2)
    
    with ProcessPoolExecutor(max_workers=10) as executor:
        for result in executor.map(run_cmd, cmds):
            print(result,'Done')

if __name__ == '__main__':
    main()



