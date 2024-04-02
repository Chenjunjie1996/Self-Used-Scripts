import os
import pandas as pd
import glob
import re


def deal_matrix_file(raw_path, target_path):
    os.system(f"mkdir -p {target_path}")
    os.system(f"cp {raw_path} {target_path}")
    os.system(f"tar -xvf {target_path}/*_matrix_10X.tar -C {target_path}")


def deal_analysis_file(raw_path, target_path):
    os.system(f"mkdir -p {target_path}")
    os.system(f"cp {raw_path}/* {target_path}")


def parse_mapfile(mapfile):
    mapfile = pd.read_csv(mapfile, sep='\t', names=["id", "path", "sample", "lims_path"])
    lims_path_list = mapfile['lims_path'].tolist()
    new_lims_path_list = []
    
    for i in lims_path_list:
        tar_path = glob.glob(f"{i}/*/*/*_matrix_10X.tar")[0]
        matrix = tar_path.split('/')[-1]
        patt = r'_matrix_10X.tar'
        sample = re.sub(patt, "", matrix)
        
        target_matrix_path = os.path.abspath(f"./match_dir/{sample}/05.count/{sample}_matrix_10X")
        deal_matrix_file(tar_path, target_matrix_path)
        
        analysis_path = glob.glob(f"{i}/*/*/06.analysis*")[0]
        target_analysis_path = os.path.abspath(f"./match_dir/{sample}/06.analysis")
        deal_analysis_file(analysis_path, target_analysis_path)
        
        new_lims_path_list.append('/'.join(target_matrix_path.split('/')[:-2]))

    mapfile['lims_path'] = new_lims_path_list

    return mapfile


def main():
    mapfile = glob.glob('./*.mapfile')[0]
    new_mapfile = parse_mapfile(mapfile)
    new_mapfile.to_csv("./new.mapfile", sep='\t', header=False, index=False)


if __name__ == '__main__':
    main()

