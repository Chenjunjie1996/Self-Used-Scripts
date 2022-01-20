import pandas as pd
from celescope.tools import utils
import pysam
from Bio.Seq import Seq
import glob
import sys
from concurrent.futures import ProcessPoolExecutor

def reversed_compl(seq):
    return str(Seq(seq).reverse_complement())

def parse_file(amplif):
    shell_file = '/'.join(amplif.split('/')[:-1])
    shell_file = f'{shell_file}/run.sh'
    sample_name = amplif.split('/')[-1]
    
    with open(shell_file,'r') as f:
        lines = (line.strip() for line in f if line.strip().startswith('--mapfile'))
        for line in lines:
            mapfile_path = line.split()[1]
            
    if mapfile_path.startswith('..'):
        mapfile_path = '/'.join(shell_file.split('/')[:-2]) + '/' + mapfile_path.split('/')[-1]
            
    with open(mapfile_path,'r') as f:
        lines = (line.strip() for line in f if sample_name in line)
        for line in lines:
            library_id = line.split()[0]
            library_path = line.split()[1]
            trans_path = line.split()[3]

    raw_reads = glob.glob(f'{library_path}/{library_id}*_R1.fastq*')[0]
    match_cell_barcodes, _match_cell_number = utils.read_barcode_file(trans_path)
    match_cell_barcodes_rev = set([reversed_compl(i) for i in match_cell_barcodes])

    read_num = 0
    total_read_num = 0
    with pysam.FastxFile(raw_reads) as fq:
        for read in fq:
            total_read_num += 1
            seq = read.sequence
            barcode = seq[9:17] + seq[33:41] + seq[57:65]
            if barcode in match_cell_barcodes_rev:
                read_num += 1
    read_percent = read_num / total_read_num
    read_percent = str(round(read_percent, 4) * 100) + '%'
    read_list = [sample_name, read_num, total_read_num, read_percent]

    return read_list

def write_to_file(results): 
    """
    res_file = open('./read_count.txt','w')
    for result in results:
        for _res in result:
            res_file.write(str(_res) + '\n')
        res_file.write('\n')
    res_file.close()
    """
    p_list = []
    row0 = ['sample_name', 'barcode_reads_count', 'total_reads', 'percent'] * len(results)
    row1 = []
    for result in results:
        for _res in result:
            row1.append(_res)
    for _row in range(len(row0)):
        p_list.append({'item':row0[_row], 'value':row1[_row]})
    p_df = pd.DataFrame(p_list)
    p_df.to_csv('./read_count.txt', sep=':', index=False, header=None)
    
    
if __name__ == '__main__':
    sample_list = []
    sample_list.append(sys.argv[1])
    try:
        sample_list.append(sys.argv[2])
    except:
        pass
        
    Res = []
    with ProcessPoolExecutor(max_workers=len(sample_list)) as executor:
        for res in executor.map(parse_file, sample_list):
            Res.append(res)
    write_to_file(Res)

