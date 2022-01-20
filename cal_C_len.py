import pandas as pd
from Bio.Seq import Seq
import sys
import os

def reversed_compl(seq):
    return str(Seq(seq).reverse_complement())

def parse(sample_file):
    sample_name = os.path.basename(sample_file)
    
    AIRR = pd.read_csv(f'{sample_file}/03.assemble/{sample_name}/outs/airr_rearrangement.tsv',sep = '\t')
    AIRR_C = AIRR[['cell_id','c_sequence_start','c_sequence_end']]
    AIRR_C.dropna(inplace=True)
    a_list = AIRR_C['c_sequence_start'].tolist()
    b_list = AIRR_C['c_sequence_end'].tolist()
    c_list = []
    
    for i in range(len(a_list)):
        c_list.append(b_list[i] - a_list[i])
    AIRR_C['C_length'] = c_list
    AVERAGE_LEN = sum(c_list)/ len(c_list)
    AVERAGE_LEN = round(AVERAGE_LEN,2)

    barcode_df = pd.read_csv(f'{sample_file}/02.convert/barcode_correspond.txt', sep='\t', index_col=1)
    barcode_dict = barcode_df.to_dict()['sgr']
    AIRR_C['cell_id'] = AIRR_C['cell_id'].apply(lambda x: reversed_compl(barcode_dict[x.split('-')[0]])+'_'+ x.split('-')[1])

    return sample_name, AVERAGE_LEN, AIRR_C, len(c_list)

if __name__ == '__main__':
    sample_file = sys.argv[1]
    sample_name, AVERAGE_LEN, AIRR_C, total_contigs = parse(sample_file)
    AIRR_C.to_csv(f'./{sample_name}_C_len.csv', sep='\t', index=False)
    with open(f'./{sample_name}_count.txt', 'w') as f:
        f.writelines("total_contigs : " + str(total_contigs) + "\n" + "average_C_len : " + str(AVERAGE_LEN))



