import glob
import os
import pandas as pd 
import numpy as np 


def test():
    """Detect whether there are empty value in each area of the filter annotations file.
    """
    out_file = open("test_result.txt", 'w')
    annotation_files = glob.glob("../*/*/03.assemble/*/outs/filtered_contig_annotations.csv")
    
    for file in annotation_files:
        df = pd.read_csv(file)
        df = df[["fwr1", "fwr1_nt", "cdr1", "cdr1_nt", "fwr2", "fwr2_nt", "cdr2", "cdr2_nt", "fwr3", "fwr3_nt", "cdr3", "cdr3_nt", "fwr4", "fwr4_nt"]]
        if np.any(df.isnull()):
            out_file.write(f"{os.path.abspath(file)}\n")
    
    out_file.close()


if __name__ == '__main__':
    test()