import glob
import os
import pandas as pd 


def test():
    """Del last nt in fwr4 region of filtered annotation files.
    """
    annotation_files = glob.glob("../*/*/04.summarize/filtered_contig_annotations.csv")
    
    for file in annotation_files:
        del_nt_file = os.path.abspath(file).replace("filtered_contig_annotations.csv", "del_nt_filtered_contig_annotations.csv")
        if not os.path.exists(del_nt_file):
            df = pd.read_csv(file)
            df["fwr4_nt"] = df["fwr4_nt"].apply(lambda x: x[:-1])
            df.to_csv(del_nt_file, sep=',', index=False)


if __name__ == '__main__':
    test()