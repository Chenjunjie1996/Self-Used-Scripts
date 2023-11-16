import pandas as pd 
import glob


def main():
    files = glob.glob("./*/*/matched_clonotypes.csv")
    for file in files:
        sample = file.split('/')[-3]
        out = f"./{sample}_match_pair_clonotypes.csv"
    
        df = pd.read_csv(file)
        df.cdr3s_aa = df.cdr3s_aa.apply(lambda x: x.split(';'))
        df = df[df['cdr3s_aa'].apply(lambda x: len(x)) == 2]
        df = df[df['cdr3s_aa'].apply(lambda x: x[0]).str.contains('IGH')]
        df = df[~df['cdr3s_aa'].apply(lambda x: x[1]).str.contains('IGH')]
        df['cdr3s_aa'] = df['cdr3s_aa'].apply(lambda x: ';'.join(x))
        df.to_csv(out, sep=',', index=False)


if __name__ == '__main__':
    main()