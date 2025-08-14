import pandas as pd
import glob
import json


def dump_dict_to_json(d, json_file):
    with open(json_file, 'w') as f:
        json.dump(d, f, indent=4)


def main():
    files = glob.glob('*/outs/fragments_corrected_dedup_count.tsv.gz')
    samples = [i.split('/')[0] for i in files]
    res = dict()

    for i in range(len(files)):
        df = pd.read_csv(files[i], header=None, sep='\t', names=["chrom", "chromStart", "chromEnd", "barcode", "count"])
        df_mt = df[df['chrom'].str.contains('chrM')]
        mt_percent = sum(df_mt['count']) / sum(df['count'])
        res[samples[i]] = mt_percent

    dump_dict_to_json(res, 'mt_percent.json')


if __name__ == '__main__':
    main()

    