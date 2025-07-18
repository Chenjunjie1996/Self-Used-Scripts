import glob
import pandas as pd


INDEX_DICT = {
    "CTCCAT": "TS-B1",
    "ATCCTC": "TS-B2",
    "ACGTCT": "TS-B3",
    "TGCGAA": "TS-B4",
    "TTCTCG": "TS-B5",
    "CTGCTA": "TS-B6",
    "CATGAT": "TS-B7",
    "TCAACT": "TS-B8",
    "AGTCCT": "TS-B9",
    "GTTGAG": "TS-B10",
    "TAGCTG": "TS-B11",
    "TCGCCA": "TS-B12",
    "GAACTC": "TS-B13",
    "TATGGT": "TS-B14",
    "CGCAAC": "TS-B15",
    "TGGCAG": "TS-B16",
    "ATGCAT": "TS-B17",
    "GACTAT": "TS-B18",
    "GTGATT": "TS-B19",
    "CTCTTG": "TS-B20",
    "AAGCGT": "TS-B21",
    "CCAACA": "TS-B22",
    "GTTGGT": "TS-B23",
    "TCTAGT": "TS-B24",
    "TATGTG": "TS-B25",
    "TGTGGC": "TS-B26",
    "GACCTG": "TS-B27",
    "TTCCGT": "TS-B28",
    "AAGGCA": "TS-B29",
    "TAGGAT": "TS-B30",
    "AACTCC": "TS-B31",
    "TCGATG": "TS-B32",
    "CTGCGT": "TS-B33",
    "GTCGGA": "TS-B34",
    "TCACAT": "TS-B35",
    "ATAGGT": "TS-B36",
    "CGTAAT": "TS-B37",
    "GCATGT": "TS-B38",
    "AACTGA": "TS-B39",
    "GCACAA": "TS-B40",
    "GCCATC": "TS-B41",
    "CAACCG": "TS-B42",
    "GTCTGG": "TS-B43",
    "TCCATT": "TS-B44",
    "CAGACC": "TS-B45",
    "ACGGAG": "TS-B46",
    "ACATCA": "TS-B47",
    "TATCCG": "TS-B48",
    "GGAGAG": "TS-B49",
    "CCAATG": "TS-B50",
    "TTCTGA": "TS-B51",
    "GTGACG": "TS-B52",
    "ATGGTG": "TS-B53",
    "ACTTGT": "TS-B54",
    "ATAGAC": "TS-B55",
    "CCTATA": "TS-B56",
    "TTAAGG": "TS-B57",
    "GATCAC": "TS-B58",
    "TAGCCT": "TS-B59",
    "AGCGCT": "TS-B60",
    "AGACGC": "TS-B61",
    "CTAAGA": "TS-B62",
    "TATCGA": "TS-B63",
    "CGCACA": "TS-B64",
    "CAAGTT": "TS-B65",
    "GAACCA": "TS-B66",
    "TACACA": "TS-B67",
    "CATTGG": "TS-B68",
    "TCATGC": "TS-B69",
    "AGGTTA": "TS-B70",
    "TCGAAT": "TS-B71",
    "TCTTGG": "TS-B72",
    "CTCTAC": "TS-B73",
    "GAGGTC": "TS-B74",
    "ACAACG": "TS-B75",
    "CAGATA": "TS-B76",
    "CAGGTA": "TS-B77",
    "TCTTAC": "TS-B78",
    "CCTGTG": "TS-B79",
    "TCGAGC": "TS-B80",
    "CTGAAT": "TS-B81",
    "ATTGGC": "TS-B82",
    "CATCTT": "TS-B83",
    "TCTCTA": "TS-B84",
    "GCGTCA": "TS-B85",
    "GTTCAT": "TS-B86",
    "AATCAG": "TS-B87",
    "CGGTGT": "TS-B88",
    "TCCGTC": "TS-B89",
    "CTCACC": "TS-B90",
    "TTGACT": "TS-B91",
    "GCCGTA": "TS-B92",
    "CGACTC": "TS-B93",
    "ATCCAA": "TS-B94",
    "TGCCAT": "TS-B95",
    "ACGATA": "TS-B96",
}


def main():
    
    df = glob.glob('*/*/*_filtered_annotations.csv')[0]
    sample = df.split('/')[0]
    df = pd.read_csv(df)
    barcode_sample = pd.read_csv('barcode_sample.xls',sep='\t')
    barcode_sample_dict = barcode_sample.set_index('barcode')['sample'].to_dict()
    
    new_index_dict = {v: k for k, v in INDEX_DICT.items()}
    df['barcode'] = df["barcode"].apply(lambda x: new_index_dict[x])
    df['barcode'] = df["barcode"].apply(lambda x: barcode_sample_dict[x])
    df.to_csv(f'./{sample}_sample_filtered_annotations.csv', sep=",", index=False)


if __name__ == '__main__':
    main()