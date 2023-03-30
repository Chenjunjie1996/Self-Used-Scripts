import pandas as pd
import json
import pysam
import argparse
import glob
import os
import io
import numpy as np
from collections import defaultdict, OrderedDict
from celescope.tools import utils
from celescope.tools.emptydrop_cr import get_plot_elements
from celescope.tools.plotly_plot import Bar_plot
from jinja2 import Environment, FileSystemLoader, select_autoescape


def fasta_line(name, seq):
    return f'>{name}\n{seq}\n'


def check_mkdir(dir_name):
    """if dir_name is not exist, make one"""
    if not os.path.exists(dir_name):
        os.system(f"mkdir -p {dir_name}")

env = Environment(
    loader=FileSystemLoader("/SGRNJ06/randd/USER/cjj/Celescope/refactor_trust/CeleScope/celescope/templates/"),
    autoescape=select_autoescape(['html', 'xml']),
)


def gen_vj_annotation_metrics(df, seqtype):
    """
    Generate vdj Annotation Metrics from contig annotations file.
    """

    def get_vj_spanning_pair():
        """
        Get Productive V-J Spanning_Pair metric from annotation file
        Return productive chain pair number. eg: TRA/TRB or IGH/IGL, IGH/IGK.
        """
        df_productive = df[df['productive'] == True]

        if seqtype == "BCR":
            df_chain_heavy = df_productive[(df_productive['chain'] == 'IGH')]
            df_chain_light = df_productive[(df_productive['chain'] == 'IGL') | (df_productive['chain'] == 'IGK')]
        else:
            df_chain_heavy = df_productive[df_productive['chain'] == 'TRA']
            df_chain_light = df_productive[df_productive['chain'] == 'TRB']

        for _df in [df_chain_heavy, df_chain_light]:
            _df.drop_duplicates(['barcode'], inplace=True)

        VJ_Spanning_Pair_Cells = pd.merge(df_chain_heavy, df_chain_light, on='barcode', how='inner')

        return VJ_Spanning_Pair_Cells.shape[0]

    metric_dict = OrderedDict()
    chains, chain_pairs = ['TRA', 'TRB'], ['TRA_TRB']

    metric_dict["Cells With Productive V-J Spanning Pair"] = get_vj_spanning_pair()

    for pair in chain_pairs:
        chain1, chain2 = pair.split('_')[0], pair.split('_')[1]
        cbs1 = set(df[(df['productive'] == True) & (df['chain'] == chain1)].barcode)
        cbs2 = set(df[(df['productive'] == True) & (df['chain'] == chain2)].barcode)
        metric_dict[f"Cells With Productive V-J Spanning ({chain1}, {chain2}) Pair"] = len(cbs1.intersection(cbs2))

    for chain in chains:
        metric_dict[f"Cells With {chain} Contig"] = len(set(df[df['chain'] == chain].barcode))
        metric_dict[f"Cells With CDR3-annotated {chain} Contig"] = len(
            set(df[(df['chain'] == chain) & (df['cdr3'] != 'None')].barcode))
        metric_dict[f"Cells With V-J Spanning {chain} Contig"] = len(
            set(df[(df['full_length'] == True) & (df['chain'] == chain)].barcode))
        metric_dict[f"Cells With Productive {chain} Contig"] = len(
            set(df[(df['productive'] == True) & (df['chain'] == chain)].barcode))

    return metric_dict


def gen_clonotypes_table(df, out_clonotypes):
    """
    Generate clonotypes.csv file
    """
    df = df[df['productive'] == True]
    df['chain_cdr3aa'] = df[['chain', 'cdr3']].apply(':'.join, axis=1)

    match_clonotypes = open(out_clonotypes, 'w')
    match_clonotypes.write('barcode\tcdr3s_aa\n')
    for cb in set(df.barcode):
        temp = df[df['barcode']==cb].sort_values(by='chain', ascending=True)
        chain_pair = ';'.join(temp['chain_cdr3aa'].tolist())
        match_clonotypes.write(f'{cb}\t{chain_pair}\n')
    match_clonotypes.close()

    df_clonetypes = pd.read_csv(out_clonotypes, sep='\t', index_col=None)
    df_clonetypes = df_clonetypes.groupby('cdr3s_aa', as_index=False).agg({'barcode': 'count'})
    df_clonetypes.rename(columns={'barcode': 'frequency'}, inplace=True)
    df_clonetypes = df_clonetypes.sort_values("frequency", ascending=False)
    sum_f = df_clonetypes['frequency'].sum()
    df_clonetypes['proportion'] = df_clonetypes['frequency'].apply(lambda x: x/sum_f)
    df_clonetypes['clonotype_id'] = [f'clonotype{i}' for i in range(1, df_clonetypes.shape[0]+1)]
    df_clonetypes = df_clonetypes.reindex(columns=['clonotype_id', 'cdr3s_aa', 'frequency', 'proportion'])
    df_clonetypes.sort_values(by='frequency', ascending=False, inplace=True)
    df_clonetypes.to_csv(out_clonotypes, sep=',', index=False)


class VDJ_calling:
    """
    cell calling for cellranger result.
    对 is_cell=False的细胞 calling回组装出多条高可信度的productive链的barcode不被判定为细胞
    条件:
    1. 满足productive=True
    2. high_confidence=True
    3. groupby barcode chain，取最高umi的链
    4. calling回的细胞合并至filtered annotation文件
    """

    def __init__(self, args):

        # in
        self.all_contig_anno = pd.read_csv(glob.glob(f"{args.cr_dir}/03.assemble/*/outs/all_contig_annotations.csv")[0])
        self.filter_contig_anno = pd.read_csv(glob.glob(f"{args.cr_dir}/03.assemble/*/outs/filtered_contig_annotations.csv")[0])
        self.all_contig_fasta = glob.glob(f"{args.cr_dir}/03.assemble/*/outs/all_contig.fasta")[0]
        self.filter_contig_fasta = glob.glob(f"{args.cr_dir}/03.assemble/*/outs/filtered_contig.fasta")[0]
        self.all_bam = glob.glob(f"{args.cr_dir}/03.assemble/*/outs/all_contig.bam")[0]
        self.match_dir = args.match_dir
        self.match_cell_barcodes, _ = utils.get_barcode_from_match_dir(self.match_dir)
        with open(f"{args.cr_dir}/02.convert/barcode_convert.json", 'r') as f:
            self.tenX_sgr = json.load(f)

        # out
        self.out_dir = f"{args.cr_dir}/07.cell_calling"
        check_mkdir(self.out_dir)

    def __call__(self):
        calling_cells, df_merge, df_match = self.cell_calling()
        self.Barcode_rank_plot(calling_cells)
        self.render_html(calling_cells, df_merge, df_match)

    def get_table_dict(self, title, table_id, df_table):
        """
        table_dict {title: '', table_id: '', df_table: pd.DataFrame}
        """
        table_dict = {}
        table_dict['title'] = title
        table_dict['table'] = df_table.to_html(
            escape=False,
            index=False,
            table_id=table_id,
            justify="center")
        table_dict['id'] = table_id
        return table_dict

    @utils.add_log
    def cell_calling(self):
        df_productive = self.all_contig_anno[self.all_contig_anno["productive"] == True]
        df_call = df_productive[df_productive["is_cell"] == False]
        df_call = df_call[df_call["high_confidence"] == True]
        df_call = df_call.sort_values("umis", ascending=False)
        df_call = df_call.groupby(["barcode", "chain"], as_index=False).head(1)
        df_call = df_call.groupby("barcode").filter(lambda x: (len(x) > 1))
        df_merge = pd.concat([self.filter_contig_anno, df_call])
        calling_cells = set(df_merge.barcode)
        calling_contigs = set(df_merge.contig_id)

        # generate fasta file
        out_filter_fasta = open(f"{self.out_dir}/filtered_contig.fasta", 'w')
        out_match_filter_fasta = open(f"{self.out_dir}/matched_contig.fasta", 'w')
        with pysam.FastxFile(self.all_contig_fasta) as f:
            for entry in f:
                name = entry.name
                seq = entry.sequence
                attrs = name.split('_')
                if name in calling_contigs:
                    convert_bc = self.tenX_sgr[attrs[0].split('-')[0]]
                    new_name = convert_bc + '_' + attrs[1] + '_' + attrs[2]
                    out_filter_fasta.write(fasta_line(new_name, seq))
                    if convert_bc in self.match_cell_barcodes:
                        out_match_filter_fasta.write(fasta_line(new_name, seq))

        out_filter_fasta.close()
        out_match_filter_fasta.close()

        # generate annotation file
        df_merge['barcode'] = df_merge['barcode'].apply(lambda x: self.tenX_sgr[x.split('-')[0]])
        df_merge['contig_id'] = df_merge['contig_id'].apply(
            lambda x: self.tenX_sgr[x.split('-')[0]] + '_' + x.split('_')[1] + '_' + x.split('_')[2])
        df_merge['is_cell'] = "TRUE"
        df_merge.to_csv(f"{self.out_dir}/filtered_contig_annotations.csv", sep=',', index=False)
        df_match = df_merge[df_merge["barcode"].isin(self.match_cell_barcodes)]
        df_match.to_csv(f"{self.out_dir}/matched_contig_annotations.csv", sep=',', index=False)

        return calling_cells, df_merge, df_match

    @utils.add_log
    def Barcode_rank_plot(self, calling_cells):
        dic_umi = defaultdict(set)

        with pysam.AlignmentFile(self.all_bam) as fh:
            for read in fh:
                cb = read.get_tag('CB')
                umi = read.get_tag('UB')
                dic_umi[cb].add(umi)

        df_umi = pd.DataFrame()
        df_umi['barcode'] = list(dic_umi.keys())
        df_umi['UMI'] = [len(dic_umi[i]) for i in dic_umi]
        df_umi = df_umi.sort_values(by='UMI', ascending=False)
        df_umi['mark'] = df_umi['barcode'].apply(lambda x: 'CB' if x in calling_cells else 'UB')
        df_umi['barcode'] = df_umi['barcode'].apply(lambda x: self.tenX_sgr[x.split('-')[0]])

        df_umi.to_csv(f"{self.out_dir}/count.txt", sep='\t', index=False)
        self.add_data(chart=get_plot_elements.plot_barcode_rank(self.count_file))

    @utils.add_log
    def render_html(self, calling_cells, df_merge, df_match):
        os.system(f"cp {self.out_dir}/../.data.json {self.out_dir}")
        fh = open(f"{self.out_dir}/.data.json")
        data = json.load(fh)
        fh.close()

        """
        Estimated Number of Cells	6,983
        Fraction Reads in Cells	28.2%
        Mean Reads per Cell	5,756
        Mean Used Reads per Cell	1,518
        Median Used TRA UMIs per Cell	4
        Median Used TRB UMIs per Cell
        """
        cells_reads, total_reads = 0, 0
        with pysam.AlignmentFile(self.all_bam) as f:
            for read in f:
                total_reads += 1
                cb = read.get_tag('CB')
                if cb in calling_cells:
                    cells_reads += 1

        data["cells_summary"]["metric_list"][0]["display"] = str(format(len(calling_cells), ','))
        data["cells_summary"]["metric_list"][1]["display"] = f'{round(cells_reads / total_reads * 100, 2)}%'
        data["cells_summary"]["metric_list"][2]["display"] = str(format(total_reads // len(calling_cells), ','))
        data["cells_summary"]["metric_list"][3]["display"] = str(format(cells_reads // len(calling_cells), ','))
        data["cells_summary"]["metric_list"][4]["display"] = np.median(df_merge[df_merge["chain"]=="TRA"].umis)
        data["cells_summary"]["metric_list"][5]["display"] = np.median(df_merge[df_merge["chain"]=="TRB"].umis)
        data["cells_summary"]["chart"] = get_plot_elements.plot_barcode_rank(f"{self.out_dir}/count.txt")

        """
        Annotation 
        Cells With Productive V-J Spanning Pair	49.5%
        Cells With Productive V-J Spanning (TRA, TRB) Pair	49.5%
        Cells With TRA Contig	74.1%
        Cells With CDR3-annotated TRA Contig	59.6%
        Cells With V-J Spanning TRA Contig	63.7%
        Cells With Productive TRA Contig	50.5%
        Cells With TRB Contig	99.3%
        Cells With CDR3-annotated TRB Contig	99.1%
        Cells With V-J Spanning TRB Contig	99.1%
        Cells With Productive TRB Contig	99.0%
        """
        metrics_dict = gen_vj_annotation_metrics(df_merge, seqtype="TCR")
        for k, v in metrics_dict.items():
            for i in data["annotation_summary"]["metric_list"]:
                if i["name"] == k:
                    i["display"] = f'{round(v / len(calling_cells) * 100, 2)}%'

        data["match_summary"]["metric_list"][0]["display"] = str(format(len(set(df_match.barcode)), ','))
        metrics_dict = gen_vj_annotation_metrics(df_match, seqtype="TCR")
        for k, v in metrics_dict.items():
            for i in data["match_summary"]["metric_list"]:
                if i["name"] == k:
                    i["display"] = f'{round(v / len(set(df_match.barcode)) * 100, 2)}%'

        gen_clonotypes_table(df_match, f"{self.out_dir}/matched_clonotypes.csv")
        title = 'Clonetypes'
        raw_clonotypes = pd.read_csv(f"{self.out_dir}/matched_clonotypes.csv", sep=',', index_col=None)
        raw_clonotypes['ClonotypeID'] = raw_clonotypes['clonotype_id'].apply(lambda x: x.strip('clonetype'))
        raw_clonotypes['Frequency'] = raw_clonotypes['frequency']
        raw_clonotypes['Proportion'] = raw_clonotypes['proportion'].apply(lambda x: f'{round(x*100, 2)}%')
        raw_clonotypes['CDR3_aa'] = raw_clonotypes['cdr3s_aa'].apply(lambda x: x.replace(';', '<br>'))

        table_dict = self.get_table_dict(
            title=title,
            table_id='clonetypes',
            df_table=raw_clonotypes[['ClonotypeID', 'CDR3_aa', 'Frequency', 'Proportion']]
        )
        data["match_summary"]["table_dict"] = table_dict

        raw_clonotypes['ClonotypeID'] = raw_clonotypes['ClonotypeID'].astype("int")
        raw_clonotypes.sort_values(by=['ClonotypeID'], inplace=True)
        Barplot = Bar_plot(df_bar=raw_clonotypes).get_plotly_div()
        data["match_summary"]["Barplot"] = Barplot

        template = env.get_template(f'html/flv_CR/base.html')
        report_html = f"{self.out_dir}/cellcalling_report.html"
        with io.open(report_html, 'w', encoding='utf8') as fh:
            html = template.render(data)
            # fh.write(html.encode('utf-8'))
            fh.write(html)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="cell calling for flv_CR result")
    parser.add_argument("--cr_dir", help="cr dir", required=True)
    parser.add_argument("--match_dir", help="match dir", required=True)
    args = parser.parse_args()
    VDJ_calling(args)()