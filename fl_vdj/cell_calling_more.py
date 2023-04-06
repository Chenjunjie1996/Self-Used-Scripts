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


CHAIN = {
	'TCR': ['TRA', 'TRB'], 
	'BCR': ['IGH', 'IGL', 'IGK']
	}
PAIRED_CHAIN = {
	'TCR': ['TRA_TRB'],
	'BCR': ['IGK_IGH', 'IGL_IGH']
}


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
    chains, chain_pairs = CHAIN[seqtype], PAIRED_CHAIN[seqtype]

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


class Filter_noise:
    """ Filter noise barcode.
    auto:
    keep paired-chain barcode when highest umi of contig >= 2 * second highest umi of contig.
    
    snr:
    https://www.nature.com/articles/s41598-017-18303-z#:~:text=The%20signal%2Dto%2Dnoise%20ratio,of%20background%20pixel%20intensity%2C%20respectively.
    keep paired-chain barcode when SNR value >= coefficient.
    The signal-to-noise ratio (SNR) of the single molecules was defined as SNR = (S − B)/σ,
    where S is the peak single molecule pixel intensity and B and 
    σ are the average and standard deviation of background pixel intensity, respectively.
    
    not_filter:
    for paired-chain barcode, keep highest tra and trb chain.
    """
    
    def __init__(self, args, df):
        self.method = args.method
        self.coeff = float(args.coeff)
        self.df = df
        self.seqtype = args.seqtype
        self.df = self.df.sort_values("umis", ascending=False)
    
    @utils.add_log
    def __call__(self):
        
        if self.method:
            bc_chain_dict = self.df.groupby("barcode")["chain"].apply(lambda x: set(x)).to_dict()
            if self.seqtype == "TCR":
                pair_chain_dict = {key: value for key, value in bc_chain_dict.items() if len(value)==2}
            else:
                pair_chain_dict = {key: value for key, value in bc_chain_dict.items() if len(value)>=2 and "IGH" in value}
            self.df = self.df[self.df["barcode"].isin(pair_chain_dict)]
            self.df = self.df.sort_values(["barcode","umis"], ascending=[False, False])
            
            # 有多条重链或轻链的barcode
            df_multi_chain = self.df.groupby("barcode").filter(lambda x: (len(x) > 2))
            # 仅有一对轻重链配对的barcode
            df_pair_chain = self.df.groupby("barcode").filter(lambda x: (len(x) == 2))
            
            if self.seqtype == "TCR":
                df_heavy = df_multi_chain[df_multi_chain['chain'] == 'TRB']
                df_light = df_multi_chain[df_multi_chain['chain'] == 'TRA']
            else:
                df_heavy = df_multi_chain[(df_multi_chain['chain'] == 'IGH')] 
                df_light = df_multi_chain[(df_multi_chain['chain'] == 'IGL') | (df_multi_chain['chain'] =='IGK')]
            light_dict = df_light.groupby("barcode")["umis"].apply(lambda x: x.tolist()).to_dict()
            heavy_dict = df_heavy.groupby("barcode")["umis"].apply(lambda x: x.tolist()).to_dict()

            if self.method == "auto":
                for chain_dict in [light_dict, heavy_dict]:
                    for k in chain_dict:
                        if len(chain_dict[k]) == 1:
                            chain_dict[k].append(0.1)

                # highest umi of contig >= coeff * second highest umi of contig
                filter_noise_light = {key for key,value in light_dict.items() if value[0]/value[1] >= self.coeff}
                filter_noise_heavy = {key for key,value in heavy_dict.items() if value[0]/value[1] >= self.coeff}
                
            elif self.method == "snr":
                for chain_dict in [light_dict, heavy_dict]:
                    for k in chain_dict:
                        if len(chain_dict[k]) == 1:
                            chain_dict[k].append(chain_dict[k][0])
                
                filter_noise_light = self.snr_filter(light_dict, self.coeff)
                filter_noise_heavy = self.snr_filter(heavy_dict, self.coeff)
                        
            filter_noise_barcode = filter_noise_light & filter_noise_heavy | set(df_pair_chain.barcode)
            self.df = self.df[self.df["barcode"].isin(filter_noise_barcode)]
            
        else:
            self.df = self.not_filter(self.df, self.seqtype)
        
        return self.df

    @staticmethod
    def not_filter(df, seqtype):
        """do not filter any noise, keep all paired-chain barcode"""
        
        if seqtype == "TCR":
            df = df.groupby(["barcode", "chain"], as_index=False).head(1)
            df = df.groupby("barcode").filter(lambda x: (len(x) > 1))
        else:
            df_chain_heavy = df[(df['chain'] == 'IGH')]
            df_chain_light = df[(df['chain'] == 'IGL') | (df['chain'] =='IGK')]
            df_chain_heavy.drop_duplicates(['barcode'], inplace=True)
            df_chain_light.drop_duplicates(['barcode'], inplace=True)
            df = pd.concat([df_chain_heavy, df_chain_light])
        
        return df
    
    @staticmethod
    def snr_filter(chain_dict, coeff):
        """ calculate SNR for each barcode
        SNR = (S − B)/σ
        
        :param chain_dict: {'AAACATCGAAGGACACCCTAATCC': [7, 3, 1],
                         'AAACATCGAATCCGTCCATCAAGT': [10, 2, 1],
                         'AAACATCGAATCCGTCCCGACAAC': [4, 3],
                         }
        :return: {AAACATCGAAGGACACCCTAATCC, AAACATCGAATCCGTCCATCAAGT}
        """
        filter_noise_barcode = set()
        for k, v in chain_dict.items():
            S, noise = v[0], v[1:]
            B = np.mean(noise)
            O = np.std(noise)
            if O == 0:
                filter_noise_barcode.add(k)
                continue
            SNR = (S - B) / O
            if SNR >= coeff:
                filter_noise_barcode.add(k)

        return filter_noise_barcode


class VDJ_calling:
    """
    cell calling for flv_CR result.
    对is_cell=False的细胞 calling回组装出多条productive链而不被判定为细胞的barcode。
    条件:
        1.满足productive=True, is_cell=False
        2.对这样的barcode中的contig进行过滤：SNR, AUTO, NOT_FILTER
        3.仅call回双链细胞，calling回的细胞合并至filtered annotation文件
        4.生成新的结果文件和报告 outdir: 07.cell_calling
    """

    def __init__(self, args):

        # in
        self.all_contig_anno = pd.read_csv(glob.glob(f"{args.cr_dir}/03.assemble/*/outs/all_contig_annotations.csv")[0])
        self.filter_contig_anno = pd.read_csv(glob.glob(f"{args.cr_dir}/03.assemble/*/outs/filtered_contig_annotations.csv")[0])
        self.all_contig_fasta = glob.glob(f"{args.cr_dir}/03.assemble/*/outs/all_contig.fasta")[0]
        self.filter_contig_fasta = glob.glob(f"{args.cr_dir}/03.assemble/*/outs/filtered_contig.fasta")[0]
        self.all_bam = glob.glob(f"{args.cr_dir}/03.assemble/*/outs/all_contig.bam")[0]
        self.match_dir = args.match_dir
        self.seqtype = args.seqtype
        
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
        df_call = Filter_noise(args, df_call)()
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
        # self.add_data(chart=get_plot_elements.plot_barcode_rank(self.count_file))

    @utils.add_log
    def render_html(self, calling_cells, df_merge, df_match):
        os.system(f"cp {self.out_dir}/../.data.json {self.out_dir}")
        fh = open(f"{self.out_dir}/.data.json")
        data = json.load(fh)
        fh.close()

        """
        Cells
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
        
        if self.seqtype == "TCR":
            data["cells_summary"]["metric_list"][4]["display"] = np.median(df_merge[df_merge["chain"]=="TRA"].umis)
            data["cells_summary"]["metric_list"][5]["display"] = np.median(df_merge[df_merge["chain"]=="TRB"].umis)
        else:
            data["cells_summary"]["metric_list"][4]["display"] = np.median(df_merge[df_merge["chain"]=="IGH"].umis)
            data["cells_summary"]["metric_list"][5]["display"] = np.median(df_merge[df_merge["chain"]=="IGL"].umis)
            data["cells_summary"]["metric_list"][6]["display"] = np.median(df_merge[df_merge["chain"]=="IGK"].umis)
            
        data["cells_summary"]["chart"] = get_plot_elements.plot_barcode_rank(f"{self.out_dir}/count.txt",  log_uniform=True)


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
        metrics_dict = gen_vj_annotation_metrics(df_merge, self.seqtype)
        for k, v in metrics_dict.items():
            for i in data["annotation_summary"]["metric_list"]:
                if i["name"] == k:
                    i["display"] = f'{round(v / len(calling_cells) * 100, 2)}%'

        """
        Match 
        Cells Match with ScRNA-seq Analysis	1,128
        Cells With Productive V-J Spanning Pair	727(64.45%)
        Cells With Productive V-J Spanning (TRA, TRB) Pair	727(64.45%)
        Cells With TRA Contig	784(69.5%)
        Cells With CDR3-annotated TRA Contig	784(69.5%)
        Cells With V-J Spanning TRA Contig	784(69.5%)
        Cells With Productive TRA Contig	784(69.5%)
        Cells With TRB Contig	1,071(94.95%)
        Cells With CDR3-annotated TRB Contig	1,071(94.95%)
        Cells With V-J Spanning TRB Contig	1,071(94.95%)
        Cells With Productive TRB Contig	1,071(94.95%)
        """
        data["match_summary"]["metric_list"][0]["display"] = str(format(len(set(df_match.barcode)), ','))
        metrics_dict = gen_vj_annotation_metrics(df_match, self.seqtype)
        for k, v in metrics_dict.items():
            for i in data["match_summary"]["metric_list"]:
                if i["name"] == k:
                    i["display"] = f'{round(v / len(set(df_match.barcode)) * 100, 2)}%'

        """
        Clonotypes table
        """
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
    parser.add_argument("--method", help="filter method", default=None)
    parser.add_argument("--coeff",
                        help="coefficient will affect auto and snr noise filter, recommend 2 for auto, 10 for snr",
                        default=2
                        )
    parser.add_argument("--seqtype", choices=['TCR', 'BCR'], help="TCR or BCR", required=True)
    args = parser.parse_args()
    VDJ_calling(args)()