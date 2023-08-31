import csv
import sys
import re
import collections
import scanpy as sc
import argparse
import plotly.graph_objects as go
#pio.templates.default = 'plotly_white'


PATTERN = re.compile(r'(\S+?)\s*"(.*?)"')


class GeneIdNotFound(Exception):
    pass


def get_properties_dict(properties_str):
    """
    allow no space after semicolon
    """
        
    if isinstance(properties_str, dict):
        return properties_str

    properties = collections.OrderedDict()
    attrs = properties_str.split(';')
    for attr in attrs:
        if attr:
            m = re.search(PATTERN, attr)
            if m:
                key = m.group(1)
                key = key.strip()
                value = m.group(2)
                value = value.strip()
                properties[key] = value

    return properties


def gtf_reader_iter():
    with open("/SGRNJ06/randd/USER/cjj/filtered_ref/hs_ensembl_99/Homo_sapiens.GRCh38.99.gtf", mode='rt') as f:
        reader = csv.reader(f, delimiter='\t')
        for i, row in enumerate(reader, start=1):
            if len(row) == 0:
                continue
            if row[0].startswith('#'):
                yield row, True, None, None
                continue
            if len(row) != 9:
                sys.exit(f"Invalid number of columns in GTF line {i}: {row}\n")
            if row[6] not in ['+', '-']:
                sys.exit(f"Invalid strand in GTF line {i}: {row}\n")
            properties = get_properties_dict(row[8])
            annotation = row[2]
            if annotation == 'exon':
                if 'gene_id' not in properties:
                    raise GeneIdNotFound(f"Property 'gene_id' not found in GTF line {i}: {row}\n")
            yield row, False, annotation, properties


def get_lncRNA_gene_set(): 
    lncRNA_gene_set = set()
    
    for _row, _is_comment, annotation, properties in gtf_reader_iter():
        if _is_comment:
            continue
        properties = dict(properties)
        if properties["gene_biotype"] == "lncRNA":
            gene_name = properties["gene_name"]
            lncRNA_gene_set.add(gene_name)
    
    return lncRNA_gene_set


def box_plot(matrix_files):
    fig = go.Figure()

    for matrix_file in matrix_files:
        adata = sc.read_10x_mtx(matrix_file,var_names='gene_symbols',)
        sample = matrix_file.split('/')[-1]
        adata.var["lncRNA"] = adata.var_names.str.upper().isin(lncRNA_gene_set)

        sc.pp.calculate_qc_metrics(
            adata,
            qc_vars=["lncRNA"], 
            percent_top=None,
            use_raw=False,
            log1p=False, 
            inplace=True
        )
    
        df = adata.obs
        fig.add_trace(go.Box(y=list(df.pct_counts_lncRNA),name=f"{sample}"))

    fig.update_layout(title={#"text":f'{sample}',
                            'y':0.98, 'x':0.40, 'xanchor': 'center', 'yanchor': 'top'},
                            yaxis_title = "pct_counts_lncRNA", font=dict(size=18,color="Black"))
    fig.write_image(f"./boxplot.pdf", scale=4, width=2560, height=1440)
    # fig.show()
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="filter 5-prime bam file")
    parser.add_argument("--matrix_files", help="matrix files used to plot, split by ','", required=True)
    args = parser.parse_args()

    lncRNA_gene_set = get_lncRNA_gene_set()
    matrix_files = args.matrix_files.split(',')
    box_plot(matrix_files)

    
    