import pandas as pd
from celescope.tools import utils
import seaborn as sns
import matplotlib.pyplot as plt
import os
import glob
import sys
import numpy as np
from PIL import Image


def parse_contig(Cellranger, trust, trans):
    match_cell_barcodes, _match_cell_number = utils.read_barcode_file(trans)
    trust_contig = glob.glob(f'{trust}/04.summarize/*_filtered_contig.csv')[0]
    trust_sample_name = os.path.basename(os.path.abspath(trust))
    trust_contig = pd.read_csv(trust_contig)
    trust_contig = trust_contig[trust_contig['productive'] == True]
    trust_contig = trust_contig[trust_contig['barcode'].isin(match_cell_barcodes)]
    trust_contig.sort_values(by='umis', ascending=False, inplace=True)
    trust_contig.drop_duplicates(subset=['barcode', 'chain'], keep='first', inplace=True)

    cr_contig = glob.glob(f'{Cellranger}/03.assemble/all/filtered_contig_annotations.csv')[0]
    cr_sample_name = os.path.basename(os.path.abspath(Cellranger))
    cr_contig = pd.read_csv(cr_contig)
    cr_contig = cr_contig[cr_contig['productive'] == True]
    cr_contig = cr_contig[cr_contig['barcode'].isin(match_cell_barcodes)]

    return cr_contig, trust_contig, cr_sample_name, trust_sample_name


def Vlnplot(cr_contig, trust_contig, cr_sample_name, trust_sample_name):
    sns.set(font_scale=1.5, style='whitegrid')
    snsFig = sns.violinplot(y=cr_contig["length"], linewidth=5, color="skyblue")
    snsFig.set(ylim=[300, 900])
    plt.xlabel(cr_sample_name)
    plt.ylabel("Length")
    plt.title("Cellranger")
    plt.savefig("./cr-vln.png", bbox_inches='tight', dpi=300)
    plt.clf()

    sns.set(font_scale=1.5, style='whitegrid')
    snsFig = sns.violinplot(y=trust_contig["length"], linewidth=5)
    snsFig.set(ylim=[300, 900])
    plt.xlabel(trust_sample_name)
    plt.ylabel("Length")
    plt.title("Trust-4")
    plt.savefig("./trust-vln.png", bbox_inches='tight', dpi=300)

def join():
    png1 = "./cr-vln.png"
    png2 = "./trust-vln.png"
    img1, img2 = Image.open(png1), Image.open(png2)
    size1, size2 = img1.size, img2.size
    joint = Image.new('RGB', (size1[0] + size2[0], size1[1]))
    loc1, loc2 = (0, 0), (size1[0], 0)
    joint.paste(img1, loc1)
    joint.paste(img2, loc2)
    joint.save('./all-vln.png')




if __name__ == '__main__':
    Cellranger, trust, trans = sys.argv[1], sys.argv[2], sys.argv[3]
    cr_contig, trust_contig, cr_sample_name, trust_sample_name = parse_contig(Cellranger, trust, trans)
    Vlnplot(cr_contig, trust_contig, cr_sample_name, trust_sample_name)
    join()
