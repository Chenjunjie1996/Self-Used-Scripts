import argparse
import itertools
import os
import re
import sys
import time
from collections import OrderedDict
from html.parser import HTMLParser
from io import StringIO
from six import ensure_binary, ensure_str

import requests
from Bio import SeqIO

VDJ_C_FEATURE_TYPES = ["C", "C-REGION"]
CHAINS_WITH_ISOTYPES = ["IGH"]
# Fields that are stored in the ref fasta header
REF_FASTA_FIELDS = ["feature_id", "display_name"]
REF_FASTA_AUX_FIELDS = [
    "record_id",
    "gene_name",
    "region_type",
    "chain_type",
    "chain",
    "isotype",
    "allele_name",
]
from typing import (
    Any,
    AnyStr,
    Dict,
    NamedTuple,
    List,
    Optional,
    Tuple,
    Union,
)


IMGT_GENEDB_URL = "http://www.imgt.org/genedb/GENElect"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--genome", help="V(D)J reference package name, e.g., my-vdj-ref", required=True
    )
    parser.add_argument(
        "--species",
        default="Homo sapiens",
        help="IMGT species to fetch; e.g. 'Homo sapiens' (note capitalization and quotes)",
    )
    args = parser.parse_args()

    queries = make_queries(args.species)

    try:
        filenames = download_files(args.species, queries)

    except requests.exceptions.RequestException as ex:
        print("Failed to download from IMGT. %s\n" % ex)
        print("Failed to download all files from IMGT. Exiting.")
        sys.exit(1)

    fid = 0
    # Write IMGT fasta to a file
    with open(args.genome + "-imgt-raw.fasta", "w") as raw_imgt_fa, open(
        args.genome + "-mkvdjref-input.fasta", "w"
    ) as mkvdjref_fa:
        for filename in filenames:
            with open(filename, "r") as htmlfile:
                fa_txt = _GetLastPreTag()
                fa_txt.feed(htmlfile.read())
                fa_txt = fa_txt.last_pre

            raw_imgt_fa.write(fa_txt + "\n")

            f = StringIO(str(fa_txt))
            for record in SeqIO.parse(f, "fasta"):
                fid += 1
                feature = make_feature(record, fid)
                if feature is None:
                    continue

                mkvdjref_fa.write(convert_vdj_feature_to_fasta_entry(feature) + "\n")
                
                
def convert_vdj_feature_to_fasta_entry(feature):
    """ Generate a fasta entry from a VdjAnnotationFeature """
    fasta_fields = [getattr(feature, f) for f in REF_FASTA_FIELDS]
    aux_fields = [getattr(feature, f) for f in REF_FASTA_AUX_FIELDS]
    hdr = (
        "|".join([x.decode() if hasattr(x, "decode") else str(x) for x in fasta_fields])
        + " "
        + "|".join([x.decode() if hasattr(x, "decode") else str(x) for x in aux_fields])
    )

    return ">%s\n%s" % (hdr, ensure_str(feature.sequence))
    
    
VdjAnnotationFeature = NamedTuple(
    "VdjAnnotationFeature",
    [
        ("feature_id", int),  # Globally unique, 1-based integer; e.g., 3
        ("record_id", bytes),  # Originating transcript or cDNA id; e.g., ENST0001
        (
            "display_name",
            bytes,
        ),  # Fully qualified name suitable for display (gene*allele), e.g., TRAV1-1*01. If the allele is absent this is just gene_name.
        ("gene_name", bytes),  # E.g., TRAV1-1
        ("region_type", bytes),  # E.g., V-REGION
        ("chain_type", Optional[bytes]),  # E.g., TR or IG
        ("chain", Optional[bytes]),  # E.g., TRA or IGK
        (
            "isotype",
            Optional[bytes],
        ),  # BCR isotype. E.g., M or G2 for mu or gamma2, respectively.
        ("allele_name", Optional[bytes]),  # Allele (according to IMGT, for example); e.g., 01.
        ("sequence", Optional[bytes]),  # Sequence of feature E.g., AGCGT
    ],
)

def make_queries(species):
    """IMGT GENE-DB queries"""
    v_genes = ["TRAV", "TRBV", "TRDV", "TRGV", "IGHV", "IGKV", "IGLV"]
    v_label = "L-PART1 V-EXON"  # let requests turn the ' ' into a '+'
    v_query = "8.1"

    d_genes = ["TRBD", "TRDD", "IGHD"]
    d_label = None
    d_query = "7.2"

    j_genes = ["TRAJ", "TRBJ", "TRDJ", "TRGJ", "IGHJ", "IGKJ", "IGLJ"]
    j_label = None
    j_query = "7.2"

    c_genes = ["TRAC", "TRBC", "TRDC", "TRGC", "IGHC"]
    c_label = None
    c_query = "14.1"

    c_genes2 = ["IGKC", "IGLC"]
    c_label2 = None
    if species == "Mus musculus":
        c_query2 = "7.2"
    else:
        c_query2 = "14.1"

    return list(
        itertools.chain(
            itertools.product(v_genes, [v_label], [v_query]),
            itertools.product(d_genes, [d_label], [d_query]),
            itertools.product(j_genes, [j_label], [j_query]),
            itertools.product(c_genes, [c_label], [c_query]),
            itertools.product(c_genes2, [c_label2], [c_query2]),
        )
    )


def download_files(species, queries):
    """Download the sequences.

    e.g. http://www.imgt.org/genedb/GENElect?query=8.1+TRBV&species=Homo+sapiens&IMGTlabel=L-PART1+V-EXON
    """
    filenames = []
    for gene, label, number in queries:
        filename = "_".join((species.replace(" ", ""), number, gene)) + ".html"
        filenames.append(filename)
        if os.path.exists(filename):
            print("Already downloaded %s, skipping" % filename)
            continue

        # Note: IMGT is sensitive to the param order
        payload = OrderedDict(
            query="%s %s" % (number, gene),
            species=species,
        )

        if label:
            payload["IMGTlabel"] = label

        r = requests.get(IMGT_GENEDB_URL, params=payload.items())

        # Get the original url (pre-redirect)
        if len(r.history) > 0:
            used_url = r.history[0].url
        else:
            used_url = r.url

        print("Downloading %s to %s ..." % (used_url, filename))

        r.raise_for_status()
        with open(filename, "w") as f:
            f.write(r.text)

        # Don't hammer the server
        time.sleep(5)
    return filenames


def make_feature(record, fid):
    """Create a VdjAnnotationFeature from a record.

    Args:
        record (str): The record to parse.
        fid (int): The feature ID.

    Returns:
        vdj_reference.VdjAnnotationFeature: The feature.
    """
    row = record.description.split("|")

    region_type = get_region_type(row[4])
    if region_type is None:
        print("Warning: Unrecognized IMGT region type: %s; skipping..." % row[4])
        return None

    chain_type = infer_imgt_vdj_chain_type(row[1])
    chain = infer_imgt_vdj_chain(row[1])

    if (
        region_type in VDJ_C_FEATURE_TYPES
        and chain in CHAINS_WITH_ISOTYPES
    ):
        isotype = infer_imgt_isotype(row[1])
    else:
        isotype = None

    gene_name = re.sub("[*].*$", "", row[1])
    allele_name = infer_imgt_allele(row[1])

    return VdjAnnotationFeature(
        feature_id=fid,
        record_id=row[0],
        display_name=row[1],
        gene_name=gene_name,
        region_type=region_type,
        chain_type=chain_type,
        chain=chain,
        isotype=isotype,
        allele_name=allele_name,
        sequence=str(record.seq).upper(),
    )


class _GetLastPreTag(HTMLParser):  # pylint: disable=abstract-method
    def __init__(self):
        super().__init__()
        self._last_data = None
        self.last_pre = None

    def handle_data(self, data):
        self._last_data = data

    def handle_endtag(self, tag):
        if tag == "pre":
            self.last_pre = self._last_data


# Parse the HTML files
def get_region_type(imgt_label):
    """ Convert IMGT labels into CR region type strings """
    if imgt_label == "L-PART1+V-EXON":
        return "L-REGION+V-REGION"
    elif imgt_label in ("J-REGION", "D-REGION"):
        return imgt_label
    elif (
        "EX" in imgt_label
        or "CH" in imgt_label
        or "CL" in imgt_label
        or "M" in imgt_label
        or imgt_label == "C-REGION"
    ):
        return "C-REGION"
    else:
        return None


def infer_imgt_vdj_chain_type(gene_name):
    """ Infer e.g., TR or IG from the IMGT gene name """
    return gene_name[0:2]


def infer_imgt_vdj_chain(gene_name):
    """ Infer e.g., TRA or IGH from the IMGT gene name """
    return gene_name[0:3]


def infer_imgt_isotype(gene_name):
    """ Infer, e.g., E from IGHE """
    if len(gene_name) <= 3:
        return None
    return re.sub("[*].*$", "", gene_name)[3:]


def infer_imgt_allele(gene_name):
    return re.sub("^.*[*]", "", gene_name) or None


if __name__ == "__main__":
    main()
