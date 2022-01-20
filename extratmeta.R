library(Seurat)
library(tidyverse)
library(argparse)
library(ggplot2)
library(stringr)
library(dplyr)

parser <- ArgumentParser()
parser$add_argument("--rds", help="scRNA rds path", required=TRUE)
parser$add_argument("--outdir", help="out dir", required= TRUE)
args <- parser$parse_args()


rds = readRDS(args$rds)
meta = rds@meta.data
out_csv = stringr::str_glue("{args$outdir}/meta.csv")
write.csv(meta,file=out_csv)

