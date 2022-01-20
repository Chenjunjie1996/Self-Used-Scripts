library(Seurat)
library(tidyverse)
library(argparse)
library(ggplot2)
library(stringr)
library(dplyr)

parser <- ArgumentParser()
parser$add_argument("--sample",help="sample name", required=TRUE)
parser$add_argument("--rds", help="scRNA rds path", required=TRUE)
parser$add_argument("--UMI_file", help='VDJ match_contigs.csv path', required=TRUE)
parser$add_argument("--outdir", help='out dir', required=TRUE)
args <- parser$parse_args()

rds <- readRDS(args$rds)
tcrt <- read_tsv(args$UMI_file)

df <- rds@meta.data
df$barcode <- rownames(df)

tcrt$shi<-paste(args$sample,"_",sep="")
tcrt$barcode<-str_c(tcrt$shi,tcrt$barcode)
tcrt <- subset(tcrt ,select=c('barcode','UMI','tag'))

df_merge <- merge(df,tcrt,by='barcode',all.x=TRUE)
rownames(df_merge) <- df_merge$barcode
df_merge <- df_merge[,-1]
df_merge[is.na(df_merge)] <- 0

out_csv = stringr::str_glue("{args$outdir}/{args$sample}_meta.txt")
write_csv(df_merge,file=out_csv)

rds@meta.data <- df_merge
outP = stringr::str_glue("{args$outdir}/{args$sample}_umap.png")
png(outP, height=1000, width=1000)
FeaturePlot(rds,features= 'UMI',cols = c("lightgrey","red"),reduction = 'umap',pt.size=1.4, label=TRUE,label.size = 12)
dev.off()

outP1 = stringr::str_glue("{args$outdir}/{args$sample}_cluster_umap.png")
png(outP1, height=1000, width=1000)
UMAPPlot(rds,group.by='seurat_clusters',pt.size=1.4,label=TRUE,label.size = 12)
dev.off()
