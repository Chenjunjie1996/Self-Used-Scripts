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
Idents(rds) <- "orig.ident"
rna <- subset(rds,idents = c(args$sample))

tcrt <- read_tsv(args$UMI_file)
tcrt$shi<-paste(args$sample,"_",sep="")
tcrt$barcode<-str_c(tcrt$shi,tcrt$barcode)
tcrt <- subset(tcrt ,select=c('barcode','UMI','tag'))

barcodes <- unique(tcrt$barcode)

df <- rna@meta.data
df$barcode <- rownames(df)

filter_df <- filter(df, barcode %in% barcodes)
res <- table(filter_df$seurat_clusters)
res <- as.data.frame(res)

out_csv = stringr::str_glue("{args$outdir}/{args$sample}_match_barcodes_celltypes_distribution.txt")
write_csv(res,file=out_csv)

meta = rna@meta.data
meta$Class = 'NA'
meta[barcodes,'Class'] = 'positive'
meta <- meta %>% drop_na(nCount_RNA)
table(meta$Class == 'positive')
rna@meta.data = meta
# rna <- RunUMAP(rna, dims = 1:20)

outP = stringr::str_glue("{args$outdir}/{args$sample}_cluster_umap.png")
png(outP, height=1000, width=1000)
UMAPPlot(rna,group.by='seurat_clusters',label=TRUE)
dev.off()

outP1 = stringr::str_glue("{args$outdir}/{args$sample}_umapplot.png")
png(outP1, height=1000, width=1000)
UMAPPlot(rna,group.by='Class',cols=c('grey','red'))
dev.off()

outP2 = stringr::str_glue("{args$outdir}/{args$sample}_assign.png")
png(outP2, height=1000, width=1000)
UMAPPlot(rna,group.by='celltype',label=TRUE,label.box=TRUE)
dev.off()

meta = rna@meta.data
outMeta = stringr::str_glue("{args$outdir}/{args$sample}_metadata.tsv")
write.table(as.data.frame(meta), file=outMeta, sep='\t')

