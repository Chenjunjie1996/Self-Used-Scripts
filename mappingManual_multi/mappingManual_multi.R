library(Seurat)
library(tidyverse)
library(argparse)
library(ggplot2)
library(stringr)
library(dplyr)

parser <- ArgumentParser()
parser$add_argument("--sample",help="sample name", required=TRUE)
parser$add_argument("--rds", help="scRNA rds path", required=TRUE)
parser$add_argument("--VDJ", help='VDJ match_contigs.csv path', required=TRUE)
parser$add_argument("--outdir", help='out dir', required=TRUE)
parser$add_argument("--contig_name", help='VDJ library name', required=TRUE)
args <- parser$parse_args()

#samplename = tail(strsplit(args$rds,split="/")[[1]],1)
#samplename = str_replace(samplename, '.rds', '')

rds <- readRDS(args$rds)

Idents(rds) <- "orig.ident"
rna <- subset(rds,idents = c(args$sample))
 
vdj <- read.table(args$VDJ, sep=',', header=T)
vdj$shi<-paste(args$sample,"_",sep="")
vdj$barcode<-str_c(vdj$shi,vdj$barcode)

cells <- subset(vdj, productive=='True')
barcodes <- unique(cells$barcode)

df <- rna@meta.data
df$barcode <- rownames(df)

filter_df <- filter(df, barcode %in% barcodes)
res <- table(filter_df$seurat_clusters)
res <- as.data.frame(res)

out_csv = stringr::str_glue("{args$outdir}/{args$contig_name}_match_barcodes_celltypes_distribution.txt")
write_csv(res,file=out_csv)

meta = rna@meta.data
meta$Class = 'NA'
meta[barcodes,'Class'] = 'T/BCR'
rna@meta.data = meta
rna <- RunUMAP(rna, dims = 1:20)

outP = stringr::str_glue("{args$outdir}/{args$contig_name}_cluster_umap.png")
png(outP, height=1000, width=1000)
UMAPPlot(rna,group.by='seurat_clusters',label=TRUE)
dev.off()

outP1 = stringr::str_glue("{args$outdir}/{args$contig_name}_umapplot.png")
png(outP1, height=1000, width=1000)
UMAPPlot(rna,group.by='Class',cols=c('grey','red'),label=TRUE)
dev.off()

outP2 = stringr::str_glue("{args$outdir}/{args$contig_name}_assign.png")
png(outP2, height=1000, width=1000)
UMAPPlot(rna,group.by='new_ident',label=TRUE,label.box=TRUE)
dev.off()

meta = rna@meta.data
outMeta = stringr::str_glue("{args$outdir}/{args$contig_name}_meta.csv")
write.csv(meta,outMeta)

