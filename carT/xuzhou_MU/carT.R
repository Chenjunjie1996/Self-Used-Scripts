library(Seurat)
library(tidyverse)
library(argparse)
library(ggplot2)
library(stringr)
library(dplyr)

parser <- ArgumentParser()
parser$add_argument("--sample",help="sample name", required=TRUE)
parser$add_argument("--rds", help="scRNA rds path", required=TRUE)
parser$add_argument("--UMI_tsne", help='umi or tsne file', required=TRUE)
parser$add_argument("--outdir", help='out dir', required=TRUE)
args <- parser$parse_args()

rds <- readRDS(args$rds)
Idents(rds) <- "orig.ident"
rds <- subset(rds,idents = c(args$sample))

if (grepl("UMI.csv", args$UMI_tsne)){
  print(args$UMI_tsne)
  cart <- read.table(args$UMI_tsne, sep=',', header=T)
  cart$shi<-paste(args$sample,"_",sep="")
  cart$barcode<-str_c(cart$shi,cart$barcode)
  #
  cells <- subset(cart, sum_UMI!=0)
  barcodes <- unique(cells$barcode)
}else{
  print(args$UMI_tsne)
  cart <- read_tsv(args$UMI_tsne)
  cart$shi<-paste(args$sample,"_",sep="")
  cart$barcode<-str_c(cart$shi,cart$barcode)
  cells <- subset(cart, tag!=0)
  barcodes <- unique(cells$barcode)
}

df <- rds@meta.data
df$barcode <- rownames(df)

filter_df <- filter(df, barcode %in% barcodes)
#
res <- table(filter_df$celltype)
res <- as.data.frame(res)

out_csv = stringr::str_glue("{args$outdir}/{args$sample}_match_barcodes_celltypes_distribution.txt")
write_csv(res,file=out_csv)

meta = rds@meta.data
meta$Class = 'NA'
meta[barcodes,'Class'] = 'positive'
meta <- meta %>% drop_na(nCount_RNA)
table(meta$Class == 'positive')
rds@meta.data = meta

out.df = stringr::str_glue("{args$outdir}/{args$sample}_metadata.tsv")
write.table(as.data.frame(rds@meta.data), file=out.df, sep='\t')

outP = stringr::str_glue("{args$outdir}/{args$sample}_cluster_umap.png")
png(outP, height=1000, width=1000)
UMAPPlot(rds,group.by='seurat_clusters',label=TRUE, label.size=8)
dev.off()

outP1 = stringr::str_glue("{args$outdir}/{args$sample}_umapplot.png")
png(outP1, height=1000, width=1000)
UMAPPlot(rds,group.by='Class',cols=c('grey','red'), label.size=8)
dev.off()

outP2 = stringr::str_glue("{args$outdir}/{args$sample}_assign.png")
png(outP2, height=1000, width=1000)
#
UMAPPlot(rds,group.by='celltype',label=TRUE, label.size=8)
dev.off()

print(res)
print(table(meta$celltype))
