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
args <- parser$parse_args()


rds <- readRDS(args$rds)
vdj <- read.table(args$VDJ, sep=',', header=T)

if (grepl("_", rownames(df)[1])){
  vdj$shi<-paste(args$sample,"_",sep="")
  vdj$barcode<-str_c(vdj$shi,vdj$barcode)
}

cells <- subset(vdj, productive=='True')
barcodes <- unique(cells$barcode)
df <- rds@meta.data
df$barcode <- rownames(df)

filter_df <- filter(df, barcode %in% barcodes)
res <- table(filter_df$celltype)
res <- as.data.frame(res)

out_csv = stringr::str_glue("{args$outdir}/{args$sample}_match_barcodes_celltypes_distribution.txt")
write_csv(res,file=out_csv)

meta = rds@meta.data
meta$Class = 'NA'
meta[barcodes,'Class'] = 'T/BCR'
meta <- meta %>% drop_na(nCount_RNA)
table(meta$Class == 'T/BCR')
rds@meta.data = meta

outP = stringr::str_glue("{args$outdir}/{args$sample}_cluster_umap.png")
png(outP, height=1000, width=1000)
UMAPPlot(rds,group.by='seurat_clusters',label=TRUE)
dev.off()

outP1 = stringr::str_glue("{args$outdir}/{args$sample}_umapplot.png")
png(outP1, height=1000, width=1000)
UMAPPlot(rds,group.by='Class',cols=c('grey','red'),label=TRUE)
dev.off()

outP2 = stringr::str_glue("{args$outdir}/{args$sample}_assign.png")
png(outP2, height=1000, width=1000)
if ('new_ident' %in% colnames(meta)){
  UMAPPlot(rds,group.by='new_ident',label=TRUE,label.box=TRUE)
}else{
  UMAPPlot(rds,group.by='celltype',label=TRUE,label.box=TRUE)
}
dev.off()

meta = rds@meta.data
outMeta = stringr::str_glue("{args$outdir}/{args$sample}_meta.csv")
write.csv(meta,outMeta)

print(res)
print(table(meta$celltype))