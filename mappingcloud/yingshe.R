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

Idents(rds) <- "orig.ident"
#Idents(rds) <- "sample"
rna <- subset(rds,idents = c(args$sample))

if (grepl("confident_count.tsv", args$VDJ)){
  vdj <- read.table(args$VDJ, sep='\t', header=T)
}else{
  vdj <- read.table(args$VDJ, sep=',', header=T)
}

# 转录组和VDJ版本不匹配情况
# vdj$barcode <- str_replace_all(vdj$barcode, "_", "")
df <- rna@meta.data
df$barcode <- rownames(df)

if (grepl("_", rownames(df)[1])){
  vdj$shi<-paste(args$sample,"_",sep="")
  vdj$barcode<-str_c(vdj$shi,vdj$barcode)
}

cells <- subset(vdj, productive=='True')
barcodes <- unique(cells$barcode)

filter_df <- filter(df, barcode %in% barcodes)
res <- table(filter_df$cluster)

res <- as.data.frame(res)

out_csv = stringr::str_glue("{args$outdir}/match_barcodes_celltypes_distribution.txt")
write_csv(res,file=out_csv)

#meta = rna@meta.data
#meta$Class = 'NA'
#meta[barcodes,'Class'] = 'T/BCR'
#meta <- meta %>% drop_na(nCount_RNA)
#table(meta$Class == 'T/BCR')
#rna@meta.data = meta
#rna <- RunUMAP(rna, dims = 1:20)
meta = rna@meta.data
meta$Class = 'NA'
meta$barcode = rownames(meta)
meta = meta %>%
  mutate(Class = if_else(barcode %in% barcodes,
                         true = "T/BCR",
                         false = "NA"))
meta = dplyr::select(meta, -c("barcode"))
meta <- meta %>% drop_na(nCount_RNA)
table(meta$Class == 'T/BCR')
rna@meta.data = meta
rna <- RunUMAP(rna, dims = 1:20)

outP = stringr::str_glue("{args$outdir}/cluster_umap.png")
png(outP, height=1000, width=1000)
UMAPPlot(rna,group.by='cluster',label=TRUE)
dev.off()

outP1 = stringr::str_glue("{args$outdir}/umapplot.png")
png(outP1, height=1000, width=1000)
UMAPPlot(rna,group.by='Class',cols=c('grey','red'),label=TRUE)
dev.off()

outP2 = stringr::str_glue("{args$outdir}/assign.png")
png(outP2, height=1000, width=1000)
UMAPPlot(rna,group.by='annot_full',label=TRUE,label.box=TRUE)
dev.off()

meta = rna@meta.data
outMeta = stringr::str_glue("{args$outdir}/meta.csv")
write.csv(meta,outMeta)