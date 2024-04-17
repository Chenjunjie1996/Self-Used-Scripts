library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
library(future)
library(ggplot2)
library(stringr)
library(dplyr)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--sample",help="sample name", required=TRUE)
parser$add_argument("--outdir", help='out dir', required=TRUE)
parser$add_argument("--species", help='human or mouse', required=TRUE)
parser$add_argument("--count", help='filtered count h5 file', required=TRUE)
parser$add_argument("--fragment", help='fragment file', required=TRUE)
parser$add_argument("--meta", help='cell qc metric file', required=TRUE)
args <- parser$parse_args()

counts <- Read10X_h5(args$count)
metadata <- read.csv(
  file = args$meta,
  header = TRUE,
  row.names = 1,
  sep='\t'
)
assay <- CreateChromatinAssay(
counts = counts,
sep = c(":", "-"),
# genome = "mm10",
fragments = args$fragment,
)

obj <- CreateSeuratObject(
  counts = assay,
  assay = 'peaks',
  project = 'ATAC',
  meta.data = metadata
)

if (grepl("mouse", args$species)){
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79, verbose = FALSE)
  seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
  genome(annotations) <- "mm10"
} else if (grepl("human", args$species)){
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86, verbose = FALSE)
  seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
  genome(annotations) <- "hg38"
} else {
  annotations1 <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79, verbose = FALSE)
  annotations2 <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86, verbose = FALSE)
  seqlevels(annotations1) <- paste0('mm10_chr', seqlevels(annotations1))
  genome(annotations1) <- "mm10"
  seqlevels(annotations2) <- paste0('GRCh38_chr', seqlevels(annotations2))
  genome(annotations2) <- "hg38"
  grl = GRangesList(annotations1, annotations2)
  annotations = unlist(as(grl, "GRangesList"))
}
# add the gene information to the object
Annotation(obj) <- annotations

obj <- NucleosomeSignal(object = obj)
obj <- TSSEnrichment(obj, fast = FALSE)

outP1 = stringr::str_glue("{args$outdir}/{args$sample}.png")
png(outP1, height=700, width=700)
TSSPlot(obj) + NoLegend() + ggtitle("Enrichment around TSS")
dev.off()

df = as.data.frame(obj@assays$peaks@positionEnrichment)
tss_score = max(colMeans(df))

out.df = stringr::str_glue("{args$outdir}/{args$sample}_TssScore.txt")
write.table(data.frame(Sample=args$sample, TssScore=tss_score), file=out.df, sep='\t')