suppressPackageStartupMessages({
    library("argparser")
    library("reticulate")
    library("SingleCellExperiment")
    library("Seurat")
})

argv <- arg_parser('')
argv <- add_argument(argv, '--indat', help='the cell rds')
argv <- parse_args(argv)

indat <- argv$indat
outdat <- gsub('rds$', 'h5ad', indat)
seurat <- readRDS(indat)
seurat@meta.data$cluster <- Idents(seurat)

sc <- import("scanpy")
adata_seurat <- sc$AnnData(
    X   = t(GetAssayData(seurat)),
    obs = seurat[[]],
    var = GetAssay(seurat)[[]]
)

if (!"gname" %in% colnames(seurat@meta.data)){
   if ("group"  %in% colnames(seurat@meta.data)){
      adata_seurat$obs$gname = adata_seurat$obs$group
   }else{
      adata_seurat$obs$gname = adata_seurat$obs$sample
}}

if (!"percent_mt" %in% colnames(seurat@meta.data)){
   if ("percent.mt"  %in% colnames(seurat@meta.data)){
      adata_seurat$obs$percent_mt = adata_seurat$obs$percent.mt
   }}

if (!"total_counts" %in% colnames(seurat@meta.data)){
   if ("nCount_RNA"  %in% colnames(seurat@meta.data)){
      adata_seurat$obs$total_counts = adata_seurat$obs$nCount_RNA
   }}


if (!"n_genes_by_counts" %in% colnames(seurat@meta.data)){
   if ("nFeature_RNA"  %in% colnames(seurat@meta.data)){
      adata_seurat$obs$n_genes_by_counts = adata_seurat$obs$nFeature_RNA
   }}


tryCatch({adata_seurat$obsm$update(X_umap = Embeddings(seurat, "umap"))}
,error=function(e){print("NO UMAP EMBEDDINGS")})

tryCatch({adata_seurat$obsm$update(X_tsne = Embeddings(seurat, "tsne"))}
,error=function(e){print("NO TSNE EMBEDDINGS")})

tryCatch({adata_seurat$obsm$update(X_pca = Embeddings(seurat, "pca"))}
,error=function(e){print("NO PCA EMBEDDINGS")})

tryCatch({adata_seurat$obsm$update(X_pca_harmony = Embeddings(seurat, "harmony"))}
,error=function(e){print("NO HARMONY EMBEDDINGS")})

tryCatch({adata_seurat$obsm$update(X_scanorama = Embeddings(seurat, "scanorama"))}
,error=function(e){print("NO SCANORAMA EMBEDDINGS")})

adata_seurat$layers$update(normalised = t(GetAssayData(seurat)))
adata_seurat$layers$update(raw = t(GetAssayData(seurat, slot="counts")))

adata_seurat$write(outdat)
