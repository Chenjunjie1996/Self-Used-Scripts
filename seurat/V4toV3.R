#! /usr/bin/env Rscript

library(Seurat)
library(argparser)

argv <- arg_parser("")
argv <- add_argument(argv,"--v4rds",help = "path of v4 rds")
argv <- add_argument(argv,"--prefix",help = "prefix")
argv <- add_argument(argv,"--outdir",help = "outdir")
argv <- add_argument(argv,"--cluster",help = "colnames of the column in meta data which storing the cluster information",default = "seurat_clusters")
argv <- add_argument(argv,"--scale_data",help = "wether transfer the scale.data to v3 object,default = 'no'",default = "no")
argv <- add_argument(argv,"--normalize",help = "whether normalize data,default = yes.",default = "yes")
argv <- add_argument(argv,"--newcol",help = "wether use new color",default = "yes")
argv <- parse_args(argv)

#transfer expression matrix and meta.data
prefix <- argv$prefix
outdir <- argv$outdir
cluster <- argv$cluster
newcol <- argv$newcol

dir.create(outdir)
PRO_v4 <- readRDS(argv$v4rds)
count <- PRO_v4@assays$RNA@counts
print(class(count))
init.meta.features  <- data.frame(row.names = rownames(count))
assay.data <- CreateAssayObject(counts = count)
Key(assay.data) <- paste0(tolower("RNA"),"_")
assay.list <- list(assay.data)
names(assay.list) <- "RNA"

meta.data <- PRO_v4@meta.data
print("Create Seurat v3 object.")
PRO_v3 <- new(Class = "Seurat",assays = assay.list,meta.data = meta.data,active.assay = "RNA",active.ident = factor("tmp"),project.name = prefix,version = packageVersion(pkg = "Seurat"))
Idents(PRO_v3) <- cluster
print("active.ident is:")
print(levels(PRO_v3))
if(nrow(PRO_v4@assays$RNA@data) != 0){
	print("tranfer data")
	PRO_v3@assays$RNA@data <- PRO_v4@assays$RNA@data
	print("finish")
}
if(nrow(PRO_v4@assays$RNA@scale.data) != 0 & argv$scale_data != "no"){
	print("transfer scale data")
	PRO_v3@assays$RNA@scale.data <- PRO_v4@assays$RNA@scale.data	
	print("finish")
}
	
#transfer pca umap and tsne information
print("transfer pca umap and tsne")
PRO_v3@reductions$pca <- PRO_v4@reductions$pca
attributes(attributes(PRO_v3@reductions$pca)$class) <- list(package = "Seurat")
print("pca")
PRO_v3@reductions$umap <- PRO_v4@reductions$umap
attributes(attributes(PRO_v3@reductions$umap)$class) <- list(package = "Seurat")
print("umap")
PRO_v3@reductions$tsne <- PRO_v4@reductions$tsne
attributes(attributes(PRO_v3@reductions$tsne)$class) <- list(package = "Seurat")
print("tsne")
if(newcol == "yes"){
	source('/SGRNJ03/randd/Integratedanalysis/script/scRNA_auto/color_protocol.R')
	clustcol <- c(color_protocol,"OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold","DarkGreen","DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4","#CDB5CD","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B","#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE","#9F79EE","#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080","#E9967A","#00CDCD","#8B4500","#DDA0DD","#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49","#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1","#4F94CD","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B","#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B","#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9","#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8","#FFDAB9")
}else{
	clustcol <- c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold","DarkGreen","DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4","#CDB5CD","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B","#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE","#9F79EE","#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080","#E9967A","#00CDCD","#8B4500","#DDA0DD","#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49","#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1","#4F94CD","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B","#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B","#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9","#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8","#FFDAB9")
}
if(argv$normalize  == "yes"){
	PRO_v3 <- NormalizeData(PRO_v3)
}
print("finish")
print(class(PRO_v3))
p_1 <- DimPlot(PRO_v3,reduction = "umap",cols = clustcol,label = T)+NoLegend()
pdf(paste0(outdir,"/",prefix,"_labumap.pdf"))
print(p_1)
dev.off()
png(paste0(outdir,"/",prefix,"_labumap.png"))
print(p_1)
dev.off()

p_2 <- DimPlot(PRO_v3,reduction = "tsne",cols = clustcol,label = T)+NoLegend()
pdf(paste0(outdir,"/",prefix,"_labutsne.pdf"))
print(p_2)
dev.off()
png(paste0(outdir,"/",prefix,"_labtsne.png"))
print(p_2)
dev.off()

p_3 <- DimPlot(PRO_v3,reduction = "umap",cols = clustcol,group.by = "sample")
pdf(paste0(outdir,"/",prefix,"_sampleumap.pdf"))
print(p_3)
dev.off()
png(paste0(outdir,"/",prefix,"_sampleumap.png"))
print(p_3)
dev.off()

p_4 <- DimPlot(PRO_v3,reduction = "tsne",cols = clustcol,group.by = "sample")
pdf(paste0(outdir,"/",prefix,"_sampletsne.pdf"))
print(p_4)
dev.off()
png(paste0(outdir,"/",prefix,"_sampletsne.png"))
print(p_4)
dev.off()
saveRDS(PRO_v3,paste0(outdir,"/",prefix,"_v4_to_v3.rds"))

