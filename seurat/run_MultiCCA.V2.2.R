#!/usr/bin/env Rscript
suppressMessages({
library(ggplot2)
library(reshape2)
library(argparser)
library(Cairo)
library(Seurat)
library(corrplot)
library(dplyr)
library(grid)
library(cowplot)
library(monocle)
library(dplyr)
library(grid)
})

argv <- arg_parser('')
argv <- add_argument(argv,"--expM", help="the cell gene express file list, split by ,")
argv <- add_argument(argv,"--datatype", help="the cell gene express file type 1 expM 0 10x expressfile, split by ,")
argv <- add_argument(argv,"--ugname", help="the sample belong which group  name list ,split by ,")
argv <- add_argument(argv,"--spname", help="the samples name list ,split by ,")
argv <- add_argument(argv,"--prefix", help="group name prefix")
argv <- add_argument(argv,"--group", help="compare group name:G1vsG2,G3vsG4; split by ,")
argv <- add_argument(argv,"--ccanumber", help="use 30 >= CCA number >= --dim",default=20)
argv <- add_argument(argv,"--dim", help="use cluster number 1:20",default=20)
argv <- add_argument(argv,"--mtfilter", help="filter cell percent.mito max",default=0.2)
argv <- add_argument(argv,"--ngenefiltermin", help="filter nGene  min",default=200)
argv <- add_argument(argv,"--ngenefiltermax", help="filter nGene  max",default=5000)
argv <- add_argument(argv,"--umifiltermin", help="filter nGene  min",default=0)
argv <- add_argument(argv,"--is_10X", help="10X or not.if not provided,default=no",default="no")
argv <- add_argument(argv,"--umifiltermax", help="filter nGene  max",default=30000)
argv <- add_argument(argv,"--resolution", help="the clust resolution",default=0.8)
argv <- add_argument(argv,"--logfcthreshold", help="Limit testing to genes which show,on average, at least X-fold difference (log-scale) between the two groups of cells.",default=0.25)
argv <- add_argument(argv,"--minpct", help="only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations.",default=0.1)
argv <- add_argument(argv,"--outdir", help="the output dir")
argv <- add_argument(argv,"--object", help="use 1 or 0 split by ,: Object1 martix/10X  and basic analysis; Object2 FindMarkers/no ; Object3 ConservedMarkers/no; object4 AverageExpression/no ; object5: trajectories/no ;object6: comparegroup/no",default="1,1,0,0,0,0")
argv <- parse_args(argv)
print(argv$expM)
expfile <- unlist(strsplit(argv$expM,split=","))
name <- unlist(strsplit(argv$spname,split=","))
groupname <- unlist(strsplit(argv$ugname,split=","))
gnumber<-unique(groupname)
group <- unlist(strsplit(argv$group,split=","))
object <- unlist(strsplit(argv$object,split=","))
object <- as.numeric(object)
datatype <- unlist(strsplit(argv$datatype,split=","))
datatype <- as.numeric(datatype)
cnum <- as.numeric(argv$dim)
CCAnum <- as.numeric(argv$ccanumber)
resol <- as.numeric(argv$resolution)
outdir <- argv$outdir
compare <- argv$prefix
is_10X <- argv$is_10X
mtfilter <- as.numeric(argv$mtfilter)
ngenefiltermin <- as.numeric(argv$ngenefiltermin)
ngenefiltermax <- as.numeric(argv$ngenefiltermax)
umifiltermin <- as.numeric(argv$umifiltermin)
umifiltermax <- as.numeric(argv$umifiltermax)
logfc <- as.numeric(argv$logfcthreshold)
minpct_u <- as.numeric(argv$minpct)
print(compare)
print(object)
print(groupname)
print(gnumber)
if(length(expfile) != length(name)){
   quit()
}
if(length(expfile) != length(groupname)){
   quit()
}

#stat
col1 <- colorRampPalette(c("#7F0000","red","red","#FF7F00","#FF7F00","yellow","yellow","cyan", "#007FFF", "blue", "#00007F"))
corrcol <- colorRampPalette(c("red","orange","blue","white","white"))
clustcol<-c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold","DarkGreen","DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4","#CDB5CD","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B","#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#E6E6FA","#FFDAB9"
)
col2<-colorRampPalette(c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold","DarkGreen","DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4","#CDB5CD","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B","#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#E6E6FA","#FFDAB9"))
######
   if (is_10X == "no"){
       RAW1 <- read.table(expfile[1],sep = "\t")
   }else{
       RAW1 <- Read10X(data.dir = expfile[1])
}
ccagenes.use <-c()
#x <- vector(mode='list', length=20)
x <- as.list(name)
#RAW1 <- Read10X(data.dir = expfile[1])
sampproj <- CreateSeuratObject(raw.data = RAW1, project = name[1], min.cells = 5)
sampproj <- RenameCells(sampproj, add.cell.id = "1")
sampproj@meta.data$sample <- groupname[1]
mito.genes1 <- grep(pattern="^mt-",x=rownames(x=sampproj@data),value=TRUE,ignore.case = TRUE)
percent.mito1 <- Matrix::colSums(sampproj@raw.data[mito.genes1,])/Matrix::colSums(sampproj@raw.data)
sampproj <- AddMetaData(object=sampproj,metadata=percent.mito1,col.name="percent.mito")
sampproj <- FilterCells(object=sampproj,subset.names=c("nGene","percent.mito","nUMI"),low.thresholds=c(ngenefiltermin,-Inf,umifiltermin),high.thresholds=c(ngenefiltermax,mtfilter,umifiltermax))
sampproj <-  NormalizeData(object = sampproj)
sampproj <-  ScaleData(object = sampproj,display.progress = F)
sampproj <- FindVariableGenes(sampproj, do.plot = F)
ccagenes.use <- head(rownames(sampproj@hvg.info), 5000)
print(name[1])
for (i in 2:length(expfile)){
   if (is_10X == "no"){
       RAW1 <- read.table(expfile[1],sep = "\t")
   }else{
       RAW1 <- Read10X(data.dir = expfile[1])
   }}
   sampproj_cc <- CreateSeuratObject(raw.data = RAW1, project = name[i], min.cells = 5)
   sampproj_cc <- RenameCells(sampproj_cc, add.cell.id = i)
   sampproj_cc@meta.data$sample <- groupname[i]
   mito.genes1 <- grep(pattern="^mt-",x=rownames(x=sampproj_cc@data),value=TRUE,ignore.case = TRUE)
   percent.mito1 <- Matrix::colSums(sampproj_cc@raw.data[mito.genes1,])/Matrix::colSums(sampproj_cc@raw.data)
   sampproj_cc <- AddMetaData(object=sampproj_cc,metadata=percent.mito1,col.name="percent.mito")
   sampproj_cc <- FilterCells(object=sampproj_cc,subset.names=c("nGene","percent.mito","nUMI"),low.thresholds=c(ngenefiltermin,-Inf,umifiltermin),high.thresholds=c(ngenefiltermax,mtfilter,umifiltermax))
   sampproj_cc <-  NormalizeData(object = sampproj_cc)
   sampproj_cc <-  ScaleData(object = sampproj_cc,display.progress = F)
   sampproj_cc <- FindVariableGenes(sampproj_cc, do.plot = F)
   

if(i==2){
	sampproj1<-sampproj_cc
}
  if(i==3){
	sampproj2<-sampproj_cc
	multiCCA_list<- list(sampproj,sampproj1,sampproj2)
}
  if(i==4){
	sampproj3<-sampproj_cc
	multiCCA_list<- list(sampproj,sampproj1,sampproj2,sampproj3)
}
 if(i==5){
	sampproj4<-sampproj_cc
	multiCCA_list<- list(sampproj,sampproj1,sampproj2,sampproj3,sampproj4)
}
 if(i==6){
        sampproj5<-sampproj_cc
        multiCCA_list<- list(sampproj,sampproj1,sampproj2,sampproj3,sampproj4,sampproj5)
}
 if(i==7){
        sampproj6<-sampproj_cc
        multiCCA_list<- list(sampproj,sampproj1,sampproj2,sampproj3,sampproj4,sampproj5,sampproj6)
}
 if(i==8){
        sampproj7<-sampproj_cc
        multiCCA_list<- list(sampproj,sampproj1,sampproj2,sampproj3,sampproj4,sampproj5,sampproj6,sampproj7)
}
 if(i==9){
        sampproj8<-sampproj_cc
        multiCCA_list<- list(sampproj,sampproj1,sampproj2,sampproj3,sampproj4,sampproj5,sampproj6,sampproj7,sampproj8)
}
 if(i==10){
        sampproj9<-sampproj_cc
        multiCCA_list<- list(sampproj,sampproj1,sampproj2,sampproj3,sampproj4,sampproj5,sampproj6,sampproj7,sampproj8,sampproj9)
}
 if(i==11){
        sampproj10<-sampproj_cc
        multiCCA_list<- list(sampproj,sampproj1,sampproj2,sampproj3,sampproj4,sampproj5,sampproj6,sampproj7,sampproj8,sampproj9,sampproj10)
}
 if(i==12){
        sampproj11<-sampproj_cc
        multiCCA_list<- list(sampproj,sampproj1,sampproj2,sampproj3,sampproj4,sampproj5,sampproj6,sampproj7,sampproj8,sampproj9,sampproj10,sampproj11)
}
 if(i==13){
        sampproj12<-sampproj_cc
        multiCCA_list<- list(sampproj,sampproj1,sampproj2,sampproj3,sampproj4,sampproj5,sampproj6,sampproj7,sampproj8,sampproj9,sampproj10,sampproj11,sampproj12)
}
 if(i==14){
        sampproj13<-sampproj_cc
        multiCCA_list<- list(sampproj,sampproj1,sampproj2,sampproj3,sampproj4,sampproj5,sampproj6,sampproj7,sampproj8,sampproj9,sampproj10,sampproj11,sampproj12,sampproj13)
}
 if(i==15){
        sampproj14<-sampproj_cc
        multiCCA_list<- list(sampproj,sampproj1,sampproj2,sampproj3,sampproj4,sampproj5,sampproj6,sampproj7,sampproj8,sampproj9,sampproj10,sampproj11,sampproj12,sampproj13,sampproj14)
}
 if(i==16){
        sampproj15<-sampproj_cc
        multiCCA_list<- list(sampproj,sampproj1,sampproj2,sampproj3,sampproj4,sampproj5,sampproj6,sampproj7,sampproj8,sampproj9,sampproj10,sampproj11,sampproj12,sampproj13,sampproj14,sampproj15)
}
 if(i==17){
        sampproj16<-sampproj_cc
        multiCCA_list<- list(sampproj,sampproj1,sampproj2,sampproj3,sampproj4,sampproj5,sampproj6,sampproj7,sampproj8,sampproj9,sampproj10,sampproj11,sampproj12,sampproj13,sampproj14,sampproj15,sampproj16)
}
   print(name[i])
   g.2 <- head(rownames(sampproj_cc@hvg.info), 5000)
#   ccagenes.use <- unique(c(ccagenes.use, g.2))
   ccagenes.use <- intersect(ccagenes.use,g.2)
}


#multiCCA_list<- list(argv$spname)
print(argv$spname)
gene_number<-length(ccagenes.use)
print(gene_number)
PRO<-RunMultiCCA(multiCCA_list, genes.use=ccagenes.use,num.cc = 30)
rm(multiCCA_list)
p1 <- DimPlot(object = PRO, reduction.use = "cca", group.by = "sample",pt.size = 0.1, do.return = TRUE)
p2 <- VlnPlot(object = PRO, features.plot = "CC1", group.by = "sample", x.lab.rot=TRUE,point.size.use=0.05,do.return = TRUE)
p3 <- MetageneBicorPlot(PRO, grouping.var = "sample", dims.eval = 1:25,display.progress = FALSE)
pdf(paste(outdir,'/',compare,'.CCA1.pdf',sep=''))
p1
dev.off()
CairoPNG(paste(outdir,'/',compare,'.CCA1.png',sep=''))
p1
dev.off()
pdf(paste(outdir,'/',compare,'.CCA2.pdf',sep=''))
p2
dev.off()
CairoPNG(paste(outdir,'/',compare,'.CCA2.png',sep=''))
p2
dev.off()
pdf(paste(outdir,'/',compare,'.PCElbowPlot.pdf',sep=''))
p3
dev.off()
CairoPNG(paste(outdir,'/',compare,'.PCElbowPlot.png',sep=''))
p3
dev.off()
print(CCAnum)
PRO <- AlignSubspace(PRO, reduction.type = "cca", grouping.var = "sample",dims.align = 1:CCAnum)
p1 <- VlnPlot(object = PRO, features.plot = "ACC1", group.by = "sample",x.lab.rot=TRUE,point.size.use=0.05,do.return = TRUE)
p2 <- VlnPlot(object = PRO, features.plot = "ACC2", group.by = "sample",x.lab.rot=TRUE,point.size.use=0.05,do.return = TRUE)
CairoPNG(paste(outdir,'/',compare,'.plot_ACC1.png',sep=''))
p1
dev.off()
CairoPNG(paste(outdir,'/',compare,'.plot_ACC2.png',sep=''))
p2
dev.off()
#####################
PRO <- FindClusters(object=PRO,reduction.type="cca.aligned",dims.use=1:cnum,resolution=resol,print.output=0,save.SNN=TRUE,force.recalc=TRUE)
#PRO_1<-PRO
#PRO <- SubsetData(object = PRO_1, ident.use = c("4"))
#PRO <-  FindVariableGenes(PRO,do.plot=F)
#PRO <- RunPCA(object=PRO,pc.genes=PRO@var.genes,do.print=TRUE,pcs.print=1:8,genes.print=5)
#PRO <- FindClusters(object=PRO,reduction.type="pca",dims.use=1:10,resolution=0.4,print.output=0,save.SNN=TRUE,force.recalc=TRUE)
#PrintFindClustersParams(object=PRO)
PRO <- RunTSNE(object=PRO,reduction.use = "cca.aligned",dims.use=1:cnum,do.fast=TRUE,check_duplicates = FALSE)
#PRO <- RunUMAP(PRO, reduction.use = "pca", dims.use = 1:cnum)
###############
#pdf(paste(outdir,'/',compare,'.CCA_DimHeatmap.pdf',sep=''))
#DimHeatmap(object =PRO, reduction.type = "cca", cells.use = 500,dim.use = 1:6, do.balanced = TRUE)
#dev.off()
############
clust<-summary(PRO@ident)
cluster_cell<-as.data.frame(clust)
cell_number<-sum(cluster_cell$clust)
print(cell_number)
pt_use<-0.6
if(cell_number>1000){
	pt_use<-0.4
}
if(cell_number>2500){
	pt_use<-0.3
}
if(cell_number>4000){
	pt_use<-0.2
}
if(cell_number>5500){
	pt_use<-0.15
}
if(cell_number>6500){
	pt_use<-0.1
}

#pt_use<-0.15
p1 <- TSNEPlot(PRO, do.return = T, pt.size = pt_use, group.by = "sample",plot.title = compare , colors.use=clustcol  )
pdf(paste(outdir,'/',compare,'.tSNE_samples.pdf',sep=''))
p1
dev.off()
CairoPNG(paste(outdir,'/',compare,'.tSNE_samples.png',sep=''))
p1
dev.off()
#p <- DimPlot(PRO, reduction.use = "umap",pt.size = pt_use,group.by = "sample",plot.title = compare,cols.use =clustcol)
#pdf(paste(outdir,'/',compare,'.UMAP_samples.pdf',sep=''))
#p
#dev.off()
#CairoPNG(paste(outdir,'/',compare,'.UMAP_samples.png',sep=''))
#p
#dev.off()
tsne_data<-(x = GetCellEmbeddings(object = PRO, reduction.type = "tsne", dims.use = 1:2))
#umap_data <-(x = GetCellEmbeddings(object = PRO, reduction.type = "umap", dims.use = 1:2))
write.table(tsne_data,file=paste(outdir,'/',compare,'.tsne_file.xls',sep=''),sep='\t',quote=F,row.names=T)
#write.table(umap_data,file=paste(outdir,'/',compare,'.umap_file.xls',sep=''),sep='\t',quote=F,row.names=T)
j<-0
clust<-summary(PRO@ident)
for (i in clust){j=j+1}
for (i in (j-1):0){PRO <- RenameIdent(object=PRO,old.ident.name=i,new.ident.name=i+1)}
#PRO_all<-PRO
#clust<-summary(PRO@ident)
cluster_num <- table(PRO@ident, PRO@meta.data[,"sample"])
freq_table <- prop.table(x=table(PRO@ident,PRO@meta.data[,"sample"]),margin=2)
labels <- paste('cluster ',rownames(freq_table),sep='')
rownames(freq_table)<-labels
rownames(cluster_num)<-labels
write.table(cluster_num,file=paste(outdir,'/',compare,'.CellsPerCluster.xls',sep=''),sep='\t',quote=F,row.names=T)
pdf(paste(outdir,'/',compare,'.PercentPerCell.pdf',sep=''))
mix<-(30/length(gnumber))
print(mix)
barplot(height=freq_table,width = mix,xlim=c(1,60),col =clustcol,legend = rownames(freq_table),args.legend = list(x = "right"),las=2,xlab="")
mtext(text="Samples", side=1, line=6.5)
dev.off()
CairoPNG(paste(outdir,'/',compare,'.PercentPerCell.png',sep=''))
barplot(height=freq_table,width = mix,xlim=c(1,60),col =clustcol,legend = rownames(freq_table),args.legend = list(x = "right"),las=2,xlab="")
mtext(text="Samples", side=1, line=6.5)
dev.off()
######
file.loc<-paste(outdir,'/',compare,'.cell_cluster.txt',sep='')
SaveClusters(object = PRO, file = file.loc)
#cluster<-PRO@ident
#cell_cluster<-as.data.frame(cluster)
#write.table(cell_cluster,file=paste(outdir,'/',name,'.cell_cluster.txt',sep=''),sep='\t',quote=F,row.names=T)
####
p <- TSNEPlot(object = PRO,pt.size = pt_use,plot.title = compare,colors.use=clustcol )
pdf(paste(outdir,'/',compare,'.tsne.pdf',sep=''))
p
dev.off()
CairoPNG(paste(outdir,'/',compare,'.tsne.png',sep=''))
p
dev.off()
####
p <- TSNEPlot(object = PRO,pt.size = pt_use,plot.title = compare,colors.use=clustcol,do.return = TRUE, no.legend = TRUE,do.label = TRUE)
pdf(paste(outdir,'/',compare,'.labtsne.pdf',sep=''))
p
dev.off()
CairoPNG(paste(outdir,'/',compare,'.labtsne.png',sep=''))
p
dev.off()
save_file=paste(outdir,'/',compare,'.diff_PRO.rds',sep='')
saveRDS(PRO,file=save_file)
########umap
#p <- DimPlot(PRO, reduction.use = "umap",pt.size = pt_use,plot.title = compare,cols.use =clustcol )
#pdf(paste(outdir,'/',compare,'.umap.pdf',sep=''))
#p
#dev.off()
#CairoPNG(paste(outdir,'/',compare,'.umap.png',sep=''))
#p
#dev.off()
###
#p <- DimPlot(PRO, reduction.use = "umap",pt.size = pt_use,plot.title = compare,cols.use =clustcol,do.return = TRUE, no.legend = TRUE,do.label = TRUE)
#pdf(paste(outdir,'/',compare,'.labumap.pdf',sep=''))
#p
#dev.off()
#CairoPNG(paste(outdir,'/',compare,'.labumap.png',sep=''))
#p
#dev.off()
##########
###########################FindConservedMarkers##############################################
#cluster.markers <- FindMarkers(object=PRO,ident.1=5,ident.2=12,min.pct=0.25)
#write.table(cluster.markers,file=paste(outdir,'/',compare,'.cluster5vscluster12.diffgenes.xls',sep=''),sep='\t',quote=F,row.names=T)
for (l in 1:j){
if(object[3]){ ###object3 ConservedMarkers
    tryCatch({
	cluster.markers <- FindConservedMarkers(PRO, ident.1 = l, grouping.var = "sample", print.bar = FALSE)
	write.table(cluster.markers,file=paste(outdir,'/',compare,'.cluster',l,'_FindConservedMarkers.xls',sep=''),sep='\t',quote=F,row.names=T)},error=function(e){cat("ERROR :",conditionMessage(e),"\n")}
    )
    }
if(object[2]){ ###object2 FindMarkers
#	cluster.markers <- FindMarkers(object=PRO,ident.1=l,min.pct=0.01,logfc.threshold = 0.01)
	cluster.markers <- FindMarkers(object=PRO,ident.1=l,min.pct=minpct_u,logfc.threshold = logfc)
	write.table(cluster.markers,file=paste(outdir,'/',compare,'.cluster',l,'_diffgenes.xls',sep=''),sep='\t',quote=F,row.names=T)
	}
}
###################################
if(object[4]){ ##object4 AverageExpression
cluster.averages_rep <- AverageExpression(object = PRO,add.ident = "sample")
write.table(cluster.averages_rep,file=paste(outdir,'/',compare,'.cluster_averages_repeat.xls',sep=''),sep='\t',quote=F,row.names=T)
##############################################
## use averages gene plot heatmap
cluster.averages <- AverageExpression(object = PRO,return.seurat = TRUE,show.progress = FALSE)
## use corrplot
cluster.averages1 <- AverageExpression(object = PRO)
#colnames(cluster.averages1)<-regions ## change cluster ID
write.table(cluster.averages1,file=paste(outdir,'/',compare,'.cluster_averages.xls',sep=''),sep='\t',quote=F,row.names=T)
M <- cor(cluster.averages1)
order.hc <- corrMatOrder(M, order = "hclust")
M.hc <- M[order.hc, order.hc]
pdf(paste(outdir,'/',compare,'.corrplot.pdf',sep=''))
corrplot(M.hc, method = "number", cl.lim = c(0, 1),tl.col = "black",col = rev(corrcol(50)))
dev.off()
CairoPNG(paste(outdir,'/',compare,'.corrplot.png',sep=''))
corrplot(M.hc, method = "number",cl.lim = c(0, 1),tl.col = "black",col = rev(corrcol(50)))
dev.off()
}  ##object4 AverageExpression


if(object[2]){


PRO.marker<-FindAllMarkers(object=PRO,only.pos=TRUE,min.pct=0.25,thresh.use=0.25)
markergene<-PRO.marker %>% group_by(cluster) %>% top_n(2,avg_logFC)
markergenetop10<-PRO.marker %>% group_by(cluster) %>% top_n(10,avg_logFC)

heatname=paste(compare,' top 10 marker genes heatmap',sep='')
cex_row = 10
cex_col = 10
if( j > 5){
        cex_row <- 5
        cex_col <- 7
}
pdf(paste(outdir,'/',compare,'.top10markergene_DOHeatmapplot.pdf',sep=''))
print(DoHeatmap(object=PRO,genes.use=markergenetop10$gene,slim.col.label=TRUE,group.label.rot=FALSE, title = heatname,group.spacing = 0.23,rotate.key = TRUE,cex.row=cex_row,cex.col=cex_col,col.low="blue",col.mid="white",col.high="red"))
dev.off()
CairoPNG(paste(outdir,'/',compare,'.top10markergene_DOHeatmapplot.png',sep=''))
print(DoHeatmap(object=PRO,genes.use=markergenetop10$gene,slim.col.label=TRUE,group.label.rot=FALSE, title = heatname,group.spacing = 0.23,rotate.key = TRUE,cex.row=cex_row,cex.col=cex_col,col.low="blue",col.mid="white",col.high="red"))
dev.off()

allmarker<-as.data.frame(markergene)
allmarkertop10<-as.data.frame(markergenetop10)
write.table(allmarker,file=paste(outdir,'/',compare,'.ALLmarkergenelist.xls',sep=''),sep='\t',quote=F,row.names=F)
write.table(allmarkertop10,file=paste(outdir,'/',compare,'.ALLmarkerTop10genelist.xls',sep=''),sep='\t',quote=F,row.names=F)
marker<-allmarker$gene
#marker <- c("Batf","Ccl5","Cxcr3","Eomes","Foxo1","Gzmb","Id2","Id3","Ifng","Il2ra","Il2rb","Il7r","Irf4","Kif2","Klrg1","Sell","Stat4","Tbx21","Tcf7","Zeb2","Bcat1","Bcl11b","Ccr7","Foxo1","Gzmb","Hif1a","Hk2","Id3","Ifng","Il2ra","Irf4","Myc","Rheb","Stat4","Stat5a","Tbx21","Tcf12","Xcl1","Fos","Il7r","Il10rb","Klf2","S1pr1","Bach2","Cxcr3","Id2","Lef1","Tcf7")
print(marker)
#marker<-c(marker,"Cd3g","Cd3d","Cd3e","Cd8a","Cd8b1","Thy1","Lck")
marker<-unique(marker)
#cluster.markers <- FindMarkers(object=PRO,ident.1=5,ident.2 =9,min.pct=0.25)
#write.table(cluster.markers,file=paste(outdir,'/',compare,'.cluster5vscluster9_diffgenes.xls',sep=''),sep='\t',quote=F,row.names=T)
################################
#marker<-intersect(rownames(PRO@data),marker)
print(marker)
#####
if(object[4]){
pdf(paste(outdir,'/',compare,'.cluster_Heatmapplot.pdf',sep=''))
print(DoHeatmap(object = cluster.averages,genes.use = marker,group.label.rot = TRUE, group.cex = 0,col.low="blue",col.mid="white",col.high="red"))
dev.off()
CairoPNG(paste(outdir,'/',compare,'.cluster_Heatmapplot.png',sep=''))
print(DoHeatmap(object = cluster.averages,genes.use = marker,group.label.rot = TRUE, group.cex = 0,col.low="blue",col.mid="white",col.high="red"))
dev.off()
}
####
g=0
pnum=1
plot_g<-c()
x_size=9
title_size=11
y_size=5
for (i  in marker){
	if(g==4){
pdf(paste(outdir,'/',compare,'.markergene_tsneFeatureplot',pnum,'.pdf',sep=''))
FeaturePlot(object=PRO,features.plot=plot_g,cols.use = c("lightgrey","red"),min.cutoff="q1",max.cutoff="q99",pt.size = pt_use)
dev.off()
CairoPNG(paste(outdir,'/',compare,'.markergene_tsneFeatureplot',pnum,'.png',sep=''))
FeaturePlot(object=PRO,features.plot=plot_g,cols.use = c("lightgrey","red"),min.cutoff="q1",max.cutoff="q99",pt.size = pt_use)
dev.off()
#pdf(paste(outdir,'/',compare,'.markergene_umapFeatureplot',pnum,'.pdf',sep=''))
#FeaturePlot(object=PRO,features.plot=plot_g,cols.use = c("lightgrey","red"),min.cutoff="q1",max.cutoff="q99",pt.size = pt_use, reduction.use="umap")
#dev.off()
#CairoPNG(paste(outdir,'/',compare,'.markergene_umapFeatureplot',pnum,'.png',sep=''))
#FeaturePlot(object=PRO,features.plot=plot_g,cols.use = c("lightgrey","red"),min.cutoff="q1",max.cutoff="q99",pt.size = pt_use, reduction.use="umap")
#dev.off()
pdf(paste(outdir,'/',compare,'.markergene_vlnplot',pnum,'.pdf',sep=''))
print(VlnPlot(object=PRO,features.plot=plot_g,x.lab.rot=TRUE,size.x.use=x_size,size.y.use=y_size,size.title.use=title_size,point.size.use=0,cols.use=clustcol,nCol=2,do.return=TRUE))
dev.off()
CairoPNG(paste(outdir,'/',compare,'.markergene_vlnplot',pnum,'.png',sep=''))
print(VlnPlot(object=PRO,features.plot=plot_g,x.lab.rot=TRUE,size.x.use=x_size,size.y.use=y_size,size.title.use=title_size,point.size.use=0,cols.use=clustcol,nCol=2,do.return=TRUE))
dev.off()
pnum=pnum+1
g=1
plot_g<-c(i)
	}else{
plot_g<-c(plot_g,i)
g=g+1
	}
}
if(g > 0){
pdf(paste(outdir,'/',compare,'.markergene_tsneFeatureplot',pnum,'.pdf',sep=''))
FeaturePlot(object=PRO,features.plot=plot_g,cols.use = c("lightgrey","red"),min.cutoff="q1",max.cutoff="q99",pt.size = pt_use)
dev.off()
CairoPNG(paste(outdir,'/',compare,'.markergene_tsneFeatureplot',pnum,'.png',sep=''))
FeaturePlot(object=PRO,features.plot=plot_g,cols.use = c("lightgrey","red"),min.cutoff="q1",max.cutoff="q99",pt.size = pt_use)
dev.off()
#pdf(paste(outdir,'/',compare,'.markergene_umapFeatureplot',pnum,'.pdf',sep=''))
#FeaturePlot(object=PRO,features.plot=plot_g,cols.use = c("lightgrey","red"), pt.size = pt_use ,min.cutoff="q1",max.cutoff="q99",reduction.use="umap")
#dev.off()
#CairoPNG(paste(outdir,'/',compare,'.markergene_umapFeatureplot',pnum,'.png',sep=''))
#FeaturePlot(object=PRO,features.plot=plot_g,cols.use = c("lightgrey","red"), pt.size = pt_use,min.cutoff="q1",max.cutoff="q99",reduction.use="umap")
#dev.off()
pdf(paste(outdir,'/',compare,'.markergene_vlnplot',pnum,'.pdf',sep=''))
print(VlnPlot(object=PRO,features.plot=plot_g,x.lab.rot=TRUE,size.x.use=x_size,size.y.use=y_size,size.title.use=title_size,point.size.use=0,cols.use=clustcol,nCol=2,do.return=TRUE))
dev.off()
CairoPNG(paste(outdir,'/',compare,'.markergene_vlnplot',pnum,'.png',sep=''))
print(VlnPlot(object=PRO,features.plot=plot_g,x.lab.rot=TRUE,size.x.use=x_size,size.y.use=y_size,size.title.use=title_size,point.size.use=0,cols.use=clustcol,nCol=2,do.return=TRUE))
dev.off()
g=0
plot_g<-c()
}	
pdf(paste(outdir,'/',compare,'.markergene_Dotplot.pdf',sep=''))
print(DotPlot(object=PRO,genes.plot=marker,plot.legend=TRUE,cols.use = c("blue","red"), x.lab.rot = TRUE))
dev.off()
CairoPNG(paste(outdir,'/',compare,'.markergene_Dotplot.png',sep=''))
print(DotPlot(object=PRO,genes.plot=marker,plot.legend=TRUE,cols.use = c("blue","red"), x.lab.rot = TRUE))
dev.off()

}
####
if(object[5]){ #boject5:trajectories

PRO_seurat <- importCDS(PRO)
barcode<-rownames(pData(PRO_seurat))
pfile<-data.frame(barcode=barcode,pData(PRO_seurat))
colnames(pfile)<-c("barcode","num_genes_expressed","UMI","na_use","sample","percent.mito","Cluster","Size_Factor")
head(pfile)
pfile$Cluster<-as.numeric(pfile$Cluster)+1  ##### change ID  
type_vec <- as.character(sort(unique(pfile$Cluster)))
cnumber<-max(pfile$Cluster)
pfile$Cluster <- as.character(pfile$Cluster)
type_cols_u <- clustcol[1:cnumber]
names(type_cols_u) <- type_vec
fdata<-data.frame(id=fData(PRO_seurat)$gene_short_name,gene_short_name=fData(PRO_seurat)$gene_short_name)
rownames(fdata)<-fData(PRO_seurat)$gene_short_name
my_cds<-newCellDataSet(exprs(PRO_seurat),phenoData=new("AnnotatedDataFrame", data = pfile),featureData=new("AnnotatedDataFrame", data =fdata),lowerDetectionLimit=0.5,expressionFamily=negbinomial.size())
my_cds <- estimateSizeFactors(my_cds)
my_cds <- estimateDispersions(my_cds)
my_cds <- detectGenes(my_cds, min_expr = 0.1)
pData(my_cds)$UMI <- Matrix::colSums(exprs(my_cds))
my_cds <- setOrderingFilter(my_cds,PRO@var.genes)
my_cds <- reduceDimension(my_cds,method = 'DDRTree')
my_cds <- orderCells(my_cds)
p<-plot_cell_trajectory(my_cds, color_by = "Cluster",cell_size = 0.3 )+ scale_color_manual(values = type_cols_u,name = "Cluster")
pdf(paste(outdir,'/',compare,'.trajectories_cluster.pdf',sep=''),width = 10, height = 10)
print(p)
dev.off()
CairoPNG(paste(outdir,'/',compare,'.trajectories_cluster.png',sep=''))
print(p)
dev.off()
p1<-plot_cell_trajectory(my_cds, color_by = "Pseudotime", cell_size = 0.3)
pdf(paste(outdir,'/',compare,'.Pseudotime_cluster.pdf',sep=''),width = 10, height = 10)
print(p1)
dev.off()
CairoPNG(paste(outdir,'/',compare,'.Pseudotime_cluster.png',sep=''))
print(p1)
dev.off()
write.table(pData(my_cds),file=paste(outdir,'/',compare,'.monocle_Pseudotime.xls',sep=''),quote=F,sep="\t",row.names=F)

p<-plot_cell_trajectory(my_cds, color_by = "Cluster",cell_size = 0.3)+ scale_color_manual(values = type_cols_u,name = "Cluster")+facet_wrap(~Cluster)
pdf(paste(outdir,'/',compare,'.trajectories_facet_cluster.pdf',sep=''),width = 10, height = 10)
print(p)
dev.off()
CairoPNG(paste(outdir,'/',compare,'.trajectories_facet_cluster.png',sep=''))
print(p)
dev.off()

p<-plot_cell_trajectory(my_cds, color_by = "Cluster",cell_size = 0.3)+ scale_color_manual(values = type_cols_u,name = "Cluster")+facet_wrap(~sample)
pdf(paste(outdir,'/',compare,'.trajectories_cluster_facet_sample.pdf',sep=''),width = 10, height = 10)
print(p)
dev.off()
CairoPNG(paste(outdir,'/',compare,'.trajectories_cluster_facet_sample.png',sep=''))
print(p)
dev.off()

p<-plot_complex_cell_trajectory(my_cds,color_by = "Cluster",show_branch_points = T,cell_size = 0.3)+facet_wrap(~sample, nrow = 1)+scale_size(range = c(0.2, 0.2))+scale_color_manual(values = type_cols_u,name = "Cluster")
pdf(paste(outdir,'/',compare,'.complex_trajectories_sample.pdf',sep=''),width = 10, height = 10)
print(p)
dev.off()
CairoPNG(paste(outdir,'/',compare,'.complex_trajectories_sample.png',sep=''))
print(p)
dev.off()

p<-plot_complex_cell_trajectory(my_cds,color_by = "Cluster",show_branch_points = T,cell_size = 0.3)+scale_size(range = c(0.2, 0.2))+scale_color_manual(values = type_cols_u,name = "Cluster")
pdf(paste(outdir,'/',compare,'.complex_trajectories.pdf',sep=''),width = 10, height = 10)
print(p)
dev.off()
CairoPNG(paste(outdir,'/',compare,'.complex_trajectories.png',sep=''))
print(p)
dev.off()

p<-plot_complex_cell_trajectory(my_cds,color_by = "Cluster",show_branch_points = T,cell_size = 0.3)+facet_wrap(~Cluster, nrow = 2)+scale_size(range = c(0.2, 0.2))+scale_color_manual(values = type_cols_u,name = "Cluster")
pdf(paste(outdir,'/',compare,'.complex_trajectories.facet.pdf',sep=''),width = 10, height = 10)
print(p)
dev.off()
CairoPNG(paste(outdir,'/',compare,'.complex_trajectories.facet.png',sep=''))
print(p)
dev.off()
################################################
my_pseudotime_de <- differentialGeneTest(my_cds,fullModelFormulaStr = "~sm.ns(Pseudotime)")
my_pseudotime_de %>% arrange(qval) %>% head() %>% select(id) -> my_pseudotime_gene
my_pseudotime_gene <- my_pseudotime_gene$id
my_pseudotime_de %>% arrange(qval) %>% head(30) %>% select(id) -> gene_to_cluster
print(gene_to_cluster)
gene_to_cluster <- gene_to_cluster$id
print(gene_to_cluster)
print(cnumber)
my_pseudotime_cluster <- plot_pseudotime_heatmap(my_cds[gene_to_cluster,],num_clusters = cnumber,cores = 8,show_rownames = TRUE,return_heatmap = TRUE)
pdf(paste(outdir,'/',compare,'.Pseudotime_gene_cluster.pdf',sep=''),width = 10, height = 10)
print(my_pseudotime_cluster)
dev.off()
} 
##### split test,marker gene exp
#PRO@meta.data$celltype.sample <- paste0(PRO@ident,"_",PRO@meta.data$sample)
#PRO <- StashIdent(PRO, save.name = "celltype")
#PRO <- SetAllIdent(PRO, id = "celltype.sample")
PRO_all<-PRO
if(object[2]){
PRO@meta.data$celltype.sample <- paste0(PRO@ident,"_",PRO@meta.data$sample)
PRO <- StashIdent(PRO, save.name = "celltype")
PRO <- SetAllIdent(PRO, id = "celltype.sample")
g=0
pnum=1
plot_g<-c()
print(marker)
for (i in marker){
    if(g==2){
    print(plot_g)
    pdf(paste(outdir,'/',compare,'.samples_tsneFeatureplot',pnum,'.pdf',sep=''))
    FeatureHeatmap(PRO,features.plot = plot_g, group.by = "sample",cols.use = c("lightgrey","red"),pt.size = pt_use, key.position = "top",reduction.use = "tsne")
    dev.off()
    CairoPNG(paste(outdir,'/',compare,'.samples_tsneFeatureplot',pnum,'.png',sep=''))
    print(FeatureHeatmap(PRO,features.plot = plot_g, group.by = "sample",cols.use = c("lightgrey","red"),pt.size = pt_use, key.position = "top",reduction.use = "tsne"))
    dev.off()
    pnum=pnum+1
    g=1
    plot_g<-c(i)
    }else{
    plot_g<-c(plot_g,i)
    g=g+1
    }
}
print(plot_g)
pdf(paste(outdir,'/',compare,'.samples_tsneFeatureplot',pnum,'.pdf',sep=''))
FeatureHeatmap(PRO,features.plot = plot_g, group.by = "sample",cols.use = c("lightgrey","red"),pt.size = pt_use, key.position = "top",reduction.use = "tsne")
dev.off()
CairoPNG(paste(outdir,'/',compare,'.samples_tsneFeatureplot',pnum,'.png',sep=''))
FeatureHeatmap(PRO,features.plot = plot_g, group.by = "sample",cols.use = c("lightgrey","red"),pt.size = pt_use, key.position = "top",reduction.use = "tsne")
dev.off()
#g_u<-unique(c(g_u))
dot_size=5
if(length(marker) > 10){
        dot_size=4.5
}
if(length(marker) > 15){
        dot_size=4
}
if(length(marker) > 20){
        dot_size=3.5
}
if(length(marker) > 25){
        dot_size=3
}
if(length(marker) > 30){
        dot_size=2.5
}
sdp <- SplitDotPlotGG(PRO_all, genes.plot = rev(marker), cols.use = clustcol, x.lab.rot = T, plot.legend = T, dot.scale = dot_size, do.return = T, grouping.var = "sample")
pdf(paste(outdir,'/',compare,'.SplitDotPlotGG.pdf',sep=''),width = 10, height = 10)
print(sdp)
dev.off()
CairoPNG(paste(outdir,'/',compare,'.SplitDotPlotGG.png',sep=''))
print(sdp)
dev.off()
}
####00000000000000000000000000000000000000000000000000
if(object[6]){
#####################################dif############################

#PRO@meta.data$celltype.sample <- paste0(PRO@ident,"_",PRO@meta.data$sample)
#PRO <- StashIdent(PRO, save.name = "celltype")
#PRO <- SetAllIdent(PRO, id = "celltype.sample")

for (cm in group){
comparegroup <- unlist(strsplit(cm,split="vs"))
g_u <- c()
gene_plot <- c()
for (l in 1:j){
	tryCatch({
	s1<-paste(l,'_',comparegroup[1],sep='')
	s2<-paste(l,'_',comparegroup[2],sep='')
	b.interferon.response <- FindMarkers(PRO, ident.1 = s1, ident.2 = s2,print.bar = FALSE)
	write.table(b.interferon.response,file=paste(outdir,'/',compare,'_',cm,'.cluster',l,'_diffexpressed.xls',sep=''),sep='\t',quote=F,row.names=T)
#	gene_plot <- head(rownames(b.interferon.response),4)
	gene_plot <- head(rownames(b.interferon.response),1)
	g_u <- c(g_u,head(rownames(b.interferon.response),2))
	pdf(paste(outdir,'/',compare,'_',cm,'.cluster',l,'diffgene_tsneFeatureplot.pdf',sep=''))
#	FeatureHeatmap(PRO,features.plot = gene_plot, group.by = "sample",cols.use = c("blue","red"),pt.size = 0.25, key.position = "top",max.exp = 3,reduction.use = "tsne")
	FeatureHeatmap(PRO,features.plot = gene_plot, group.by = "sample",cols.use = c("blue","red"),pt.size = 0.8, key.position = "top",reduction.use = "tsne")
	dev.off()
	CairoPNG(paste(outdir,'/',compare,'_',cm,'.cluster',l,'diffgene_tsneFeatureplot.png',sep=''))
	FeatureHeatmap(PRO,features.plot = gene_plot, group.by = "sample",cols.use = c("blue","red"),pt.size = 0.8, key.position = "top",reduction.use = "tsne")
	dev.off()
	pdf(paste(outdir,'/',compare,'_',cm,'.cluster',l,'diffgene_umapFeatureplot.pdf',sep=''))
#	FeatureHeatmap(PRO,features.plot = gene_plot, group.by = "sample", cols.use = c("blue","red"),pt.size = 0.25, key.position = "top",max.exp = 3,reduction.use = "umap")
	FeatureHeatmap(PRO,features.plot = gene_plot, group.by = "sample", cols.use = c("blue","red"),pt.size = 0.8, key.position = "top",reduction.use = "umap")
	dev.off()
	CairoPNG(paste(outdir,'/',compare,'_',cm,'.cluster',l,'diffgene_umapFeatureplot.png',sep=''))
	FeatureHeatmap(PRO,features.plot = gene_plot, group.by = "sample", cols.use = c("blue","red"),pt.size = 0.8, key.position = "top",reduction.use = "umap")
	dev.off()
	},error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
}
tryCatch({
g_u<-unique(c(g_u))
dot_size=5
if(length(g_u) > 10){
	dot_size=4.5
}
if(length(g_u) > 15){
	dot_size=4
}
if(length(g_u) > 20){
	dot_size=3.5
}
if(length(g_u) > 25){
	dot_size=3
}
if(length(g_u) > 30){
	dot_size=2.5
}
print(g_u)
sdp <- SplitDotPlotGG(PRO_all, genes.plot = rev(g_u), cols.use = clustcol, x.lab.rot = T, plot.legend = T, dot.scale = dot_size, do.return = T, grouping.var = "sample")
pdf(paste(outdir,'/',compare,'_',cm,'.SplitDotPlotGG.pdf',sep=''),width = 10, height = 10)
sdp
dev.off()
CairoPNG(paste(outdir,'/',compare,'_',cm,'.SplitDotPlotGG.png',sep=''))
sdp
dev.off()
},error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
##################################################################
}
}
#####################################################0000
#save_file=paste(outdir,'/',compare,'.diff_PRO.rds',sep='')
#saveRDS(PRO,file=save_file)
#CairoSVG(paste(outdir,'/',name,'_Region.svg',sep=''))
#p
#dev.off()
