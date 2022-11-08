#!/usr/bin/env Rscript
suppressMessages({
library(ggplot2)
library(reshape2)
library(argparser)
library(Seurat)
library(corrplot)
library(dplyr)
library(grid)
library(cowplot)
library(dplyr)
library(grid)
library(tidyverse)
})

argv <- arg_parser('')
argv <- add_argument(argv,"--expM", help="the cell gene express file list, split by ,")
argv <- add_argument(argv,"--ugname", help="the sample belong which group  name list ,split by ,")
argv <- add_argument(argv,"--spname", help="the samples name list ,split by ,")
argv <- add_argument(argv,"--prefix", help="group name prefix")
argv <- add_argument(argv,"--rds",help="exist rds")
argv <- add_argument(argv,"--group", help="compare group name:G1vsG2,G3vsG4; split by ,")
argv <- add_argument(argv,"--dim", help="use dim 1:20",default=20)
argv <- add_argument(argv,"--mtfilter", help="filter cell percent.mito max",default=0.2)
argv <- add_argument(argv,"--ngenefiltermin", help="filter nGene  min",default=200)
argv <- add_argument(argv,"--ngenefiltermax", help="filter nGene  max",default=5000)
argv <- add_argument(argv,"--umifiltermin", help="filter nGene  min",default=0)
argv <- add_argument(argv,"--umifiltermax", help="filter nGene  max",default=30000)
argv <- add_argument(argv,"--resolution", help="the clust resolution",default=0.8)
argv <- add_argument(argv,"--is_10X", help="10X or not.if not provided,default=no",default="no")
argv <- add_argument(argv,"--logfcthreshold", help="Limit testing to genes which show,on average, at least X-fold difference (log-scale) between the two groups of cells.",default=0.25)
argv <- add_argument(argv,"--type_marker_tsv",help="cell type marker tsv")
argv <- add_argument(argv,"--ident_tsv",help="cluster identity tsv")
argv <- add_argument(argv,"--minpct", help="only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations.",default=0.1)
argv <- add_argument(argv,"--outdir", help="the output dir")
argv <- add_argument(argv,"--step",help="steps to do.all steps are 123",default="1")
argv <- add_argument(argv,"--remove_contamination",help="remove doublets/unknown",default="no")
argv <- parse_args(argv)
print(argv$expM)
expfile <- unlist(strsplit(argv$expM,split=","))
type_marker_tsv <- argv$type_marker_tsv
rds <- argv$rds
step <- as.character(argv$step)
ident_tsv <- argv$ident_tsv
name <- unlist(strsplit(argv$spname,split=","))
groupname <- unlist(strsplit(argv$ugname,split=","))
gnumber<-unique(groupname)
group <- unlist(strsplit(argv$group,split=","))
cnum <- as.numeric(argv$dim)
resol <- as.numeric(argv$resolution)
remove_contamination <- argv$remove_contamination
outdir <- argv$outdir
dir.create(outdir)
compare <- argv$prefix
mtfilter <- as.numeric(argv$mtfilter)
ngenefiltermin <- as.numeric(argv$ngenefiltermin)
ngenefiltermax <- as.numeric(argv$ngenefiltermax)
umifiltermin <- as.numeric(argv$umifiltermin)
umifiltermax <- as.numeric(argv$umifiltermax)
logfc <- as.numeric(argv$logfcthreshold)
minpct_u <- as.numeric(argv$minpct)
is_10X <- argv$is_10X
print(compare)
print(groupname)
print(gnumber)
if(length(expfile) != length(name)){
   quit()
}
if(length(expfile) != length(groupname)){
   quit()
}

#check if all_data exist:wq
#checkall_data <- function(){
    #if (!exists("PRO")){
       # print ("Import RDS file...")
        #PRO <- readRDS(rds)
        #print ("Import done.")
    #}
        #return (all_data)
#}

#stat
col1 <- colorRampPalette(c("#7F0000","red","red","#FF7F00","#FF7F00","yellow","yellow","cyan", "#007FFF", "blue", "#00007F"))
corrcol <- colorRampPalette(c("red","orange","blue","white","white"))
clustcol<-c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold","DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4","#CDB5CD","DarkGreen","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B","#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE","#9F79EE","#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080","#E9967A","#00CDCD","#8B4500","#DDA0DD","#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49","#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1","#4F94CD","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B","#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B","#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9","#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8","#FFDAB9")
col2<-colorRampPalette(c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold","DarkGreen","DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4","#CDB5CD","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B","#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#E6E6FA","#FFDAB9"))
######
library_number <- length(expfile)
                                                     
if (grepl(1,step)) {
   if (is_10X == "no"){
       RAW1 <- read.table(expfile[1],sep = "\t",header=TRUE,row.names=1)
   }else{
       RAW1 <- Read10X(data.dir = expfile[1])
}
PRO1 <- CreateSeuratObject(counts = RAW1, project = name[1], min.cells = 5)
PRO1 <- RenameCells(PRO1, add.cell.id = name[1])
PRO1@meta.data$sample <- groupname[1]
#mito.genes1 <- grep(pattern="^mt-",x=rownames(PRO1@assays[["RNA"]]),value=TRUE)
#percent.mito1 <- Matrix::colSums(PRO1@assays[["RNA"]][mito.genes1,])/Matrix::colSums(PRO1@assays[["RNA"]])
#PRO1 <- AddMetaData(object=PRO1,metadata=percent.mito1,col.name="percent.mito")
PRO1[["percent.mt"]] <- PercentageFeatureSet(PRO1, pattern = "^MT-")
VlnPlot(object = PRO1,features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
PRO1 <-  NormalizeData(object = PRO1)
PRO <-  ScaleData(object = PRO1)

if (library_number>1){ 
    for (i in 2:length(expfile)){
        if (is_10X == "no"){
          RAW1 <- read.table(expfile[i],sep = "\t",header=TRUE,row.names=1)
        }else{
          RAW1 <- Read10X(data.dir = expfile[i])
        }
        PRO1 <- CreateSeuratObject(counts = RAW1, project = name[i], min.cells = 5)
        PRO1 <- RenameCells(PRO1, add.cell.id = name[i])
        PRO1@meta.data$sample <- groupname[i]
        mito.genes1 <- grep(pattern="^mt-",x=rownames(PRO1@assays[["RNA"]]),value=TRUE)
        percent.mito1 <- Matrix::colSums(PRO1@assays[["RNA"]][mito.genes1,])/Matrix::colSums(PRO1@assays[["RNA"]])
        PRO1 <- AddMetaData(object=PRO1,metadata=percent.mito1,col.name="percent.mito")
        VlnPlot(object = PRO1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
        PRO1 <-  NormalizeData(object = PRO1)
        PR01 <-  ScaleData(object = PRO1) 
        PRO <- merge(x=PRO,y=PRO1)}}
PRO <-  NormalizeData(object = PRO)
PRO <-  ScaleData(object = PRO)
PRO <- FindVariableFeatures(object = PRO)
genes.use<- head(HVFInfo(object = PRO),2000)
PRO <- RunPCA(object=PRO,features = VariableFeatures(object = PRO))
PRO <- FindNeighbors(PRO, reduction = "pca", dims = 1:cnum)
PRO <- FindClusters(PRO,resolution = resol, algorithm = 1)
PRO <- RunTSNE(object=PRO,dims.use=1:cnum,do.fast=TRUE,check_duplicates = FALSE)
PRO <- RunUMAP(PRO, reduction = "pca", dims = 1:cnum)
saveRDS(PRO, file = paste(outdir,'/',compare,'.diff_PRO.rds',sep=''))
P<-DimPlot(object = PRO, reduction = "umap",label = TRUE)
print(P)
ggplot2::ggsave(filename = paste(outdir,'/',compare,'.labumap.pdf',sep=''))
P1<-DimPlot(object = PRO, reduction = "umap",group.by = "sample")
print(P1)
ggplot2::ggsave(filename = paste(outdir,'/',compare,'.sampleumap.pdf',sep=''))

cluster_num <- table(PRO@active.ident, PRO@meta.data[,"orig.ident"])
write.table(cluster_num,file=paste(outdir,'/',compare,'.CellsNumberPerCluster.xls',sep=''),sep='\t',quote=F,row.names=T)
freq_table <- prop.table(x=table(PRO@active.ident,PRO@meta.data[,"orig.ident"]),margin=2)
write.table(freq_table,file=paste(outdir,'/',compare,'.CellsPerCluster.xls',sep=''),sep='\t',quote=F,row.names=T)
mix<-(30/length(gnumber))
if(mix<1.5){
        mix=1.5}
print(mix)
pdf(file = paste(outdir,'/',compare,'.PercentPerCell.pdf',sep=''))
barplot(height=freq_table,width = mix,xlim=c(1,60),legend = rownames(freq_table),args.legend = list(x = "right"),las=2,xlab="",col =clustcol)
dev.off()
#labels <- paste('cluster',rownames(freq_table),sep='')
#rownames(freq_table)<-labels
#rownames(cluster_num)<-labels
#dat <- melt(freq_table,varnames=c("type","sample"),value.name = "percent")
#dat$type <- as.character(dat$type)
#p <- ggplot(dat, aes(x = sample, y = percent,fill = type )) + geom_bar(stat = 'identity') + coord_flip()
subc<-levels(x=PRO@active.ident)
for (l in subc){
   cluster.markers <- FindMarkers(object=PRO,ident.1=l,min.pct=0.1,logfc.threshold = 0.25)
   write.table(cluster.markers,file=paste(outdir,'/',l,'_diffgenes.xls',sep=''),sep='\t',quote=F,row.names=T)}
dev.off()
}

if (grepl(2,step))
{
rename_dir<-paste0(outdir,"/rename/")
dir.create(rename_dir)
marker_file <- read.table(type_marker_tsv,header = TRUE,sep="\t",row.names = 1,stringsAsFactors=FALSE)
cell_name <- rownames(marker_file)
n_cell_name <- length(cell_name)
clusters <- sort(unique(PRO@active.ident))
setwd(rename_dir)
c = 0
for (cluster in clusters){
  index = 0
  for (cell in cell_name){
    index = index + 1
    features = unlist(strsplit(marker_file[index,1],","))
    for (feature in features){
      tryCatch({
        dat <- FindMarkers(object=PRO,features=feature,ident.1=cluster,min.pct = 0,logfc.threshold = -Inf)
        dat$cell_type <- cell
        dat$cluster <- cluster
        dat <- rownames_to_column(dat,var="gene")
        if (c==0){
          all_dat <- dat
          c = c + 1
        } else {
          all_dat <- rbind(all_dat,dat)
          }
        }
        ,error=function(e){print(paste0(feature," not found!")) })
    }
  }
}
all_dat <- mutate(all_dat,pct.diff=pct.1-pct.2)
write_tsv(all_dat,"type_marker_exp.tsv")
exp  <- all_dat
a <- group_by(exp,cluster,cell_type)
as <- summarize(a,avg_pct.diff=mean(pct.diff),avg_logfc=mean(avg_log2FC))
as1 <- group_by(ungroup(as),cluster)
as1 <- mutate(as1,pct_rank = rank(avg_pct.diff),
              logfc_rank= rank(avg_logfc),total_rank=pct_rank+logfc_rank)
as2 <- as1 %>% ungroup %>% group_by(cluster) %>%
  filter(total_rank==max(total_rank)) %>% arrange(as.numeric(cluster))
write_tsv(as2,"auto_cluster_type.tsv")
}

if (grepl(3,step))
{
PRO <- readRDS(rds)
cell_ident_file <- read.table(ident_tsv,header = TRUE,sep="\t",stringsAsFactors=FALSE)
current_ident <- cell_ident_file[,1]
new_ident <- cell_ident_file[,2]

if (remove_contamination=="Y"){
  bool <- grepl("unknown|doublet",new_ident)
  current_ident <- current_ident[!bool]
  new_ident <- new_ident[!bool]
  PRO <- subset(PRO,idents=current_ident)
}

if (TRUE %in% duplicated(current_ident))
{
  stop ("duplicated cluster names")
}
#c_ident <- 0:(length(current_ident)-1)
#c_ident <- paste0("C",c_ident)
#new_ident <- paste(c_ident,new_ident,sep="_")
PRO@active.ident <- plyr::mapvalues(x = PRO@active.ident, from = current_ident, to = new_ident)
PRO[["new_ident"]] <- Idents(object = PRO)

P<-DimPlot(object = PRO, reduction = "umap",label = TRUE,repel = TRUE)
print(P)
ggplot2::ggsave(filename = paste(outdir,'/',compare,'.labumap.pdf',sep=''))

P<-DimPlot(object = PRO, reduction = "umap",group.by = "sample")
print(P)
ggplot2::ggsave(filename = paste(outdir,'/',compare,'.sampleumap.pdf',sep=''))

cluster_num <- table(PRO@active.ident, PRO@meta.data[,"orig.ident"])
write.table(cluster_num,file=paste(outdir,'/',compare,'.CellsNumberPerCluster.xls',sep=''),sep='\t',quote=F,row.names=T)
freq_table <- prop.table(x=table(PRO@active.ident ,PRO@meta.data[,"orig.ident"]),margin=2)
write.table(freq_table,file=paste(outdir,'/',compare,'.CellsPerCluster.xls',sep=''),sep='\t',quote=F,row.names=T)

mixed=length(unique(PRO@meta.data$orig.ident))
mix<-(30/mixed)
if(mix<1.5){
        mix=1.5}
print(mix)
pdf(file = paste(outdir,'/',compare,'.PercentPerCell.pdf',sep=''))
barplot(height=freq_table,width = mix,xlim=c(1,60),legend = rownames(freq_table),args.legend = list(x = "right"),las=2,xlab="",col =clustcol)
dev.off()
print(mix)

saveRDS(PRO, file = paste(outdir,'/',compare,'.diff_PRO.rds',sep=''))

subc<-levels(x=PRO@active.ident)
for (l in subc){
   cluster.markers <- FindMarkers(object=PRO,ident.1=l,min.pct=0.1,logfc.threshold = 0.25)
   write.table(cluster.markers,file=paste(outdir,'/',l,'_diffgenes.xls',sep=''),sep='\t',quote=F,row.names=T)}
dev.off()
}


