suppressMessages({library(tidyverse)
  library(Seurat)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library("ggsci")
  library("gridExtra")})
library(Seurat)
library(monocle)
library(tidyverse)
library(patchwork)
library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)
library(tibble)
library(GSVA)
library(limma)
library(GSEABase)
library(presto)
library(DESeq2)
library(clusterProfiler)
library(pheatmap)
library(Seurat)
library(org.Hs.eg.db)
library(enrichplot)

#rds = readRDS('total.rds')
rds = readRDS("Tcells.rds")

{#~ read rds and add car  ####
rds = readRDS('total.rds')
GC0 = read.table('220309010_404_FJ_filtered_UMI.csv',sep=',',header=T)
GC1 = read.table('220303012CD_FJ_filtered_UMI.csv',sep=',',header=T)
GC2 = read.table('220305023ABC_FJ_filtered_UMI.csv',sep=',',header=T)
GC0$shi<-paste('GC_0',"_",sep="")
GC0$barcode<-str_c(GC0$shi,GC0$barcode)
GC0_cells <- subset(GC0, sum_UMI!=0)
GC0_cells <- unique(GC0_cells$barcode)
GC1$shi<-paste('GC_1',"_",sep="")
GC1$barcode<-str_c(GC1$shi,GC1$barcode)
GC1_cells <- subset(GC1, sum_UMI!=0)
GC1_cells <- unique(GC1_cells$barcode)
GC2$shi<-paste('GC_2',"_",sep="")
GC2$barcode<-str_c(GC2$shi,GC2$barcode)
GC2_cells <- subset(GC2, sum_UMI!=0)
GC2_cells <- unique(GC2_cells$barcode)
cells = c(GC0_cells,GC1_cells,GC2_cells)
rds1 = readRDS("GC_1.2_T.rds")
meta = rds@meta.data
meta$car = 'CAR-'
meta$barcode = rownames(meta)
meta = meta %>%
  mutate(car = if_else(barcode %in% cells,
                         true = "CAR+",
                         false = "CAR-"))
meta = dplyr::select(meta, -c("barcode"))
meta <- meta %>% drop_na(nCount_RNA)
table(meta$car == 'CAR+')
rds@meta.data=meta
Idents(rds) <- "cluster"
rds <- subset(rds,idents = c("CD8Teff","CD8Tex","Erythrocytes","MPs","NaiveT","Neutrophils","NK","NKT","PlasmaCells","Platelets"))
saveRDS(rds, "total.rds")

Idents(rds) <- "cluster"
rds1 <- subset(rds,idents = c("CD8Teff","CD8Tex","NaiveT","NKT"))
rds1 <-  NormalizeData(object = rds1)
rds1 <-  ScaleData(object = rds1)
rds1 <- FindVariableFeatures(object = rds1)
genes.use<- head(HVFInfo(object = rds1),2000)
rds1 <- RunPCA(object=rds1,features = VariableFeatures(object = rds1))
rds1 <- FindNeighbors(rds1, reduction = "pca", dims = 1:20)
rds1 <- FindClusters(rds1,resolution = 0.8, algorithm = 1)
rds1 <- RunTSNE(object=rds1,dims.use=1:20,do.fast=TRUE,check_duplicates = FALSE)
rds1 <- RunUMAP(rds1, reduction = "pca", dims = 1:20)
saveRDS(object = rds1, file = "Tcells.rds")
}

{#~ umap and bar plot  ####
rds = readRDS('total.rds')
Tcells = readRDS("Tcells.rds")

Idents(rds) <- "sample"
rds0 = subset(rds,idents='GC_0')
rds1 = subset(rds,idents='GC_1')
rds2 = subset(rds,idents='GC_2')

UMAPPlot(rds, group.by='cluster')
UMAPPlot(rds0, group.by='cluster')
UMAPPlot(rds1, group.by='cluster')
UMAPPlot(rds2, group.by='cluster')

UMAPPlot(rds, group.by='car',cols = c("grey", "red"))
UMAPPlot(rds0, group.by='car',cols = c("grey", "red"))
UMAPPlot(rds1, group.by='car',cols = c("grey", "red"))
UMAPPlot(rds2, group.by='car',cols = c("grey", "red"))

meta = rds@meta.data
meta0 = rds0@meta.data
meta1 = rds1@meta.data
meta2 = rds2@meta.data
# barplot
PP <- ggplot(data=meta, aes(x= sample, fill= car))+
  geom_bar(width=0.6, alpha=0.8 , colour="slategray",size=1 )+
  theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14),
        axis.text.x = element_text(size = 16,color="black"),axis.text.y = element_text(size = 14,color="black"),
        legend.text = element_text(size=15)) +scale_fill_manual(values=c("lightskyblue","#dc143c"))+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
print(PP)

PP <- ggplot(data=meta, aes(x= cluster, fill= car))+
  geom_bar(width=0.6, alpha=0.8 , colour="slategray",size=1 )+
  theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14),
        axis.text.x = element_text(size = 16,color="black"),axis.text.y = element_text(size = 14,color="black"),
        legend.text = element_text(size=15)) +scale_fill_manual(values=c("lightskyblue","#dc143c"))+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
print(PP)

PP <- ggplot(data=meta2, aes(x= cluster, fill= car))+
  geom_bar(width=0.6, alpha=0.8 , colour="slategray",size=0.1 )+
  theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14),
        axis.text.x = element_text(size = 16,color="black"),axis.text.y = element_text(size = 14,color="black"),
        legend.text = element_text(size=15)) +scale_fill_manual(values=c("lightskyblue","#dc143c"))+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
print(PP)


}

{#~ cellphone db  ####
library(psych)
library(qgraph)
library(igraph)
library(tidyverse)
library(ggplot2)
library(Seurat)
library(SeuratData)

write.table(as.matrix(rds@assays$RNA@data), 'cellphonedb_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(rds@meta.data), rds@meta.data[,'cluster', drop=F])
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unkown" #  ϸ???????в?????NA
write.table(meta_data, 'cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)

netf<- "count_network.txt"
mynet <- read.delim(paste0("./count_network.txt"), check.names = FALSE)
table(mynet$count)
mynet %>% filter(count>0) -> mynet  # ?????ᱨ??
head(mynet)
net<- graph_from_data_frame(mynet,)
plot(net)

allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB",
            "#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3",
            "#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
            "#40E0D0","#5F9EA0","#FF1493",
            "#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
            "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347",
            "#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")


karate_groups <- cluster_optimal(net)
coords <- layout_in_circle(net, order =
                             order(membership(karate_groups)))  # ???????粼??
E(net)$width  <- E(net)$count/10  # ?ߵ?Ȩ?أ???ϸ??
plot(net, edge.arrow.size=.1, 
     edge.curved=0,
     vertex.color=allcolour,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     layout = coords,
     vertex.label.cex=.7) 

net2 <- net  # ????һ?ݱ???

for (i in 1: length(unique(mynet$SOURCE)) ){
  E(net)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$color <- allcolour[i]
}  # ?Ⲩ????˭?и??õĽ????????? 

plot(net, edge.arrow.size=.1, 
     edge.curved=0,
     vertex.color=allcolour,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     layout = coords,
     vertex.label.cex=.7) 


plot(net, edge.arrow.size=.1, 
     edge.curved=0.2, # ֻ?ǵ???????????
     vertex.color=allcolour,
     vertex.frame.color="#555555",
     vertex.label.color="black",
     layout = coords,
     vertex.label.cex=1) 


dev.off()

length(unique(mynet$SOURCE)) # ?鿴??Ҫ???ƶ?????ͼ???Է??㲼??
par(mfrow=c(3,4), mar=c(.3,.3,.3,.3))

for (i in 1: length(unique(mynet$SOURCE)) ){
  net1<-net2
  E(net1)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$color <- allcolour[i]
  
  plot(net1, edge.arrow.size=0.2, 
       edge.curved=0.4,
       vertex.color=allcolour,
       vertex.frame.color="#555555",
       vertex.label.color="black",
       layout = coords,
       vertex.label.cex=1.2) 
  
}



dev.off()
length(unique(mynet$SOURCE))
par(mfrow=c(3,4), mar=c(.3,.3,.3,.3))

for (i in 1: length(unique(mynet$SOURCE)) ){
  net1<-net2
  E(net1)$count <- ""
  E(net1)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$count  <- E(net2)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$count  # ?ʼ???ʩ
  
  E(net1)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$color <- allcolour[i]
  
  plot(net1, edge.arrow.size=.1, 
       edge.curved=0.2,
       edge.label = E(net1)$count, # ???Ʊߵ?Ȩ??
       vertex.color=allcolour,
       vertex.frame.color="#555555",
       vertex.label.color="black",
       layout = coords,
       vertex.label.cex=1.2
  ) 
  
}
dev.off()


}

{#~ monocle barplot  ####
  meta=read.table("./HSMMmeta.xls")
  meta$State = as.character(meta$State)
  meta$state = NA
  meta = meta %>%
    mutate(state = if_else(State == 3,
                           true = "false",
                           false = "true")) 
  
  
  PP <- ggplot(data=meta, aes(x= State, fill= sample))+
    geom_bar(width=0.6, alpha=0.8 , colour="slategray",size=1 )+
    theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14),
          axis.text.x = element_text(size = 16,color="black"),axis.text.y = element_text(size = 14,color="black"),
          legend.text = element_text(size=15))+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
  print(PP)
}

{#~ gsva ####
  logFCcutoff=0
  adjPvalueCutoff=1
  library(GSVA)
  library(limma)
  library(GSEABase)
  
  meta = rds@meta.data
  #write.table(meta,file="tmpmeta_car.xls",sep = '\t',quote=F,row.names=T)
  meta=read.table("./tmpmeta_car.xls")
  rds@meta.data = meta
  
  Idents(rds) <- "sample"
  rt <- AverageExpression(rds)
  write.table(rt,"averageexpression.xls",sep="\t",quote=F)
  
  rt=read.table("./averageexpression.xls",sep="\t",header=T,check.names=F,row.names=1)
  rt=as.matrix(rt)
  dimnames=list(rownames(rt),colnames(rt))
  mat=matrix(as.numeric(as.matrix(rt)),nrow=nrow(rt),dimnames=dimnames)
  mat=avereps(rt)
  mat=normalizeBetweenArrays(mat)
  gmtFile = getGmt('./h.all.v7.2.symbols.gmt')
  
  gsvaOut=gsva(mat,
               gmtFile,
               min.sz=10,
               max.sz=500,
               verbose=TRUE,
               parallel.sz=1)
  gsvaOut=rbind(id=colnames(gsvaOut),gsvaOut)
  write.table(gsvaOut,file="gsvaOut.txt",sep="\t",quote=F,col.names=F)
  rt=read.table("gsvaOut.txt",sep="\t",header=T,check.names=F,row.names=1)
  
  fit=lmFit(rt)
  fit=eBayes(fit)
  all=topTable(fit, number=Inf,adjust.method="holm")
  all=rbind(id=colnames(all),all)
  write.table(all,file="all.txt",sep="\t",quote=F,col.names=F)
  
  diff <- topTable(fit, number=Inf,
                   p.value=adjPvalueCutoff, adjust="holm", lfc=logFCcutoff)
  diffName=row.names(diff)
  diff=rbind(id=colnames(diff),diff)
  write.table(diff,file="diff.xls",sep="\t",quote=F,col.names=F)
  
  hmExp=rt[diffName,]
  hmExp=rbind(id=colnames(hmExp),hmExp)
  write.table(hmExp,file="heatmap.txt",sep="\t",quote=F,col.names=F)
  
  
  rt=read.table("heatmap.txt",sep="\t",header=T,row.names=1,check.names=F)
  rt <- rt[c(1,2,4,3,6,5)]
  library(pheatmap)
  
  tiff(file="heatmap.tiff",
       width = 80,            #ͼƬ?Ŀ???
       height =140,            #ͼƬ?ĸ߶?
       units ="cm",
       compression="lzw",
       bg="white",
       res=500)
  fontsize_new = 15
  fontsize_new <- as.numeric(fontsize_new)
  
  pheatmap(rt,
           color = colorRampPalette(c( "MistyRose", "white","FireBrick"))(10),
           cluster_cols =F,
           fontsize = fontsize_new,
           fontsize_row=20,
           fontsize_col=20,angle_col = 45,
           cellwidth = 50,
           cellheight = 50,
           legend_labels = 50
  )
  dev.off()
}

{#~ gsea ####
  
# 单样本 CAR+ VS CAR-
rds = readRDS("Tcells.rds")
Idents(rds) <- "sample"
rds0 = subset(rds,idents='GC_0')
rds1 = subset(rds,idents='GC_1')
rds2 = subset(rds,idents='GC_2')

rds = rds2
exp <- rds@assays$RNA@data
dim(exp)
exp[1:4,1:4]
table(rds$car)
Idents(rds) <- 'car'
deg = FindMarkers(object = rds,ident.1 = 'CAR+',min.pct = 0.1,logfc.threshold = 0.1)
head(deg)
geneList=deg$avg_log2FC
names(geneList)= toupper(rownames(deg))
geneList=sort(geneList,decreasing = T)
head(geneList)

gmtfile ='./h.all.v7.2.symbols.gmt'
geneset <- read.gmt( gmtfile )
length(unique(geneset$term))
egmt <- GSEA(geneList, TERM2GENE=geneset, 
             minGSSize = 1,
             pvalueCutoff = 0.99,
             verbose=FALSE)
head(egmt)
egmt@result 
gsea_results_df <- egmt@result 
rownames(gsea_results_df)


m1=as.data.frame(egmt@result)

gseaplot2(egmt,geneSetID = c('HALLMARK_OXIDATIVE_PHOSPHORYLATION','HALLMARK_E2F_TARGETS','HALLMARK_MYC_TARGETS_V1'),
          pvalue_table=T,rel_heights = c(3, 1, 2),color = c('palevioletRed','lightskyblue','darkCyan'))

gseaplot2(egmt,geneSetID = c('HALLMARK_TGF_BETA_SIGNALING','HALLMARK_KRAS_SIGNALING_UP','HALLMARK_INTERFERON_GAMMA_RESPONSE'),
          pvalue_table=T,rel_heights = c(3, 1, 2),color = c('palevioletRed','lightskyblue','darkCyan'))


# 样本间
Idents(rds) <- 'car'
rds <- subset(rds,idents = c("CAR+"))
exp <- rds@assays$RNA@data
dim(exp)
exp[1:4,1:4]
table(rds$car)
Idents(rds) <- 'sample'
# , min.pct = 0,logfc.threshold = 0
deg = FindMarkers(object = rds,ident.1 = 'GC_2', ident.2 = 'GC_1' ,min.pct = 0.1,logfc.threshold = 0.1)
#deg = FindAllMarkers(rds)
head(deg)

#deg%>%
#  filter(cluster %in% c('GC_0_ZL'))->GC0
#deg%>%
#  filter(cluster %in% c('GC_1_ZL'))->GC1
#deg%>%
#  filter(cluster %in% c('GC_2_ZL'))->GC2

geneList=deg$avg_log2FC
names(geneList)= toupper(rownames(deg))
geneList=sort(geneList,decreasing = T)
head(geneList)

gmtfile ='./h.all.v7.2.symbols.gmt'
geneset <- read.gmt( gmtfile )
length(unique(geneset$term))
#egmt <- GSEA(geneList, TERM2GENE=geneset)
egmt <- GSEA(geneList, TERM2GENE=geneset, 
             minGSSize = 1,
             pvalueCutoff = 0.99,
             verbose=FALSE)
head(egmt)
egmt@result 
gsea_results_df <- egmt@result 
rownames(gsea_results_df)

m1=as.data.frame(egmt@result)

gseaplot2(egmt,geneSetID = c('HALLMARK_TNFA_SIGNALING_VIA_NFKB','HALLMARK_ESTROGEN_RESPONSE_EARLY','HALLMARK_UV_RESPONSE_DN'),
          pvalue_table=T,rel_heights = c(3, 1, 2),color = c('palevioletRed','lightskyblue','darkCyan'))

gseaplot2(egmt,geneSetID = c('HALLMARK_INTERFERON_ALPHA_RESPONSE','HALLMARK_INTERFERON_GAMMA_RESPONSE','HALLMARK_E2F_TARGETS'),
          pvalue_table=T,rel_heights = c(3, 1, 2),color = c('palevioletRed','lightskyblue','darkCyan'))

# 特殊
gseaplot2(egmt,geneSetID = c('HALLMARK_TNFA_SIGNALING_VIA_NFKB','HALLMARK_ALLOGRAFT_REJECTION'),
          pvalue_table=T,rel_heights = c(3, 1, 2),color = c('palevioletRed','lightskyblue'))
}

{#~ gene ####
  rds = readRDS("Tcells.rds")
  Idents(rds) <- "sample"
  rds0 = subset(rds,idents='GC_0')
  rds1 = subset(rds,idents='GC_1')
  rds2 = subset(rds,idents='GC_2')
  
  Idents(rds2) <- 'car'
  deg = FindMarkers(object = rds2,ident.1 = 'CAR+', min.pct = 0.1,logfc.threshold = 0.1)
  
  my_levels <- c("CAR-","CAR+")
  rds2@active.ident <- factor(x = rds2@active.ident, levels = my_levels)
  VlnPlot(rds2, pt.size=0,features = c("TNFRSF9","GZMK","STMN1","RGS1","TUBA1B","NME1"),cols=c('grey','red'))
  write.table(deg,file="C:/Users/admin/Desktop/daily work/GITHUB/徐州医科大施明/20221109/gene/单样本差异表达基因/gc2_diff_genes.xls",sep='\t',quote=F,row.names=T)
  
  # heatmap
  rds = readRDS("Tcells.rds")
  meta = rds@meta.data
  meta$celltype_car <-str_c(meta$cluster, "_", meta$car)
  rds@meta.data = meta
  Idents(rds) <- "celltype_car"
  markers = FindAllMarkers(rds)
  markers_top5 = markers%>%
    group_by(cluster)%>%
    top_n(n = 5,wt = avg_log2FC) 
  rds_scale = NormalizeData(rds)
  rds_scale = ScaleData(rds_scale)
  DoHeatmap(object = rds_scale,features = unique(markers_top5$gene),size=4)
  
  # 亚型的差异表达基因在每个样本中表达情况
  Idents(rds) <- "cluster"
  NaiveT = subset(rds, idents=c('NaiveT'))
  NKT = subset(rds, idents=c('NKT'))
  CD8Tex = subset(rds, idents=c('CD8Tex'))
  CD8Teff = subset(rds, idents=c('CD8Teff'))
  
  Idents(NaiveT) <- "car"
  markers = FindMarkers(NaiveT, ident.1 = 'CAR+')
  Idents(NaiveT) <- "sample"
  VlnPlot(NaiveT, pt.size=0, features = c('TNFRSF9','CXCL10','HLA-DRA','BIRC3','HLA-DQA1','TRAF1'))
  
  Idents(NKT) <- "car"
  markers = FindMarkers(NKT, ident.1 = 'CAR+')
  Idents(NKT) <- "sample"
  VlnPlot(NKT, pt.size=0, features = c('CD27','KLRC1','CDK2AP2','CDK2AP2','SIT1','COTL1'))
  
  Idents(CD8Tex) <- "car"
  markers = FindMarkers(CD8Tex, ident.1 = 'CAR+', min.pct = 0.1,logfc.threshold = 0.1)
  Idents(CD8Tex) <- "sample"
  VlnPlot(CD8Tex, pt.size=0, features = c('TNFRSF9','IFI6','STMN1','OCIAD2','LAG3','ISG15'))
  
  Idents(CD8Teff) <- "car"
  markers = FindMarkers(CD8Teff, ident.1 = 'CAR+')
  Idents(CD8Teff) <- "sample"
  VlnPlot(CD8Teff, pt.size=0, features = c('TNFRSF9','STAT1','BIRC3','HLA-DQA1','IFI6','TNFSF10'))
  
  # 差异表达基因变化
  rds = readRDS('Tcells.rds')
  meta = rds@meta.data
  Idents(rds) <- "sample"
  markers = FindAllMarkers(rds)
  
  GC10 = FindMarkers(object = rds,ident.1 = 'GC_1', ident.2 = 'GC_0')
  GC20 = FindMarkers(object = rds,ident.1 = 'GC_2', ident.2 = 'GC_0')
  GC21 = FindMarkers(object = rds,ident.1 = 'GC_2', ident.2 = 'GC_1')
  
  GC10$gene = rownames(GC10)
  GC20$gene = rownames(GC20)
  GC21$gene = rownames(GC21)
  
  GC10 = dplyr::select(GC10, -3,-4)
  GC20 = dplyr::select(GC20, -3,-4)
  GC21 = dplyr::select(GC21, -3,-4)
  
  colnames(GC10) <- c("p_val_10", "avg_log2fc_10","p_val_adj_10","gene")
  colnames(GC20) <- c("p_val_20", "avg_log2fc_20","p_val_adj_20","gene")
  colnames(GC21) <- c("p_val_21", "avg_log2fc_21","p_val_adj_21","gene")
  
  merge1 = merge(GC10, GC21, by='gene')
  merge2 = merge(merge1, GC20, by='gene')
  
  write.table(merge2,file="C:/Users/admin/Desktop/daily work/GITHUB/徐州医科大施明/20221109/gene/差异基因表达量变化/gene_exp_change.txt",sep="\t",quote=F)
  
  Idents(rds) <- "sample"
  # 一直上调
  VlnPlot(rds, pt.size=0, features = c('CD247','CCL5','CMC1','MALAT1','CCL4','HBB',"GZMK","RGS1","KLRC1"))
  # 先上后下
  VlnPlot(rds, pt.size=0, features = c('GZMB','PRF1','GZMH','IFI6','LAG3','ISG15','IFITM1','IFITM2','CTSW'))
  # 先下后上
  VlnPlot(NaiveT, pt.size=0, features = c('HSPA8','RPL7A','RPL3','RPL7','RPS23','RPS12'))
  # 一直下
  VlnPlot(rds, pt.size=0, features = c('TUBB','LDHA','TUBA1B','TPI1','ENO1','MKI67', "PKM", "TOP2A", "PCNA"))
}