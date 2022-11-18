library('Seurat')
library('ggplot2')
library("ggsci")
library('tidyverse')
library('monocle')
library('argparser')
library('stringr')
options(warn=-1)

argv <- arg_parser('')
argv <- add_argument(argv,"--rds",help="rds")
argv <- parse_args(argv)

rds <- readRDS(argv$rds)

data <- as.sparse(rds@assays$RNA@counts)
pd <- new('AnnotatedDataFrame', data = rds@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
#Construct monocle cds
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              expressionFamily = negbinomial.size())
HSMM<-monocle_cds

HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

##使用monocle选择的高变基因
disp_table <- dispersionTable(HSMM)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
HSMM <- setOrderingFilter(HSMM, disp.genes)

outP1 = stringr::str_glue('./dispersion_plot.png')
png(outP1, height=1000, width=1000)
plot_ordering_genes(HSMM)
dev.off()

HSMM <- reduceDimension(HSMM, max_components = 2, method = 'DDRTree')
HSMM <- orderCells(HSMM)

outP1 = stringr::str_glue('./pseudotime.png')
png(outP1, height=1000, width=1000)
plot_cell_trajectory(HSMM, color_by = "Pseudotime",cell_size = 1) +
  theme(legend.text = element_text(size=10),legend.background = element_rect(colour='black', size=1),
        legend.key.size = unit(50, "pt"))
dev.off()

outP1 = stringr::str_glue('./state.png')
png(outP1, height=1000, width=1000)
plot_cell_trajectory(HSMM, color_by = "State",cell_size = 1.5) + 
  guides(colour = guide_legend(override.aes = list(size=4)))+
  theme(legend.text = element_text(size=15))+ scale_color_npg()
dev.off()

outP1 = stringr::str_glue('./celltype.png')
png(outP1, height=1000, width=1000)
plot_cell_trajectory(HSMM, color_by = "cluster",cell_size = 1.5) + 
  guides(colour = guide_legend(override.aes = list(size=4)))+
  theme(legend.text = element_text(size=15))+ scale_color_npg()
dev.off()

outP1 = stringr::str_glue('./sample.png')
png(outP1, height=1000, width=1000)
plot_cell_trajectory(HSMM, color_by = "sample",cell_size = 1) + 
  guides(colour = guide_legend(override.aes = list(size=4)))+
  theme(legend.text = element_text(size=15))+ scale_color_npg()
dev.off()

meta = HSMM@phenoData@data
write.table(meta,file="./HSMMmeta.xls",sep = '\t',quote=F,row.names=T)

outP1 = stringr::str_glue('./state_barplot.png')
png(outP1, height=1000, width=1000)
ggplot(data=meta, aes(x= State, fill= cluster))+
    geom_bar(width=0.6, alpha=0.8 , colour="slategray",size=1 )+
    theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14),
          axis.text.x = element_text(size = 16,color="black"),axis.text.y = element_text(size = 14,color="black"),
          legend.text = element_text(size=15))+ theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black"))
dev.off()

#outP1 = stringr::str_glue('./car.png')
#png(outP1, height=1000, width=1000)
#plot_cell_trajectory(HSMM, color_by = "car",cell_size = 1.5) + 
#  guides(colour = guide_legend(override.aes = list(size=4)))+
#  theme(legend.text = element_text(size=15))+ scale_color_npg()
#dev.off()
