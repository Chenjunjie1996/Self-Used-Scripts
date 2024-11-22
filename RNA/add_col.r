library(Seurat)
library(dplyr)

rds = readRDS('result/Modules/1.ClustAnno/CellLabel/Main/v3/P23072705_v4_to_v3.rds')
group <- read.table(opt$group, sep='\t', head=F)
colnames(group) <- c('sample', 'group1','group2')
groupData <- data.frame(sample = rds@meta.data$sample) %>% left_join(group, by='sample')
groupData1 <- groupData[,2]
groupData2 <- groupData[,3]
names(groupData1) <- colnames(rds)
names(groupData2) <- colnames(rds)
rds$group1 <- as.character(groupData1)
rds$group2 <- as.character(groupData2)
saveRDS(rds, file = "Main_newgroup.rds")
