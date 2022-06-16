library(SingleR)
library(Seurat)
library(celldex)
library(argparser)
library(stringr)

color1 <- c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold",
  "DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4",
  "#CDB5CD","DarkGreen","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B",
  "#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE",
  "#9F79EE","#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080",
  "#E9967A","#00CDCD","#8B4500","#DDA0DD","#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49",
  "#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1","#4F94CD",
  "#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B",
  "#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B",
  "#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9",
  "#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000",
  "#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC",
  "#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8","#FFDAB9")

argv <- arg_parser('')
argv <- add_argument(argv,"--rds", help="seurat rds file")
argv <- add_argument(argv,"--sample", help="sample name")
argv <- add_argument(argv,"--species", help="human or mouse")
argv <- add_argument(argv,"--outdir", help="out dir")
argv <- parse_args(argv)

if (argv$species == "mouse") {
    # ref = MouseRNAseqData()
    ref = readRDS('/SGRNJ03/randd/cjj/Script/singleR/MouseRNAseqData.rds')
} else {
    # ref = HumanPrimaryCellAtlasData()
    ref = readRDS('/SGRNJ03/randd/cjj/Script/singleR/HumanPrimaryCellAtlasData.rds')
}

# out
out.image = stringr::str_glue("{argv$outdir}/{argv$sample}_cell_type.png")
out.rds = stringr::str_glue("{argv$outdir}/{argv$sample}.rds")

rds = readRDS(argv$rds)

df.data = GetAssayData(rds)
pred.cluster = SingleR(test = df.data, ref = ref, clusters=Idents(rds),
    labels = ref$label.main)
celltype <-data.frame(ClusterID=rownames(pred.cluster),celltype=pred.cluster$labels,stringsAsFactors = F)
rds[['celltype']]<- celltype$celltype[match(Idents(rds), celltype$ClusterID)]

out.df = stringr::str_glue("{argv$outdir}/{argv$sample}_metadata.tsv")
write.table(as.data.frame(rds@meta.data), file=out.df, sep='\t')

png(out.image, height=1000, width=1000)
DimPlot(rds,reduction="tsne", group.by="celltype", cols=color1, pt.size=1.4, label=TRUE, label.size = 8)
dev.off()

saveRDS(rds, out.rds)


