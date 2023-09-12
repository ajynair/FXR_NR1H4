# 2022 March 15 created 
# this is the analysis of human sc/sn spatial seq data from 2021 Andrews

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)


markerGenes <- c("ALB","APOE",
                 "ANXA4","EPCAM", "CFTR",#"KRT7","KRT8",'KRT9','KRT19','SOX9',
                 "COL4A1","LAMA2","CCBE1","HGF", #CALD1
                 "PTPRB","PECAM1","FLT1",
                 'UPK3B','MSLN',
                 "IL7R","THEMIS","FYN",
                 "NKG7",
                 "TOP2A","MKI67",
                 "EBF1","IGKC","IGHA1",
                 "MZB1","JCHAIN",
                 "HBB","HBA1","SNCA",
                 "LGALS3","CCR2",
                 # "Clec4f","Vsig4",
                 "CD163","VSIG4",'DMXL2' )


#******************************************************#******************************************************
#******************************************************reading the shared seurat objects
#******************************************************#******************************************************


# reading the shared seurat objects
# cholangiocytes
seurat <- readRDS(file="data/2021Andrews/SeuratObjects/Chol_harmony_cleaned.RDS")
DimPlot(seurat,label = T)
colnames(seurat@meta.data)
DimPlot(seurat,label = T,group.by = "orig.ident")
DimPlot(seurat,label = T,group.by = "donor")
DimPlot(seurat,label = T,group.by = "sample")
table(seurat@meta.data$sample)

DotPlot(seurat, features = markerGenes) + RotatedAxis()



#******************************************************#******************************************************
#******************************************************PLOTS 
#******************************************************#******************************************************
library(RColorBrewer)
library("ggsci")
library("scales")

# cholangiocytes
data.combined <- readRDS(file="data/2021Andrews/SeuratObjects/Chol_harmony_cleaned.RDS")
DimPlot(data.combined,label = T)
# ggsave(filename = paste0("results/figs/hum_liv_2021Andrews/hum_chol_subtypes.pdf"))

# (Gpbar1 (TGR5), Nr1i2 (PXR), Nr1i3 (CAR), Vdr (VDR), S1pr2 (S1PR2) and Rorc (RORgT) in cholangiocytes and maybe also Yap and B1-integrin expression. 
genes <- c(unlist(strsplit("NR1H4, GPBAR1, NR1I2, NR1I3, VDR, S1PR2, RORC, YAP1, ITGB1",split = ", ")))
genes <- c(unlist(strsplit("CYP27A1, CYP7A1, CYP3A4",split = ", ")))
# Gpbar1, Nr1i2, Nr1i3, VdR, S1pr2, Rorc, Itgb1
# markers required
DefaultAssay(data.combined) <- "RNA"
DotPlot(data.combined, features = markerGenes) + RotatedAxis()
genes <- intersect(genes, rownames(data.combined@assays$RNA@data))
for(gene in genes){
  FeaturePlot(object = data.combined, features = c(gene), max.cutoff = "q95", pt.size = 2) + scale_color_viridis_c()
  ggsave(filename = paste0("results/figs/hum_liv_2021Andrews/hum_chol_gene_",gene,".pdf"))
}
# data.combined@assays$
FeaturePlot(object = data.combined, features = c(""), max.cutoff = "q95", pt.size = 2) + scale_color_viridis_c()

VlnPlot(data.combined,features = c("NR1H4"))
genes <- c(unlist(strsplit("NR1H4, YAP1, ITGB1",split = ", ")))
genes <- intersect(genes, rownames(data.combined@assays$RNA@data))
for(gene in genes){
  VlnPlot(data.combined,features = c(gene))
  ggsave(filename = paste0("results/figs/hum_liv_2021Andrews/vln_hum_chol_gene_",gene,".pdf"))
}


# wLiver
data.combined <- readRDS(file="data/2021Andrews/SeuratObjects/Integrated_with_Subannotations.rds")
DimPlot(data.combined,label = T,group.by = "Manual_Annotation")
# ggsave(filename = paste0("results/figs/hum_liv_2021Andrews/hum_wLiver_celltypes.pdf"))
DotPlot(data.combined, group.by = "Manual_Annotation", features = markerGenes) + RotatedAxis()
# ggsave(filename = paste0("results/figs/hum_liv_2021Andrews/hum_wLiver_celltype_markers.pdf"))

genes <- c(unlist(strsplit("NR1H4, GPBAR1, NR1I2, NR1I3, VDR, S1PR2, RORC, YAP1, ITGB1",split = ", ")))
genes <- intersect(genes, rownames(data.combined@assays$RNA@data))
for(gene in genes){
  FeaturePlot(object = data.combined, features = c(gene), pt.size = 0.5, order = T) + scale_color_viridis_c()
  ggsave(filename = paste0("results/figs/hum_liv_2021Andrews/hum_wLiver_gene_",gene,".pdf"))
}
# data.combined@assays$
FeaturePlot(object = data.combined, features = c(""), max.cutoff = "q95", pt.size = 2) + scale_color_viridis_c()

VlnPlot(data.combined,features = c("NR1H4"), group.by = "Manual_Annotation")
genes <- c(unlist(strsplit("NR1H4, YAP1, ITGB1",split = ", ")))
genes <- intersect(genes, rownames(data.combined@assays$RNA@data))
for(gene in genes){
  VlnPlot(data.combined,features = c(gene), group.by = "Manual_Annotation")
  ggsave(filename = paste0("results/figs/hum_liv_2021Andrews/vln_hum_wLiver_gene_",gene,".pdf"))
}
