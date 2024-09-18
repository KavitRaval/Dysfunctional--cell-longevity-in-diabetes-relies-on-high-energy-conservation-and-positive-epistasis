library(Seurat)
library(Matrix)
library(ggplot2)

# Load the data
sample = readRDS("T1D_T2D_20220428.rds")

# Select the Beta cells
subsample = subset(sample, cells=colnames(sample)[sample$cell_type == "Beta"])
DefaultAssay(subsample) = "RNA"

# Add module score
genes = c("DDIT3","ATF3","PPP1R15A","RICTOR","NFE2L2")
subsample = AddModuleScore(subsample, features=list(genes), name="ModuleScore")
VlnPlot(subsample, features="ModuleScore1", pt.size=0, group.by="disease_state") + guides(fill="none") +
  labs(x=NULL) + theme(axis.text.x=element_text(angle=0, hjust=0.5))

# Select the specific cells
index1 = subsample$disease_state == "AAB"
index2 = subsample$ModuleScore1 > 0
index3 = subsample$ModuleScore1 <= 0
cell_group1 = colnames(subsample)[index1 & index2]
cell_group2 = colnames(subsample)[index1 & index3]
length(cell_group1)
length(cell_group2)
degs = FindMarkers(subsample, ident.1=cell_group1, ident.2=cell_group2, logfc.threshold=0.1)
write.csv(degs, "degs_MS_pos_vs_neg_in_AAB.csv")

index1 = subsample$disease_state == "Control"
index2 = subsample$ModuleScore1 > 0
index3 = subsample$ModuleScore1 <= 0
cell_group1 = colnames(subsample)[index1 & index2]
cell_group2 = colnames(subsample)[index1 & index3]
length(cell_group1)
length(cell_group2)
degs = FindMarkers(subsample, ident.1=cell_group1, ident.2=cell_group2, logfc.threshold=0.1)
write.csv(degs, "degs_MS_pos_vs_neg_in_Control.csv")

index1 = subsample$disease_state == "T1D"
index2 = subsample$ModuleScore1 > 0
index3 = subsample$ModuleScore1 <= 0
cell_group1 = colnames(subsample)[index1 & index2]
cell_group2 = colnames(subsample)[index1 & index3]
length(cell_group1)
length(cell_group2)
degs = FindMarkers(subsample, ident.1=cell_group1, ident.2=cell_group2, logfc.threshold=0.1)
write.csv(degs, "degs_MS_pos_vs_neg_in_T1D.csv")

index1 = subsample$disease_state == "T2D"
index2 = subsample$ModuleScore1 > 0
index3 = subsample$ModuleScore1 <= 0
cell_group1 = colnames(subsample)[index1 & index2]
cell_group2 = colnames(subsample)[index1 & index3]
length(cell_group1)
length(cell_group2)
degs = FindMarkers(subsample, ident.1=cell_group1, ident.2=cell_group2, logfc.threshold=0.1)
write.csv(degs, "degs_MS_pos_vs_neg_in_T2D.csv")
