
# author : Ceren Pajanoja
# date   : September 2022

# Script definition: Analysis for Figure 3 and Sup Figure 3

# Merge all stages (after SoupX and Doublet)

library(dplyr)
library(Seurat)
library(pheatmap)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(RColorBrewer)
library(tidyr)
library(viridis)

# Merge all stages, filter and then downstream analysis 
nc.data <- merge(nc.dataH5,c(nc.data1, nc.data4, nc.data7), add.cell.ids = c("HH5","1som","4som","7som"))

#:::::::::::::::::::::::::::::::::::::::    FIRST STEPS   ::::::::::::::::::::::::::::::::::::::::::::
mito.genes1 <- grep(pattern = "^ND1$|^ND2$|^COX1$|^COX2$|^ATP8$|^ATP6$|^COX3$|^ND3$|^ND4L$|^ND4$|^ND5$|^ND6$|^CYB$", x = rownames(x = nc.data), value = TRUE)
mito.genes2 <- grep(pattern = '^MT-', x = rownames(x = nc.data), value = TRUE)
mito.genes <- c(mito.genes1,mito.genes2)
percent.mt <- Matrix::colSums(nc.data[mito.genes, ])/Matrix::colSums(nc.data)
nc.data <- AddMetaData(object = nc.data, metadata = percent.mt, col.name = "percent.mt")
# cell cycle difference
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
nc.data <- CellCycleScoring(nc.data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
nc.data$CC.Difference <- nc.data$S.Score - nc.data$G2M.Score

# check counts, genes, and percent mito distribution before filtering
VlnPlot(nc.data, c("nFeature_RNA","nCount_RNA","percent.mt"), group.by = "orig.ident", pt.size = 0,
        cols = c('HH5' = '#e6c029', '1som' = '#fa7148', '4som' = '#4cd2ff','7som' = '#2f7cff'))
# filter below
nc.data <- subset(nc.data, subset = nFeature_RNA > 1500 & percent.mt < 1 & nFeature_RNA < 4000)
VlnPlot(nc.data, c("nFeature_RNA","nCount_RNA","percent.mt"), group.by = "orig.ident", pt.size = 0,
              cols = c('HH5' = '#e6c029', '1som' = '#fa7148', '4som' = '#4cd2ff','7som' = '#2f7cff'))

# Perform normalization, PCA and UMAP
DefaultAssay(object = nc.data) <- "RNA"
nc.data<- NormalizeData(nc.data)
nc.data<- FindVariableFeatures(nc.data)
nc.data<- ScaleData(object = nc.data, verbose = FALSE,vars.to.regress =c("percent.mt","CC.Difference"))
nc.data <- RunPCA(object = nc.data, npcs = 30, verbose = FALSE)
ElbowPlot(nc.data)
nc.data <- FindNeighbors(object = nc.data, dims = 1:15)
nc.data <- FindClusters(object = nc.data, resolution = 0.2, method = "pca")
nc.data<- RunUMAP(object = nc.data, reduction = "pca", 
                        dims = 1:15)
#set idents into order 
my_levels <- c("HH5","1som","4som","7som")
nc.data@active.ident <- factor(x = nc.data@active.ident, levels = my_levels)
nc.data@meta.data[["old.ident"]] <-Idents(nc.data)
set.seed(42)
DimPlot(nc.data, reduction = "umap", label = F, label.size = 6, pt.size= 1, group.by = "old.ident", cols = c('HH5' = '#e6c029', '1som' = '#fa7148', '4som' = '#4cd2ff','7som' = '#2f7cff'))

#:::::::::::::::::::::::::::::::::::   FIND DE GENES   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Find marker genes for each cluster
nc.data.markers <- FindAllMarkers(object = nc.data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.2)
# Heatmap of the top 15 markers in each cluster
top8 <- nc.data.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
DoHeatmap(object = nc.data, features = top8$gene) + scale_fill_gradientn(colours = rev(brewer.pal(11,"RdBu"))) + NoLegend()
write.csv(nc.data.markers, file="DE_genes_4stages.csv", row.names = T)

#::::::::::::::::::::::::::::::::::::  RENAME CLUSTERS BASED ON MARKERS  :::::::::::::::::::::::::::::::::::::::::
nc.data[["Cell.Type"]] <- nc.data$seurat_clusters
Idents(nc.data) <- "Cell.Type"
current.origcluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
new.origcluster.ids <- c("Ectoderm", "Ectoderm", "Ectoderm", "Endoderm", "Mesoderm", "Mesoderm", "Ectoderm", "Mesoderm", "Mesoderm", "Notochord")
# embed them in
nc.data@meta.data$Cell.Type <- plyr::mapvalues(x = nc.data@meta.data$Cell.Type, from = current.origcluster.ids, to = new.origcluster.ids)
nc.data@active.ident <- nc.data$Cell.Type

DimPlot(nc.data, reduction = "umap", label = F, pt.size= 1, cols = c('Ectoderm' = 'mediumseagreen', 'Endoderm' = 'tan1', 'Mesoderm' = 'violet','Notochord' = 'mediumpurple3'))

# Feature Plots: 
#Ectoderm
patchwork::wrap_plots(FeaturePlot(nc.data, features=c("TFAP2A","SOX2"), combine=FALSE,order= TRUE,pt.size = 1.1,ncol=2)) & theme() & scale_color_gradient(low = "lightgray", high = "mediumseagreen")
#Mesoderm
patchwork::wrap_plots(FeaturePlot(nc.data, features=c("ALX1","TWIST1"), combine=FALSE,order= TRUE,pt.size = 1.1,ncol=3)) & theme() & scale_color_gradient(low = "lightgray", high = "violet")
#Endoderm
patchwork::wrap_plots(FeaturePlot(nc.data, features=c("SOX17"), combine=FALSE,order= TRUE,pt.size = 1.1,ncol=3)) & theme() & scale_color_gradient(low = "lightgray", high = "tan1")
#Notochord
patchwork::wrap_plots(FeaturePlot(nc.data, features=c("CHRD"), combine=FALSE,order= TRUE,pt.size = 1.1,ncol=3)) & theme() & scale_color_gradient(low = "lightgray", high = "mediumpurple3")
#Pluri
patchwork::wrap_plots( FeaturePlot(nc.data, features=c("NANOG","KLF4","Pou5f3"), combine=FALSE,order= TRUE,pt.size = 1.1,ncol=3)) & theme() & scale_color_gradient(low = "lightgray", high = "deepskyblue3")
plot <- DotPlot(nc.data,col.min=0 ,features = c("TFAP2A","DLX5","CLDN1","SOX2","NES","MYCN","SOX17","KRT7","HHEX","ALX1","TWIST1","HAND2","PITX2","CHRD","TBXT"))+ scale_colour_gradient2(low = "#1c0333", mid = "#bec4dc", high = "#f7aa3e")+ scale_size(range = c(1.1, 8)) & theme(text = element_text(size = 10, face = "bold"))
plot + coord_flip()

#::::::::::::::::::::::::::::::::::::  SUBSET CLUSTERS   :::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Subset only Ectoderm cluster
nc.data_filtered <- subset(nc.data, idents = c("Ectoderm"))
DefaultAssay(object = nc.data_filtered) <- "RNA"
nc.data_filtered<- NormalizeData(nc.data_filtered)
nc.data_filtered<- FindVariableFeatures(nc.data_filtered)
nc.data_filtered<- ScaleData(object = nc.data_filtered, verbose = FALSE,vars.to.regress =c("percent.mt","CC.Difference"))
nc.data_filtered <- RunPCA(object = nc.data_filtered, npcs = 30, verbose = FALSE)
ElbowPlot(nc.data_filtered)
nc.data_filtered <- FindNeighbors(object = nc.data_filtered, dims = 1:15)
nc.data_filtered <- FindClusters(object = nc.data_filtered, resolution = 0.2, method = "pca")
nc.data_filtered<- RunUMAP(object = nc.data_filtered, reduction = "pca", dims = 1:10)

set.seed(42)
DimPlot(nc.data_filtered, reduction = "umap", label = F, label.size = 6, pt.size= 1, group.by = "old.ident", cols = c('HH5' = '#e6c029', '1som' = '#fa7148', '4som' = '#4cd2ff','7som' = '#2f7cff'))

#::::::::::::::::::::::::::::::::::::  CALCULTE GENE POSITIVE CELLS IN EACH CLUSTER   :::::::::::::::::::::::::::::::
Idents(nc.data_filtered) <- "orig.ident"
# Name each stage as a variable 
HH5  <-subset(nc.data_filtered, idents = c("HH5"))
Som1 <-subset(nc.data_filtered, idents = c("1som"))
Som4 <-subset(nc.data_filtered, idents = c("4som"))
Som7 <-subset(nc.data_filtered, idents = c("7som"))

myList <- list(HH5,Som1,Som4,Som7)
geneList <- c("NANOG","KLF4","Pou5f3","TFAP2A","SOX2")
data <- data.frame()
counter =1
for (i in myList) {
  for(k in geneList) {
    pc <- length(which(FetchData(i, vars = k) > 0))
    tc <- as.numeric(table(i@active.ident))
    data[counter,k] <- round(pc/tc,3)
  }
  counter <-counter+1
}
data$Stages <- c("HH5","Som1","Som4","Som7")
data <- data  %>%
  gather("goi","counts",-"Stages") 
#plot bars showing cell counts for each gene positive cell
ggplot(data, aes(fill=goi, y=counts, x=Stages)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis(discrete = T, option = "D") +
  ggtitle("Gene Counts / Cluster") +
  theme_bw()+
  xlab("Cluster IDs") +
  ylab("Cell Counts") +
  theme(
    axis.title.x = element_text(color="black", size=16, face="bold"),
    axis.title.y = element_text(color="black", size=16, face="bold"),
    axis.text.x  = element_text(size = 15),
    axis.text.y = element_text(size = 13))

#::::::::::::::::::::::::::::::::::::  CALCULTE MODULE SCORES   :::::::::::::::::::::::::::::::::::::::
myMarkers_PE  <- read.csv("gene_lists/Pan_ecto_genes.csv")
myMarkers_NC   <- read.csv("gene_lists/Neural_crest_genes.csv")
myMarkers_N   <- read.csv("gene_lists/Neural_genes.csv")
temp_labels <- nc.data_filtered@meta.data %>%
  group_by(old.ident) %>%
  tally()
nc.data_filtered <- AddModuleScore(object = nc.data_filtered,features = list(as.character(myMarkers_PE$x)),ctrl = 5, name = 'pe')
nc.data_filtered <- AddModuleScore(object = nc.data_filtered,features = list(as.character(myMarkers_NC$x)),ctrl = 5, name = 'nc')
nc.data_filtered <- AddModuleScore(object = nc.data_filtered,features = list(as.character(myMarkers_N$x)),ctrl = 5, name = 'n')

patchwork::wrap_plots( FeaturePlot(nc.data_filtered3, features="pe1", repel = TRUE,order= TRUE,pt.size = 1, combine=FALSE)) &
  theme() & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
patchwork::wrap_plots( FeaturePlot(nc.data_filtered3, features="pe1", split.by = "old.ident",repel = TRUE,order= TRUE,pt.size = 1, combine=FALSE)) &
  theme() & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
patchwork::wrap_plots( FeaturePlot(nc.data_filtered3, features="nc1", repel = TRUE,order= TRUE,pt.size = 1, combine=FALSE)) &
  theme() & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
patchwork::wrap_plots( FeaturePlot(nc.data_filtered3, features="nc1", split.by = "old.ident",repel = TRUE,order= TRUE,pt.size = 1, combine=FALSE)) &
  theme() & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
patchwork::wrap_plots( FeaturePlot(nc.data_filtered3, features="n1", repel = TRUE,order= TRUE,pt.size = 1, combine=FALSE)) &
  theme() & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
patchwork::wrap_plots( FeaturePlot(nc.data_filtered3, features="n1", split.by = "old.ident",repel = TRUE,order= TRUE,pt.size = 1, combine=FALSE)) &
  theme() & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))

# Extract cells with high scores below :::::
cell_ids <- as.data.frame(colnames(as.data.frame(nc.data_filtered3@assays[["RNA"]]@data)))
colnames(cell_ids) <-"cell_ids"
groupz <- as.data.frame(nc.data_filtered3@meta.data[["old.ident"]])
colnames(groupz) <-"group"
# Choose which scores you want to compare 
a <- as.data.frame(nc.data_filtered3@meta.data[["nc1"]])
colnames(a) <- "module"
b <- as.data.frame(nc.data_filtered3@meta.data[["pe1"]])
colnames(b) <- "panecto"
# Merge them into dataset, then remove variables
data <- do.call("cbind",list(cell_ids,groupz,a,b))

# counts prior any condition cell counts from each stage
whole_c <-data %>%
  group_by(group) %>% 
  tally()
# cell counts for module only
data_module <- data[which(data$module > 0.4),]
module_c <-data_module %>%
  group_by(group) %>% 
  tally()
# 2nd condition to count
data_filt <- data_module[which(data_module$panecto > 0),]
filt_c <-data_filt %>%
  group_by(group) %>% 
  tally()
# create data frame with all this information
data_counts <- data.frame(
  stages = whole_c$group,
  total = whole_c$n,
  module = module_c$n,
  filt = filt_c$n,
  percent = paste0(round((filt_c$n/module_c$n)*100,digits=1),"%")
)

