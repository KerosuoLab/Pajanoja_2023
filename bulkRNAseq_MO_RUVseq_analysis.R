# author : Brad Olinger and Ceren Pajanoja
# date   : September 2022

# Script definition: RNAseq Analysis with RUVseq pipeline for Figure 4 and Sup Figure 4

#install biocversion
BiocManager::install("BiocVersion")
#Install if first time
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")
BiocManager::install("RUVSeq")
BiocManager::install("EDASeq")
BiocManager::install("DESeq2")
BiocManager::install("cluster.datasets")
BiocManager::install("devtools")
BiocManager::install("ggbiplot")
BiocManager::install("heatmaply")
#possible this instead
install_github("vqv/ggbiplot")
BiocManager::install("org.Gg.eg.db")
install.packages("dplyr")
install.packages("org.Gg.eg.db")

#Enable Libraries
library(RUVSeq)
library(EDASeq)
library(RColorBrewer)
library(DESeq2)
library(cluster.datasets)
library(devtools)
library(ggbiplot)
library(heatmaply)
library(dplyr)
source("Translate.R")
library("org.Gg.eg.db")
#Weighted scatterplot
library(datasets)
library(ggplot2)
library(rlang)

#Loading Data_____________________________

####Important Values
#percentile value for Id of control genes with CV (coefficient of variance)
#If threshold is 5 counts in 3 samples and you want 1000 HKG
percentile <- 0.05542
ER_Cutoff <- 17033
sam_Count <- 3
count_Count <- 5
# Change this variable to reflect the number of control genes
control_Gene_Count <- 1000
#Change these for applying only the p value threshold (log-fold change limits are optional)
p_Value <- 0.05
lfc_Pos <- 0
lfc_Neg <- 0

for (i in 1){
  #Load and view data
  raw_counts <- read.table("RawCountFile.txt",header=TRUE,row.names=1)
  #Filter Data by expression level
  filter <- apply(raw_counts, 1, function(x) length(x[x>count_Count])>=sam_Count)
  filtered <- raw_counts[filter,]
  #Create S4 class using EDASeq package (for full use of plotting capabilities with ESASeq)
  x <- read.table("groupData.txt",header=TRUE,row.names=1)
  x[,1] <- as.factor(x[,1])
  set <- newSeqExpressionSet(as.matrix(filtered),
                             phenoData = data.frame(x, row.names=colnames(filtered)))
  set_unnormalized <- set
  #Boxplots of relative log expression (RLE = log-ratio of read count to median read count across sample)
  #PCA
  colors <- brewer.pal(6, "Set2")
  colors <- c("dimgray", "dodgerblue3", "gold2", "violetred2")
  plotRLE(set, main="Unnormalized RLE", outline=FALSE, ylim=c(-1.25, 1.25), col=colors[x$Sample], las = 3, ylab = "Relative Log Expression", cex.axis = 1, cex.lab = 1)
  plotPCA(set, col=colors[x$Sample], cex=0.7, theme_size = 5, main = "Unnormalized PCA")
  #Usind EDA seq to normalize the data by upper-quartile normalization
  set <- betweenLaneNormalization(set, which="upper")
  #plotRLE(set, main="Normalized RLE", outline=FALSE, ylim=c(-1.25, 1.25), col=colors[x$Sample], las = 3, ylab = "Relative Log Expression", cex.axis = 1, cex.lab = 1)
  plotPCA(set, main="Normalized PCA", col=colors[x$Sample], cex=1.6)
  ###__________Obtaining List of HKG by Coefficient of Variance
  #transpose for ID of HKG
  allGenes <- t(counts(set_unnormalized))
  #Identify and separate HKG
  #SD of each gene, then AM, then SD/AM = CV
  library(matrixStats)
  counts_sd <- colSds(as.matrix(allGenes))
  counts_am <- colMeans2(as.matrix(allGenes))
  counts_cv <- (counts_sd/counts_am)
  # Add the CV values to the bottom of the original counts file
  # sort in ascending order, ID bottom 2nd percentile
  counts_cv <- rbind(allGenes, counts_cv)
  counts_sorted <- counts_cv[,order(counts_cv[nrow(counts_cv),])]
  counts_2p <- counts_sorted[,counts_sorted[nrow(counts_sorted),] <= quantile(counts_sorted[nrow(counts_sorted),], probs = percentile  )]
  #Remove the bottom row (this is only the CV)
  counts_2p <- counts_2p[-c(nrow(counts_2p)),]
  head(counts_2p)
  #rename for RUVseq
  spikes <- counts_2p
  #show the number of genes that fall below the 2nd percentile CV, write to table (only if you need this)
  ncol(spikes)
  #create a list of all DE genes (those that are not HKG)
  genes <- counts_sorted[,-c(1:ncol(spikes))]
  #remove CV value at th bottom
  genes <- genes[-c(nrow(counts_2p)),]
  #check the number of genes in each category
  ncol(allGenes) #total genes that passed theshhold
  ncol(genes) #total non-HKG
  ncol(spikes) #total HKG
  #Re-assign spikes to variable for creating a set of only control genes
  HKG_CV <- t(spikes)
  #transpose again and extract only names for RUVseq
  spikes <- t(spikes)
  spikes <- rownames(spikes)
  genes <- t(genes)
  genes <- rownames(genes)
  totalGenes <- t(allGenes)
  #__________________END List of HKG by CV_____________

}
###______End Initial Setup________

for (i in 1){
  #RUVs using replicate/negative control samples
  differences <- makeGroups(x$Sample)
  set3 <- RUVs(set, spikes, k=9, differences)
  pca_Title <- paste("K",9, sep="")
  par(mfrow = c(1, 1))
  plotPCA(set3, col=colors[x$Sample], cex=0.7, main = pca_Title)
  plotRLE(set3, outline=FALSE, ylim=c(-0.5, 0.5), col=colors[x$Sample], las = 3, ylab = "Relative Log Expression", cex.axis = 1, cex.lab = 1)

  #DGE Analysis using DESeq2 *after* RUVs
  dds <- DESeqDataSetFromMatrix(countData = counts(set3),
                                colData = pData(set3),
                                design = as.formula(~ W_1 + W_2 + W_3 + W_4
                                                    + W_5 + W_6 + W_7 + W_8
                                                    + W_9 + Sample))

  dds <- DESeq(dds)
  dds <- dds[which(mcols(dds)$betaConv),]
  res_Klf4 <- results(dds, contrast=c("Sample", "KLF_E", "C"))
  res_Nanog <- results(dds, contrast=c("Sample", "Nanog_E", "C"))
  res_PouV <- results(dds, contrast=c("Sample", "PouV_E", "C"))

  plotMA(res_Klf4, ylim=c(-2,2))
  plotMA(res_PouV, ylim=c(-2,2))
  plotMA(res_Nanog, ylim=c(-2,2))

  #Initialize Translation Tables
  load_translation_tables()
  #Saving table of DEG
  #Nanog DEG
  res_Nanog_DEG <- res_Nanog[ which(res_Nanog$log2FoldChange>=lfc_Pos
                                    & res_Nanog$padj < p_Value), ]
  res_Nanog_DEG_2 <- res_Nanog[ which(res_Nanog$log2FoldChange<=lfc_Neg
                                      & res_Nanog$padj < p_Value), ]
  res_Nanog_DEG <- rbind(res_Nanog_DEG, res_Nanog_DEG_2)
  gene_Names <- rownames(res_Nanog_DEG)
  nanog_sym <- translated_genes(gene_Names)
  error_message <- check_error(gene_Names)
  print(paste("Genes with conflicting translations: ",error_message, sep=""))
  res_Nanog_DEG <- cbind(res_Nanog_DEG, nanog_sym)

  #PouV DEG
  res_PouV_DEG <- res_PouV[ which(res_PouV$log2FoldChange>=lfc_Pos
                                  & res_PouV$padj < p_Value), ]
  res_PouV_DEG_2 <- res_PouV[ which(res_PouV$log2FoldChange<=lfc_Neg
                                    & res_PouV$padj < p_Value), ]
  res_PouV_DEG <- rbind(res_PouV_DEG, res_PouV_DEG_2)
  gene_Names <- rownames(res_PouV_DEG)
  pouv_sym <- translated_genes(gene_Names)
  error_message <- check_error(gene_Names)
  print(paste("Genes with conflicting translations: ",error_message, sep=""))
  res_PouV_DEG <- cbind(res_PouV_DEG, pouv_sym)

  #Klf4 DEG
  res_Klf4_DEG <- res_Klf4[ which(res_Klf4$log2FoldChange>=lfc_Pos
                                  & res_Klf4$padj < p_Value), ]
  res_Klf4_DEG_2 <- res_Klf4[ which(res_Klf4$log2FoldChange<=lfc_Neg
                                    & res_Klf4$padj < p_Value), ]
  res_Klf4_DEG <- rbind(res_Klf4_DEG, res_Klf4_DEG_2)
  gene_Names <- rownames(res_Klf4_DEG)
  klf4_sym <- translated_genes(gene_Names)
  error_message <- check_error(gene_Names)
  print(paste("Genes with conflicting translations: ",error_message, sep=""))
  res_Klf4_DEG <- cbind(res_Klf4_DEG, klf4_sym)

  #Continue if you want to add base mean, max sample
  exp_col <- as.numeric(c(3,4,5,6,11,12,13,14,18,19,20,21))
  ctr_col <- as.numeric(c(1,2,7,8,9,10,15,16,17))
  #Ordering and saving to file
  #Klf4
  nrow_Klf4 <- nrow(res_Klf4_DEG)
  res_Klf4_DEG_Final <- NULL
  res_Klf4_DEG_Final <- matrix(, ncol = 7, nrow = nrow_Klf4)
  colnames(res_Klf4_DEG_Final) <- c("Gene Symbol", "padj", "l2fc",
                                    "MaxRawCount", "Raw Mean",
                                    "Exp Mean", "Control Mean")
  #order by padj for easy analysis
  res_Klf4_DEG <- res_Klf4_DEG[order(res_Klf4_DEG$"padj"),]
  rownames(res_Klf4_DEG_Final) <- rownames(res_Klf4_DEG)
  row_Names_Klf4 <- rownames(res_Klf4_DEG_Final[])
  #fill gene symbol, padj, lfc, max raw count, and mean counts
  res_Klf4_DEG_Final[,1] <- res_Klf4_DEG[,7]
  res_Klf4_DEG_Final[,2] <- res_Klf4_DEG[,6]
  res_Klf4_DEG_Final[,3] <- res_Klf4_DEG[,2]
  for (i in 1:nrow(res_Klf4_DEG_Final)) {
    res_Klf4_DEG_Final[i,4] <- max(raw_counts[row_Names_Klf4[i],])
  }
  for (i in 1:nrow(res_Klf4_DEG_Final)) {
    res_Klf4_DEG_Final[i,5] <- rowMeans(raw_counts[row_Names_Klf4[i],])
  }
  for (i in 1:nrow(res_Klf4_DEG_Final)) {
    res_Klf4_DEG_Final[i,6] <- rowMeans(raw_counts[row_Names_Klf4[i],exp_col])
  }
  for (i in 1:nrow(res_Klf4_DEG_Final)) {
    res_Klf4_DEG_Final[i,7] <- rowMeans(raw_counts[row_Names_Klf4[i],ctr_col])
  }

  #PouV
  nrow_PouV <- nrow(res_PouV_DEG)
  res_PouV_DEG_Final <- NULL
  res_PouV_DEG_Final <- matrix(, ncol = 7, nrow = nrow_PouV)
  colnames(res_PouV_DEG_Final) <- c("Gene Symbol", "padj", "l2fc",
                                    "MaxRawCount", "Raw Mean",
                                    "Exp Mean", "Control Mean")
  #order by padj for easy analysis
  res_PouV_DEG <- res_PouV_DEG[order(res_PouV_DEG$"padj"),]
  rownames(res_PouV_DEG_Final) <- rownames(res_PouV_DEG)
  row_Names_PouV <- rownames(res_PouV_DEG_Final[])
  #fill gene symbol, padj, lfc
  res_PouV_DEG_Final[,1] <- res_PouV_DEG[,7]
  res_PouV_DEG_Final[,2] <- res_PouV_DEG[,6]
  res_PouV_DEG_Final[,3] <- res_PouV_DEG[,2]
  for (i in 1:nrow(res_PouV_DEG_Final)) {
    res_PouV_DEG_Final[i,4] <- max(raw_counts[row_Names_PouV[i],])
  }
  for (i in 1:nrow(res_PouV_DEG_Final)) {
    res_PouV_DEG_Final[i,5] <- rowMeans(raw_counts[row_Names_PouV[i],])
  }
  for (i in 1:nrow(res_PouV_DEG_Final)) {
    res_PouV_DEG_Final[i,6] <- rowMeans(raw_counts[row_Names_PouV[i],exp_col])
  }
  for (i in 1:nrow(res_PouV_DEG_Final)) {
    res_PouV_DEG_Final[i,7] <- rowMeans(raw_counts[row_Names_PouV[i],ctr_col])
  }

  #Nanog
  nrow_Nanog <- nrow(res_Nanog_DEG)
  res_Nanog_DEG_Final <- NULL
  res_Nanog_DEG_Final <- matrix(, ncol = 7, nrow = nrow_Nanog)
  colnames(res_Nanog_DEG_Final) <- c("Gene Symbol", "padj", "l2fc",
                                     "MaxRawCount", "Raw Mean",
                                     "Exp Mean", "Control Mean")
  #order by padj for easy analysis
  res_Nanog_DEG <- res_Nanog_DEG[order(res_Nanog_DEG$"padj"),]
  rownames(res_Nanog_DEG_Final) <- rownames(res_Nanog_DEG)
  row_Names_Nanog <- rownames(res_Nanog_DEG_Final[])
  #fill gene symbol, padj, lfc
  res_Nanog_DEG_Final[,1] <- res_Nanog_DEG[,7]
  res_Nanog_DEG_Final[,2] <- res_Nanog_DEG[,6]
  res_Nanog_DEG_Final[,3] <- res_Nanog_DEG[,2]
  for (i in 1:nrow(res_Nanog_DEG_Final)) {
    res_Nanog_DEG_Final[i,4] <- max(raw_counts[row_Names_Nanog[i],])
  }
  for (i in 1:nrow(res_Nanog_DEG_Final)) {
    res_Nanog_DEG_Final[i,5] <- rowMeans(raw_counts[row_Names_Nanog[i],])
  }
  for (i in 1:nrow(res_Nanog_DEG_Final)) {
    res_Nanog_DEG_Final[i,6] <- rowMeans(raw_counts[row_Names_Nanog[i],exp_col])
  }
  for (i in 1:nrow(res_Nanog_DEG_Final)) {
    res_Nanog_DEG_Final[i,7] <- rowMeans(raw_counts[row_Names_Nanog[i],ctr_col])
  }

  #Save Tables
  write.table(res_PouV_DEG_Final, "res_PouV_DEG_Final.txt", sep="\t", col.names = NA, row.names = TRUE)
  write.table(res_Klf4_DEG_Final, "res_Klf4_DEG_Final.txt", sep="\t", col.names = NA, row.names = TRUE)
  write.table(res_Nanog_DEG_Final, "res_Nanog_DEG_Final.txt", sep="\t", col.names = NA, row.names = TRUE)
  #Separate code to normalize data
  vsd <- vst(dds, blind=FALSE)
}

load_translation_tables()
ensembl_2 <- rev_translated_genes(gene_sym_2)
error_message <- rev_check_error(gene_sym_2)
print(paste("Genes with conflicting translations: ",error_message, sep=""))

heatmap2_data <- matrix(, nrow = 90, ncol = 21)
colnames(heatmap2_data) <- c("KC1", "KC3", "KE1", "KE2", "KE3", "KE4",
                             "NC1", "NC2", "NC3", "NC4", "NE1", "NE2",
                             "NE3", "NE4", "PC1", "PC2", "PC3", "PE1",
                             "PE2", "PE3", "PE4")
rownames(heatmap2_data) <- c(gene_sym_2)

#Fill the table
for (i in 1:21){
  for (val in 1:90){
    heatmap2_data[val,i] <- as.numeric(assay(vsd)[ensembl_2[val],i])
  }
}

###Simplifying the heatmap in order to only have 1 column for the controls/Exp groups
heatmap4_data <- matrix(,ncol = 4, nrow = 90)
colnames(heatmap4_data) <- c("Control", "Klf4", "Nanog", "PouV")
rownames(heatmap4_data) <- rownames(heatmap2_data)
#fill the mean control values first
for (i in 1:90){
  heatmap4_data[i,1] <- mean(heatmap2_data[i,1:9])
}
#fill the remaining KD values
for (val in 1:90){
  heatmap4_data[val,2] <- mean(heatmap2_data[val,10:13])
  heatmap4_data[val,3] <- mean(heatmap2_data[val,14:17])
  heatmap4_data[val,4] <- mean(heatmap2_data[val,18:21])
}
heatmaply(heatmap4_data,
          na.rm = TRUE,
          #dendrogram = "row",
          scale = "row",
          #col_side_colors = c(rep(0, 5), rep(1, 5), rep(2,5)),
          #col_side_colors = FP_NF_Collabel,
          NA.value = "grey50",
          #Rowv=FALSE,
          Colv=FALSE,
          xlab = "", ylab = "",
          main = "90 Gene Expression Levels -- All DE in all 3 KD groups",
          #scale = "column",
          margins = c(200,100,40,20),
          grid_color = "white",
          grid_width = 0.00001,
          titleX = FALSE,
          #hide_colorbar = TRUE,
          branches_lwd = 0.1,
          label_names = c("Gene", "KD-Morpholino", "Log Expression"),
          fontsize_row = 8, fontsize_col = 8,
          labCol = colnames(heatmap4_data),
          labRow = rownames(heatmap4_data),
          heatmap_layers = theme(axis.line=element_blank()),
          hide_colorbar = FALSE,
          colorbar_len = 0.3,
          colorbar_thickness = 30,
          colorbar_xanchor = if (row_dend_left) "right" else "left",
          colorbar_yanchor = "bottom",
          colorbar_xpos = if (row_dend_left) -0.1 else 1.1,
          colorbar_ypos = 0,
          file="90_Gene_Expression_by_Embryo_Commonto3kd_B.html",
          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
            low = "blue",
            high = "red",
            midpoint = 0,
            limits = c(-1, 1)
          ))

