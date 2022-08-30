# author : Ceren Pajanoja
# date   : September 2022

# Script definition: Analysis for bulk RNAseq in Figure 2, Figure 4, Sup Fig2 and Sup Fig4

library(DESeq2)
library(biomaRt)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(readr)
library(AnnotationHub)
library(ggrepel)
library(tidyr)
library(tidyverse)
library(dplyr)
library(plotly)
library(ggpubr)
library(stringr)
library(Rmisc)

# Read in the dataset
rawData <- read.table("RawCountFile.txt",header=TRUE,row.names=1)  #raw counts here
gData <- read.table("sampletable.txt",header=TRUE)  #sample types here
# Filter: Requiring more than 5 read in at least more than 3 samples
filter <- apply(rawData, 1, function(x) length(x[which(x>5)])>=3)
countsData <- as.data.frame(rawData[filter,])

# Get annotation function  
species<-"chicken"
annotateENSEMBL <- function(ids,species){
  biomart_species<-list(mouse="mmusculus_gene_ensembl",human="hsapiens_gene_ensembl",rat="rnorvegicus_gene_ensembl",dog="cfamiliaris_gene_ensembl",pig="sscrofa_gene_ensembl",fly="dmelanogaster_gene_ensembl",frog="xtropicalis_gene_ensembl",chicken="ggallus_gene_ensembl")
  biomart_organism<-biomart_species[[species]]
  print (paste("using",biomart_organism))
  ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = biomart_organism, host= "http://nov2020.archive.ensembl.org")
  filters = listFilters(ensembl)
  attributes = listAttributes(ensembl)
  ensemblIDs <- data.frame(sapply(strsplit(ids, ".", fixed = TRUE), "[", 1))
  colnames(ensemblIDs)<-"ensembl_gene_id"
  annotation_symbols = getBM(attributes = c("ensembl_gene_id","entrezgene_id","external_gene_name","description"), filters = "ensembl_gene_id", values = ensemblIDs$ensembl_gene_id, mart = ensembl, uniqueRows = TRUE,useCache = FALSE)   
  annotation <- merge (ensemblIDs, annotation_symbols, by.x="ensembl_gene_id", by.y="ensembl_gene_id", all.x = TRUE)
  if (!is.null(annotation)) {
    print ("List succesfully annotated")
  } else {
    print ("There was an error during the annotation")
  }
  return(annotation)
}

# Find DE genes  
ddsHTSeq<-DESeqDataSetFromMatrix(countData=countsData,colData=gData, design=~condition)
ddsHTSeq$condition <- relevel(ddsHTSeq$condition, ref = "HH5")
dds <- DESeq(ddsHTSeq)
res_lrt=results(dds)
res_lrt2 <- res_lrt[order(res_lrt $padj),]
res_lrt_d <- as.data.frame(res_lrt2)
res_lrt_d $ensembl_gene_id <- rownames(res_lrt2) 
res_lrt_annotation <- annotateENSEMBL(res_lrt_d $ensembl_gene_id,species)
res_lrt_final <-merge (res_lrt_d , res_lrt_annotation, by.x="ensembl_gene_id", by.y="ensembl_gene_id", all.x = TRUE)
res_lrt_final$external_gene_name[is.na(res_lrt_final$external_gene_name)] <- as.character(res_lrt_final$ensembl_gene_id[is.na(res_lrt_final$external_gene_name)]) ## added
row.names(res_lrt_final)<- make.unique(res_lrt_final$external_gene_name, sep="-") #
res_lrt_final <- res_lrt_final[order(res_lrt_final$padj),]
res1final <- res_lrt_final
write.table(res1final,file=paste("DESeq2_Results.txt",sep=""),sep="\t",col.names=NA)

## pca by group
png('DE/DESeq2_PCA.png', width = 1000, height = 800, units = "px")
plotPCA(rld, intgroup=c("condition"))
dev.off()
## pca by sample
png("DE/DESeq2_PCA_by_sample.png", width = 1000, height = 800, units = "px")
plotPCA(rld, intgroup=c("sampleName"))
dev.off()
## clustering
hmcol <- colorRampPalette(brewer.pal(9,"Reds"))(50)
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),paste(condition,sampleName, sep=" | "))
png('DE/DESeq2_heatmaps_samplebysample.png', width = 1000, height = 800, units = "px")

#Use results below for input to PANTHER or DAVID:::
res_UP  <- res1final %>%
  filter(pvalue < 50e-3) %>%
  filter(log2FoldChange > 0.75)
write.table(res_UP,file=paste("Res_UP.txt",sep=""),sep="\t",col.names=NA)
res_DOWN  <- res1final %>%
  filter(pvalue < 50e-3) %>%
  filter(log2FoldChange < -0.75)
write.table(res_DOWN,file=paste("Res_DOWN.txt",sep=""),sep="\t",col.names=NA)

# Get VST normalized counts for further plots
vsd <- vst(dds, blind=FALSE)
vsd_mat <- as.data.frame(assay(vsd))
vsd_mat$ensembl_gene_id <- rownames(vsd_mat)
Vsdannotation <- annotateENSEMBL(vsd_mat$ensembl_gene_id,species)
Vsdfinal <-merge (vsd_mat, Vsdannotation, by.x="ensembl_gene_id", by.y="ensembl_gene_id", all.x = TRUE)
Vsdfinal$external_gene_name[is.na(Vsdfinal$external_gene_name)] <- as.character(Vsdfinal$ensembl_gene_id[is.na(Vsdfinal$external_gene_name)]) ## added
row.names(Vsdfinal)<- make.unique(Vsdfinal$external_gene_name, sep="-") # added
write.table(Vsdfinal,file=paste("Vst_counts.txt",sep=""),sep="\t",col.names=NA)
##

selectGene <- c("SOX10","PAX7","NANOG","Pou5f3","KLF4")
dataSub <- Vsdfinal[rownames(Vsdfinal) %in% selectGene, ]
dataSub <- gather(dataSub,"TP","Values",-external_gene_name)
dataSub <- cbind(dataSub,str_split_fixed(dataSub$TP,"_",2))
colnames(dataSub)[4:5]=c("stage","replicate")
df <- dataSub %>% 
  select(-stage) %>%
  mutate(stages = c(rep("HH5",times= dim(dataSub)[1]*4), rep("HH6", times= dim(dataSub)[1]*4), rep("HH7",times= dim(dataSub)[1]*4),
                    rep("2ss",times= dim(dataSub)[1]*4), rep("3ss",times= dim(dataSub)[1]*4), rep("4ss",times= dim(dataSub)[1]*4),
                    rep("5ss",times= dim(dataSub)[1]*4), rep("6ss",times= dim(dataSub)[1]*3), rep("7ss",times= dim(dataSub)[1]*4),
                    rep("8ss",times= dim(dataSub)[1]*3), rep("9ss",times= dim(dataSub)[1]*4), rep("10ss",times= dim(dataSub)[1]*3)))
df=summarySE(df, measurevar="Values", groupvars=c("stages","external_gene_name"))
df <- df[c(46:60,6:45,1:5),] 
df$stages <- factor(df$stages, levels = c("HH5", "HH6", "HH7","1ss","2ss","3ss","4ss","5ss","6ss","7ss","8ss","9ss","10ss"))

# Plot vst counts with std above and below mean
ggplot(df5_c, aes(stages, Values, group = external_gene_name, color = external_gene_name)) +
  geom_line() +
  geom_point() +
  facet_wrap( ~ external_gene_name) +
  labs(title = "Gene expressions", x = "Developmental Stages", y = "Vst Counts") +
  theme(plot.title = element_text(hjust = 0.5, size = 10), strip.text = element_text(size = 15),axis.title.y = element_text(size = 20),axis.title.x = element_text(size = 20),axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14) ) +
  geom_ribbon(aes(ymax = Values + sd, ymin = Values - sd),alpha = 0.5,fill = "grey70",colour=NA)
