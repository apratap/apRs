source("http://bioconductor.org/biocLite.R")
library("Biobase")
library("pasilla")
library("DESeq2")
data("pasillaGenes")
countData <- counts(pasillaGenes)
colData <- pData(pasillaGenes)[,c("condition","type")]
colData


dds <- DESeqDataSetFromMatrix(countData=countData, colData= colData, design = ~ condition)
colData(dds)$condition <- factor(colData(dds)$condition, levels=c("untreated","treated"))

?plotPCA

library("VennDiagram")
install.packages("VennDiagram")
venn.plot <- venn.diagram(list(A = 1:150, B = 121:170), filename="Venn_2set_simple.tiff",
                          fill = c("#377EB8", "#E41A1C"),
                          alpha = c(0.5, 0.5),
                          cex = 2,cat.fontface = 4,lty =3, fontfamily=3

venn.diagram
