source("gene_exp_analysis.R")
library(plyr)
library(reshape)
library(ggplot2)

gene_count_file = "~/apratap_bt/projects/PCBC_integrative_analysis/data/summarized_expression_calls.tsv"
#read the gene counts
gene_counts <- read.table( gene_count_file, sep="\t",header=TRUE )
#ignore a sample with lowly expressed reads
gene_counts$TRA00010815 <- NULL
#get only the counts
m <- as.matrix(gene_counts[2:ncol(gene_counts)])
rownames(m) <- as.character(gene_counts$tracking_id)
#genes which have 0 epxression > 20% of samples
genes_not_uniformly_expressed <- which( apply(m,1,function(x) {sum(x==0)/ncol(m)}) > .20)
m <- m[-genes_not_uniformly_expressed,]
#keep only those genes which have > .2 variation across the samples
drop_genes <- which(apply(m,1,sd) < .2)
m <-  m[-drop_genes,]
#log transform
m <- log2(m+1)

p <- apply(m,1,function(x) t.test(x[1:4],x[5:7])$p.value )
m_filt <- m[ (p < .05),]

write.table(m_filt,"~/apratap_bt/projects/PCBC_integrative_analysis/data/filtered_summarized_exp_values.tsv",
            sep="\t",row.names=T,
            col.names=T)


#correlation
source("gene_exp_analysis.R")
pairs(m_filt, upper.panel=panel.correlation, lower.panel=panel.smooth, main = title )
#heatmap
draw_heatmap(m_filt)
p.adj <- p.adjust(p,method="fdr")
m_filt_p.adj <- m[(p.adj<.05),]


dim(m_filt_p.adj)
#correlation
pairs(m_filt_p.adj, upper.panel=panel.correlation, lower.panel=panel.smooth, main = title)
#heatmap
draw_heatmap(m_filt_p.adj)
m_filt_p.adj


#expression distribution
m_df <- as.data.frame(m)
m_df$genes <- rownames(m_df)
rownames(m_df) <- NULL
m_df_melted <- melt(m_df,id=c('genes'))



################
#DESeq analysis
#################
library("DESeq2")

#setup the parallel environment
options(mc.cores=8)

#get the annotation
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
mm10_txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
showMethods(class=class(mm10_txdb))
mm10_refGenes_by_exons <- exonsBy(mm10_txdb,by="gene")


#get a bam file list for counting
data_base_dir <- "~/projects/PCBC_integrative_analysis/data/"
bamFilesList <- dir(data_base_dir,pattern="accepted_hits.bam",recursive=T)
bamFilesList_obj <- BamFileList(bam_files_list,yieldSize=2000000,index=character())


#restrict the counting to primary alignment of paired reads
flag <- scanBamFlag(isNotPrimaryRead=FALSE,isProperPair=TRUE)
#selecting the params in the bam needed for counting
bam_params <- c("rname","pos","cigar")
params <- ScanBamParam(what=bam_params,flag=flag)


#do the counting
counts <- summarizeOverlaps(features=mm10_refGenes_by_exons,
                            reads = bamFileList_obj,
                            ignore.strand=TRUE,
                            single.end = FALSE,
                            params = params)

#change the colnames (shorten)
new_colnames <- gsub("^.*(TRA.+?)/.+$", '\\1',colnames(counts),perl=T)
colnames(counts) <- new_colnames
save.image(file="~/projects/PCBC_integrative_analysis/data/diff_exp_analysis/R.data")



##load the saved R objects
load("~/apratap_bt/projects/PCBC_integrative_analysis/data/diff_exp_analysis/R.data")

counts_mat <- assays(counts)$counts
counts_mat <- (counts_mat[,-3]) #remove the sample with low read counts #TRA00010815
studyDesign <- data.frame(row.names = colnames(counts_mat),
                          condition = c("ablated","ablated","control","ablated","ablated",
                                        "control","control"),
                          libType   = rep('paired-end',ncol(counts_mat))
)


#DESeq 1 based
#calc the size factor for norm
cds <- newCountDataSet(counts_mat,studyDesign$condition)
cds <- estimateSizeFactors(cds)
sizeFactors(cds)
#normalize counts
head(counts(cds,normalized=TRUE))
#estimate dispersion
cds <- estimateDispersions(cds)
plotDispEsts(cds) # manual inspection
#test for differential expression
diff_exp_test <- nbinomTest(cds,"ablated","control")
plotMA(diff_exp_test)
hist(diff_exp_test$pval)
hist(diff_exp_test$padj)
sum(diff_exp_test$padj < .10, na.rm=T)



##########
#DESeq2 based
##########
#relevelling : refer the DESeq2 vignette for more details
studyDesign$condition <- relevel(studyDesign$condition,ref="control")


diff_exp_dataset  <- DESeqDataSetFromMatrix(countData = counts_mat,
                                         colData = studyDesign,
                                         design = ~ condition)
diff_exp_result <- DESeq(diff_exp_dataset)
res <- as.data.frame(results(diff_exp_result))
res$entrez_gene_id <- rownames(res)



#merge annotation
library("org.Mm.eg.db")
k <- keys(org.Mm.eg.db,keytype="ENTREZID")
keytypes(org.Mm.eg.db)
mm10_gene_names <- select(org.Mm.eg.db, keys=k, cols=c("SYMBOL","GENENAME"), keytype="ENTREZID")
dim(mm10_gene_names)
res <- merge(res,mm10_gene_names,by.x='entrez_gene_id', by.y="ENTREZID")

#order by p.adj
res <- res[order(res$padj),]


#get the diff exp genes with p.adj < .10
sig_genes <- subset(res, padj < .10)
sig_genes

#generate a table of sig genes at 
for_syn_table <- sig_genes[c('SYMBOL','GENENAME','log2FoldChange','padj')]
 x <- apply(for_syn_table,1,function(x) cat(paste(x,collapse='|'),"\n"))


#heatmap for sig genes
#getnorm matrix for heatmap
diff_exp_dataset <- estimateSizeFactors(diff_exp_dataset)
counts_norm_mat <- counts(diff_exp_dataset,normalized=TRUE)
sigGenes_counts_norm_mat <- counts_norm_mat[ sig_genes$entrez_gene_id, ]
#change the rownames of matrix to display gene names on heatmap
new_row_names <- unlist(lapply(rownames(sigGenes_counts_norm_mat), function(x) sig_genes[sig_genes$entrez_gene_id == x,]$SYMBOL))
rownames(sigGenes_counts_norm_mat) <- new_row_names


#heatmap col
source("gene_exp_analysis.R")
png(filename="~/apratap_bt/projects/PCBC_integrative_analysis/data/diff_exp_analysis/DESeq2_sigGenes.png",
     width=10,
     height=12,
     units="in",
     res = 600
     )
draw_heatmap(log2(sigGenes_counts_norm_mat+1),scale=T,labRow=T,key=T)
dev.off()


paste(sig_genes$SYMBOL,collapse=",")

######
#cuffdiff results
#####
install.packages("cummeRbund")

cuffdiff_significant_genes <-c("1500011K16Rik 1600014K23Rik 1700113A16Rik 1810043H04Rik 6330416G13Rik A2m Actr3b Adamts5 Adm Adm2 Ahsa1 Aloxe3 Ampd3 Amph Anxa4 Atf4 Atf6b Atp6v0d2 Atrip Bach1 Bag5 Baiap2l1 Bak1 Bbs5 Bcar1 Bhlhe40 Bnip3 Bola1 Btaf1 Car12 Car14 Ccdc90a Ccl3 Ccng2 Cdc42ep4 Cdkn1b Cdkn2aipnl Chst15 Cited1 Cmtm3 Col27a1 Cox4i2 Cpox Creb5 Crim1 Ctgf Ctla2a Ddit4 Ddt Dixdc1 Dnajb2 Dpp6 Dusp4 Dync1li1 Dynll2 Ear1 Edn1 Efnb3 Egln3 Egr1 Eif4h Emb Emp2 Endov Eno2 Epor Erc2 Eri2 Erich1 Errfi1 F3 Fam162a Fam189b Fam71f2 Fbxl14 Folr2 Fosl2 Fst Gabarapl1 Galk1 Gatsl3 Gtf2h2 Gtpbp6 Gzmd H2-Aa Hand2 Hephl1 Hist1h4j Hk2 Hlx Hsd17b2 Hsd3b6 Htra1 Ier3 Ift43 Igfbp3 Igsf11 Il17rd Il1rl2 Inf2 Ints1 Ipo13 Ippk Irak3 Jph1 Jpx Kansl2 Katna1 Kbtbd11 Kit Klhdc8b Klhl5 Krt18 Krt19 Krt8 Lin28a Lipg Lpcat1 Ly6c2 Ly75 Maff Map3k11 Mark1 Mcat Mcm3 Mcm5 Mcph1 Med9 Mex3c Mfap5 Mgat4b Mir208a Mir677 Mmp15 Mpl Msi1 Mt1 Mt2 Mta2 Mtrf1l Mttp Nabp1 Ndrg1 Nf2 Nfatc4 Nrn1 Nubp2 P4ha1 P4ha2 Pde10a Pdxp Pf4 Pfkfb3 Phactr2 Plekha6 Pln Polr2f Ppp1r13l Ppp1r3c Ppp5c Prdm1 Prl3d2 Prl3d3 Prl7d1 Prl8a2 Prl8a9 Purg Rangap1 Rapgef2 Rasgrp1 Reck Rgs6 Rhox6 Rhox9 Ric3 Rimklb Rmrp Rnaseh2c Rnf103 Rtel1 Runx1 Rxfp1 S100a6 Sec61a1 Selenbp1 Sema3c Sema6d Serpinb9b Serpinb9f Serpinb9g Serpine1 Sftpd Sgk1 Shisa3 Slc2a3 Slc39a14 Slc52a2 Slc6a2 Slc7a3 Slc7a5 Slco2a1 Sprr2g Sprr2h Srgn Srsf12 Star Tfpt Tmem242 Tmem41a Tnfrsf1b Tns3 Trabd Trappc6a Trib3 Ttc8 Ttll7 Usp53 Vegfa Vkorc1 Vldlr Vta1 Wfdc2 Wnk4 Wsb1 Zbtb26 Zc3h4 Zfp395 Zfp398 Zfp41 Zfp414 Zfp691 Zfyve26 Zmat2")
cuffdiff_significant_genes <- unlist(strsplit(cuffdiff_significant_genes, " " ))
names(gene_counts)
m_cuffdff <- subset(gene_counts , tracking_id %in% cuffdiff_significant_genes)
m_cuffdff <- as.matrix(m_cuffdff[2:ncol(m_cuffdff)])

#genes which have 0 epxression > 20% of samples
genes_not_uniformly_expressed <- which( apply(m_cuffdff,1,function(x) {sum(x==0)/ncol(m)}) > .20)
m_cuffdff <- m_cuffdff[-genes_not_uniformly_expressed,]


#correlation
pairs(m_cuffdff, upper.panel=panel.correlation, lower.panel=panel.smooth, main = title)
#heatmap
draw_heatmap(m_cuffdff)
m_filt_p.adj




###################
#testing
###################


biocLite("parathyroidSE")
biocLite("pasilla")
library("pasilla")
library("parathyroidSE")

data("parathyroidGenesSE")
data("pasillaGenes")
countData <- counts(pasillaGenes)
pData(pasillaGenes)
head(countData)
se <- parathyroidGenesSE
counts_mat


ls("package:org.Mm.eg.db")
class(org.Mm.egGENENAME)
columns(org.Mm.eg.db)
keytypes(org.Mm.eg.db)
keys(org.Mm.eg.db,keytype="SYMBOL")


