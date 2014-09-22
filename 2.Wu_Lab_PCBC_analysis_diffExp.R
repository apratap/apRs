############
#Diff Expression
###########

load("~/projects//PCBC_integrative_analysis//analysis//analysis.RData")

source("http://bioconductor.org/biocLite.R")
library("DESeq2")
library("VennDiagram")
library("ggplot2")
library("samr")
library("gplots")
library("grDevices")
library("RColorBrewer")


##tophat mapping based
tophat_read1_counts_file <- "~/projects//PCBC_integrative_analysis/analysis//read1/tophat_mapped_read1_htseq_counts.txt"
tophat_read2_counts_file <- "~/projects//PCBC_integrative_analysis/analysis//read2/tophat_mapped_read2_htseq_counts.txt"
tophat_pairedread_counts_file <- "~/projects//PCBC_integrative_analysis/analysis/paired/tophat_pairedread_htseq_counts.txt"

star_read1_counts_file <- "~/projects//PCBC_integrative_analysis/analysis//read1/star_mapped_read1_htseq_counts.txt"
star_read2_counts_file <- "~/projects//PCBC_integrative_analysis/analysis//read2/star_mapped_read2_htseq_counts.txt"
star_pairedread_counts_file <- "~/projects//PCBC_integrative_analysis/analysis/paired/star_pairedread_htseq_counts.txt"



read_counts <- function(counts_file){
  counts <- read.table(counts_file,header=T)
  rownames(counts) <- counts$gene
  counts$gene <- NULL
  counts$TRA00010815 <- NULL
  counts
}


diffExp <- function(counts_file,min_padj=0.10){
  counts = read_counts(counts_file)
  studyDesign <- data.frame(row.names = colnames(counts),
                            condition = c("ablated","ablated","control","ablated","ablated",
                                          "control","control"))
  #changing the level so that log2 ratio makes more sense and are w.r.t control
  studyDesign$condition <- factor(studyDesign$condition, level=c("control", "ablated"))
  
  dds <- DESeqDataSetFromMatrix(countData=counts,
                                colData = studyDesign,
                                design = ~ condition)
  cat('Running DESeq2\n')
  #DESeq2 based differential expression
  dds <- DESeq(dds)
  res <- as.data.frame(results(dds))
  
  #filter out genes with padj == NA
  res <- res[complete.cases(res),]
  diffExp_genes <- res[res$padj < min_padj, ]
  write.table(res,file=paste0(counts_file,'.difExp_results.txt'),quote=F, sep="\t")
  write.table(diffExp_genes,file=paste0(counts_file,'.difExp_genes_pvalue_lt.10.tsv'),quote=F, sep="\t")
  
  
  #Volcano plot
  cat('Volcano plot \n')
  res['significance'] = 'not-sig'
  res$significance[res$padj < min_padj] =  paste('padj <',min_padj)
  ggplot(data=res,aes(x=log2FoldChange, y=-log(padj), color=significance)) + geom_jitter(size=.8)
  ggsave(paste0(counts_file,'.volcano_plot.png'),width=8,height=8,units="in")
  
  #MA plot
  cat('plotting the MA plot\n')
  png(paste0(counts_file,'.DE_MA_plot.png'))
  plotMA(dds)
  dev.off()
  
  #regularized log transformation
  cat('creating the regularized log2 counts\n')
  rld <- rlogTransformation(dds, blind=TRUE)
  write.table(assay(rld),file=paste0(counts_file,'.regularizedLog_counts.txt'), quote=F, sep="\t")
  #variance stabilized transformation
  cat('creating the variance stabilized counts\n')
  vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
  write.table(assay(rld),file=paste0(counts_file,'.varstabilized_counts.txt'), quote=F, sep="\t")
  
  #PCA plot
  cat('PCA plot \n')
  #save the lattice plot
  png(paste0(counts_file,'.PCA_plot.png'),width=8,height=6,units="in",res=150)
  names <- rownames(colData(vsd))
  pcaplot <- plotPCA(vsd,ntop=1000)
  pcaplot = update(pcaplot, panel = function(x, y, ...) {
    lattice::panel.xyplot(x, y, ...);
    lattice::ltext(x=x, y=y, labels=names, pos=1, offset=1, cex=0.6)
  })
  print(pcaplot)
  dev.off()
  
  return(res[res$padj < min_padj,])
}


#tophat read 1 counts
tophat_read1_de_res_flt <- diffExp(tophat_read1_counts_file)
#tophat read 2 counts
tophat_read2_de_res_flt <- diffExp(tophat_read2_counts_file)
#tophat pairedread counts
tophat_pairedread_de_res_flt <- diffExp(tophat_pairedread_counts_file)

venn.diagram(list(tophat_read_1 = row.names(tophat_read1_de_res_flt), 
                  tophat_read_2 = row.names(tophat_read2_de_res_flt),
                  tophat_pairedreads = row.names(tophat_pairedread_de_res_flt)),
             fill = c("#377EB8", "#E41A1C", 'lightgreen'),
             alpha = c(0.5, 0.5,0.5),
             cex = 2,cat.fontface = 4,lty =3,
             filename="/Volumes//home/apratap/projects//PCBC_integrative_analysis/analysis/tophat_result_comparisons.tiff",
             height=4,width=6,units="in",res=80)



################
##STAR Mapped
################
#star read 1
star_read1_de_res_flt <- diffExp(star_read1_counts_file)
#star read 2
star_read2_de_res_flt <- diffExp(star_read2_counts_file)
#star pairedreads
star_pairedreads_de_res_flt <- diffExp(star_pairedread_counts_file)

venn.diagram(list(star_read_1 = row.names(star_read1_de_res_flt), 
                  star_read_2 = row.names(star_read2_de_res_flt),
                  star_pairedreads=row.names(star_pairedreads_de_res_flt)),
             fill = c("#377EB8", "#E41A1C","lightgreen"),
             alpha = c(0.5, 0.5,0.5),
             cex = 2,cat.fontface = 4,lty =3,
             filename="/Volumes//home/apratap/projects//PCBC_integrative_analysis/analysis/star_result_comparisons.tiff",
             height=4,width=6,units="in",res=80)

save.image(file="~/projects//PCBC_integrative_analysis/analysis/analysis.RData")

#compare tophat and star read1 mapped differentially exp genes
venn.diagram(list(star_read_1 = row.names(star_read1_de_res_flt), 
                  tophat_read_1 = row.names(tophat_read1_de_res_flt)),
                  fill = c("#377EB8", "#E41A1C"),
             alpha = c(0.5, 0.5),
             cex = 2,cat.fontface = 4,lty =3,
             filename="/Volumes//home/apratap/projects//PCBC_integrative_analysis/analysis/star_vs_tophat_read1_comparisons.tiff",
             height=4,width=6,units="in",res=80)

save.image(file="~/projects//PCBC_integrative_analysis/analysis/analysis.RData")



tophat_read1_counts <- read_counts(tophat_read1_counts_file)
star_read1_counts <- read_counts(star_read1_counts_file)
tophat_read2_counts <- read_counts(tophat_read2_counts_file)
star_read2_counts <- read_counts(star_read2_counts_file)
tophat_pairedread_counts <- read_counts(tophat_pairedread_counts_file)
star_pairedread_counts <- read_counts(star_pairedread_counts_file)
R_pairedread_counts <- read_counts("~/projects//PCBC_integrative_analysis/analysis//paired/pairedread_Rbased_counts.tsv")


##heatmap of DiffEXP genes using STAR mapped read1
norm_read_counts <- read.table("~/projects//PCBC_integrative_analysis/analysis//read1//star_mapped_read1_htseq_counts.txt.varstabilized_counts.txt")

diffExp_genes <- norm_read_counts[rownames(star_read1_de_res_flt),]

temp_x <- lapply(rownames(star_read1_de_res_flt), function(x){cat(paste0('|`',x,'`|'),"\n")})
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

heatmap.2( as.matrix(diffExp_genes),
           ,col= myPalette(100)
           ,key=TRUE
           ,keysize=1
           ,density.info="none"
           ,trace = "none"
           ,scale="none"
           ,margins=c(9,9)
)




###############
#Non paramteric analysis with SAMR (Significance Analysis of Microarrays)
###############

samples = c('TRA00010811' = 2,
            'TRA00010812' = 2,
            'TRA00010813' = 2,
            'TRA00010814' = 2,
            'TRA00010815' = 1,
            'TRA00010816' = 1,
            'TRA00010817' = 1,
            'TRA00010818' = 1)


sample_names = sapply(colnames(star_read1_counts), function(x){samples[[x]] })

samr_out <- SAMseq(star_read1_counts, sample_names, 
              geneid=rownames(star_read1_counts),
              genenames=rownames(star_read1_counts),
              resp.type='Two class unpaired',
              fdr.output=.10)
#names(samr_out)
#names(samr_out[['siggenes.table']])
#head(samr_out[['samr.obj']])

genes_up <- samr_out[['siggenes.table']][['genes.up']]
genes_down <- samr_out[['siggenes.table']][['genes.lo']]

write.table(genes_up,file=paste0(star_read1_counts_file,'nonparametric_diffExp_genes_padj_lt.10.tsv'), quote=F, sep="\t", row.names=F)

head(genes_up,n=100)
cat(genes_up[,1])