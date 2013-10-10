#!/usr/bin/env Rscript

#########################
#By: Abhishek Pratap
#Sage Bionetworks
#######################


#load the modules
library(getopt,quietly=TRUE)
library(gplots,quietly=TRUE)
library("RColorBrewer")


SPEC =  matrix ( 
     	       	 c( 
		    'verbose'		, 'v' , 2, "integer",
		    'help'     		, 'h' , 0, "logical",
		    'gene_count'	, 'g' , 1, "character",
		    'project_name'	, 'p'  , 1, "character"
		    ),ncol=4, byrow=TRUE);

opt = getopt(SPEC);


if ( ! is.null(opt$help)){
   cat("Usage : gene_count_correlation_analysis.R -g <gene_count_file> -p <project_name> \n\n");
   q(status=1)
}



panel.correlation <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs( cor(x, y,method="spearman") )
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}


panel.smooth <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
                          cex = 1, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
    points(x, y, pch = pch, col = col, bg = bg, cex = cex)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
        lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
            col = col.smooth, ...)
}


draw_heatmap <- function(m){

  #keep only the genes which have expression in atleast one of the sample
  genes_with_no_expression_in_any_sample <- apply(m,1,function(x) all(x == 0))
  m <- m[!(genes_with_no_expression_in_any_sample),]
  
  #scaling genes across experiments
  mat.scaled <- t(scale(t(m)))
  
  c <- cor(t(mat.scaled),method="spearman")
  dist.corr <- as.dist(1-c)
  clust.agg.complete <- hclust(dist.corr,method="complete",members=NULL)
  
  # Cuts the tree and creates color vector for clusters.
#   trimmed_tree_clusters <- cutree(clust.agg.complete,h=max(clust.agg.complete$height)/1.5)
#   mycolhc <- rainbow(length(unique(trimmed_tree_clusters)),start=0.1,end=0.9)
#   mycolhc <- mycolhc[as.vector(trimmed_tree_clusters)]
#     myheatcol <- greenred(75)
  
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

  heatmap.2( mat.scaled,Rowv=as.dendrogram(clust.agg.complete)
             ,col= myPalette(100)
             ,key=TRUE
             ,keysize=1
             ,density.info="none"
             ,trace = "none"
             ,labRow = ""
             ,margins=c(9,9)
  )
}



########
# MAIN
########
gene_count_file <- opt$gene_count
project_name   <- opt$project
project_name <- sub("\\s+", "", project_name)

#gene_count_file = "~/apratap_bt/projects/PCBC_integrative_analysis/data/summarized_expression_calls.tsv"
#read the gene counts
gene_counts <- read.table( gene_count_file, sep="\t",header=TRUE )

#get only the counts
m <- as.matrix(gene_counts[2:ncol(gene_counts)])

#genes which have 0 epxression > 20% of samples
genes_not_uniformly_expressed <- which( apply(m,1,function(x) {sum(x==0)/ncol(m)}) > .20)
m <- m[-genes_not_uniformly_expressed,]

#keep only those genes which have > .2 variation across the samples
drop_genes <- which(apply(m,1,sd) < .2)
m <-  m[-drop_genes,]


###########
#draw the correlation
##########
title <- paste("Project_",project_name,"_Gene_Expression_Correlation");
corr_map_file <- paste("Project_",project_name,"_Gene_Expression_Correlation.png",sep="")
png(corr_map_file, width=1300, height=1000)
pairs(m, upper.panel=panel.correlation, lower.panel=panel.smooth, main = title, size=12)
dev.off()


#############
#heat map
#############
heatmap_file <- paste("Project_",project_name,"_Gene_Expression_heatmap.png",sep="")
png(heatmap_file, width=1300, height=1000)
draw_heatmap(m)
dev.off()



