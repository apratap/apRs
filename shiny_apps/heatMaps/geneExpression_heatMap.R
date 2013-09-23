require("Heatplus")


#gene expression heatmap logic


get_geneExpression_heatMap <- function(m,...){
  
  #change the data type to integer
  m <- apply(m,2,as.numeric)
  
  #removing those genes which dont vary much across the samples
  # so any gene with SD < .2 across the samples will be dropped 
  drop_genes <- which(apply(m,1,sd) < .2)
  m <-  m[-drop_genes,]
  
  
  #scaling genes across experiments
  mat.scaled <- t(scale(t(m)))
  
  #custom functions for heatmap plus
  corrdist <- function(x) as.dist( 1 - cor(t(x),method="spearman"))
  hclust.avl = function(x) hclust(x, method="average")
  
  #v3
  reg3 <- regHeatmap(mat.scaled,
                     legend=2,
                     dendrogram = list(clustfun=hclust.avl,
                                       distfun = corrdist))
#                      xlab = "PCBC Samples",
#                      ylab = "Selected genes")
  plot(reg3)
}
  