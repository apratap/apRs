#require("Heatplus")
require("pheatmap")
#library("R.cache")
#source("pheatmap_mod.R")

#gene expression heatmap logic
get_geneExpression_heatMap <- function(m,annotation,...){
  
  #change the data type to integer
  m <- apply(m,2,as.numeric)
  
  #removing those genes which dont vary much across the samples
  # so any gene with SD < .2 across the samples will be dropped 
  drop_genes <- which(apply(m,1,sd) < .2)
  m <-  m[-drop_genes,]
  
  #scaling genes across experiments
  mat.scaled <- t(scale(t(m)))
  
  
  pheatmap(mat.scaled,
           scale="none",
           annotation = annotation,
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           clustering_method = "average",
           fontsize_col = 4,
           border_color = NA
  )
  

}


#   memoizedCall(what=pheatmap,
#                mat.scaled,
#                scale="none",
#                annotation = annotation,
#                clustering_distance_rows = "correlation",
#                clustering_distance_cols = "correlation",
#                clustering_method = "average",
#                fontsize_col = 4,
#                border_color = NA,
#                envir=parent.frame(),
#                force=FALSE,
#                sources=NULL,
#                dirs=NULL,
#                verbose=FALSE
#                )
  


  
  #custom functions for heatmap plus
#   corrdist <- function(x) as.dist( 1 - cor(t(x),method="spearman"))
#   hclust.avl = function(x) hclust(x, method="average")
#   
  #v3
#   reg3 <- regHeatmap(mat.scaled,
#                      dendrogram = list(clustfun=hclust.avl,
#                                        distfun = corrdist),
#                      legend = FALSE
#                      )
#   
#   par(mar=c(2,1,1,.5)+0.1,
#       fin=c(10,10),
#       pin=c(8,8))
#   print(paste0('Margin',par()$mar))
#   print(paste0('MGP', par()$mgp))
#   print(paste0('FIN', par()$fin))
#   print(paste0('PIN',par()$pin))
#   plot(reg3)
  