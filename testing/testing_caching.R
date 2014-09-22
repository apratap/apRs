library(R.cache)

for (kk in 1:5){
  
  cat(sprintf('Iteration #%d\n', kk))
  res <- evalWithMemoization({
    cat("Evaluating epxression...")
    a <- 1
    b <- 2
    c <- 4
    Sys.sleep(1)
    cat("done\n")
    b
  })
  print(res)
  
  stopifnot(a == 1 && b == 2 && c ==4)
  
  #cleanup
  rm(a,b,c)
}

memoizedCall(
            what = plot,
            x = rnorm(100),
            y = rnorm(100))

library("memoise")

a <- function(x) runif(1)
replicate(10,a())
b <- memoise(a)
replicate(10,b())

c <- memoise(function(x) {Sys.sleep(1); runif(1)})
system.time(print(c()))
system.time(print(c()))
forget(c)


#create a matrix of dummy values in R
m1 <- matrix(rnorm(1000,mean=5,sd=4),ncol=20)
m2 <- matrix(rnorm(1000,mean=30,sd=3),ncol=20)
m <- cbind(m1,m2)
rownames(m) <- paste0('row_',seq(nrow(m)))
colnames(m) <- paste0('col_',seq(ncol(m)))
dim(m)

test_rownames <- paste0('testrow_',seq(nrow(m)))

system.time(hcr <- hclust(dist(m)))
system.time(hcc <- hclust(dist(t(m))))
system.time(ddr <- as.dendrogram(hcr))
system.time(ddc <- as.dendrogram(hcc))

#basic
png('test_heatmap.png')
system.time(heatmap(m))
dev.off()


#with useRaster
#Ref: http://blog.revolutionanalytics.com/2011/07/paul-murrell-on-incorporating-images-in-r-charts.html
png('test_heatmap.png')
system.time(heatmap(m,useRaster=TRUE))
dev.off()

#with caching
png('test_heatmap.png')
system.time(heatmap(m,Rowv=ddr,Colv=ddc))
dev.off()

#no clustering and useRaster = True
png('test_heatmap.png')
system.time(heatmap(m,Rowv=ddr,Colv=ddc,useRaster=TRUE ))
dev.off()


source("~/dev/PCBC_HeatMap_app/memoised_pheatmap.R")
source("~/dev/PCBC_HeatMap_app/geneExpression_heatMap.R")

#first call : without caching and drawing Row dendrogram
png('test_heatmap.png')
rownames(m) <- NULL
system.time(memoised_pheatmap(m))
dev.off()


system.time(get_geneExpression_heatMap(m,clustering_distance_rows="euclidean",clustering_distance_cols="euclidean"))
system.time(get_geneExpression_heatMap(m,clustering_distance_rows="correlation",clustering_distance_cols="correlation"))
system.time(get_geneExpression_heatMap(m,clustering_distance_rows="manhattan",clustering_distance_cols="manhattan"))



#second call : should use cached clustering results
png('test_heatmap.png')
system.time(memoised_pheatmap(m))
dev.off()

#third call  : should use cached clustering results and dont draw Row dendogram
png('test_heatmap.png')
system.time(memoised_pheatmap(m,drawRowD=F))
dev.off()

library(grid)
memoised_pheatmap

