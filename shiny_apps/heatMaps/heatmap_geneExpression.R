######################
#get the normalized gene counts from synapse
######################

#load the modules
library(synapseClient)
require(Heatplus)

#login to synapse
synapseLogin()

#get the PCBC gene normlz counts from synapse
syn_geneNormCounts <- synGet('syn1968267')

#get the metadata from the synapse about the above sampled
syn_metadata <- synGet('syn2024470')


#get the pathways from MSigDB
#get the MsigDB object
MSIGDB<-synGet("syn2227979")
load(MSIGDB@filePath) #available as MSigDB R object



#read in the file
geneNormCounts <- read.table(syn_geneNormCounts@filePath,header=T,sep='\t')

dim(geneNormCount)
head(geneNormCounts)
names(geneNormCounts)

#select genes in a given pathway
select_pathwayGenes <- MSigDB$C2.CP.KEGG$KEGG_HEMATOPOIETIC_CELL_LINEAGE
selected_geneNormCounts <- subset(geneNormCounts, symbol %in% select_pathwayGenes)


# get only the normalized gene counts
# eliminate the first 3 cols to get rid of the annotation
m <- as.matrix(selected_geneNormCounts[4:ncol(selected_geneNormCounts)])

#change the data type to integer
m <- apply(m,2,as.numeric)


#removing those genes which dont vary much across the samples
# so any gene with SD < .2 across the samples will be dropped 
drop_genes <- which(apply(m,1,sd) < .2)
m <-  m[-drop_genes,]

#log transform
m <- apply(m,1,log2)


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
plot(reg3)





#v1
# reg1 <- regHeatmap(mat.scaled)
# plot(reg1)

#v2
# reg2 <- regHeatmap(mat.scaled,legend=2,col=heat.colors,breaks=-3:3)
# plot(reg2)



#removing those genes which dont vary much across the samples
# so any gene with SD < .2 across the samples will be dropped 
#drop_genes <- which(apply(m,1,sd) < .2)
#m <-  m[-drop_genes,]

#NO more need to do this as the above step will take care of removing genes with SD = 0 or all with expr value = 0
#ignore those rows which have all the values == 0
# these are genes with no expression
#genes_with_no_expression_in_any_sample <- apply(m,1, function(x) all(x==0))
#m <- m[!(genes_with_no_expression_in_any_sample),]


#sampling for testing
#sampled_rows <- sample(seq(1:nrow(m)),200,replace = F)
#sampled_m <- m[sampled_rows,]

#get the HUGO gene names for the sampled genes
sampled_genes_names  <- gene_annotation[sampled_rows,][,2]






MSigDB$C2.CP.KEGG$KEGG_HEMATOPOIETIC_CELL_LINEAGE
