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
#get the meta about PCBC samples
METADATA_ID <- 'syn2024470'
query <- sprintf('select * from entity where parentId=="%s"', METADATA_ID)
metadata <- synQuery(query)
dim(metadata)

cols_to_be_deleted = c('entity.benefactorId', 'entity.concreteType', 'entity.createdByPrincipalId', 
                       'entity.createdOn', 'entity.createdByPrincipalId', 'entity.id', 
                       'entity.modifiedOn', 'entity.modifiedByPrincipalId', 'entity.noteType', 
                       'entity.versionLabel', 'entity.versionComment', 'entity.versionNumber', 
                       'entity.parentId', 'entity.description')
metadata <- metadata[,!names(metadata) %in% cols_to_be_deleted]

#remove the prefix 'entity.' from the df col names
names(metadata) <- gsub('entity.','',names(metadata))

head(metadata)


apply(metadata,2,unique)

#get the pathways from MSigDB
#get the MsigDB object
MSIGDB<-synGet("syn2227979")
load(MSIGDB@filePath) #available as MSigDB R object

MSigDB$C2.CP.KEGG$KEGG_HEMATOPOIETIC_CELL_LINEAGE
"GP1BA"     "GP1BB"     "EPO"       "IL9R"      "CD33"      "TNF" 

#read in the file
geneNormCounts <- read.table(syn_geneNormCounts@filePath,header=T,sep='\t')


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


test_string <- "C,C,C,C C"
unlist(strsplit(test_string,split=c('[\\s,\\n\\r)]'),perl=T))


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
