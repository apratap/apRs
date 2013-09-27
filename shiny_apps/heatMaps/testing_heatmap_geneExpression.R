######################
#get the normalized gene counts from synapse
######################

#load the modules
library(synapseClient)
require(Heatplus)
library(gplots)
library(RColorBrewer)
library(pheatmap)

#login to synapse
synapseLogin()

#get the PCBC gene normlz counts from synapse
syn_geneNormCounts_v1 <- synGet('syn1968267',version=1)
syn_geneNormCounts_v2 <- synGet('syn1968267',version=2)
syn_geneNormCounts_v3 <- synGet('syn1968267',version=3)



#read in the file
geneNormCounts <- read.table(syn_geneNormCounts_v3@filePath,header=T,sep='\t')

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
                       'entity.parentId', 'entity.description', 'entity.eTag')
metadata <- metadata[,!names(metadata) %in% cols_to_be_deleted]

#remove the prefix 'entity.' from the df col names
names(metadata) <- gsub('entity.','',names(metadata))


apply(metadata,2,unique)

#get the pathways from MSigDB
#get the MsigDB object
MSIGDB<-synGet("syn2227979")
load(MSIGDB@filePath) #available as MSigDB R object



#select genes in a given pathway
select_pathwayGenes <- MSigDB$C2.CP.KEGG$KEGG_HEMATOPOIETIC_CELL_LINEAGE
selected_geneNormCounts <- subset(geneNormCounts, symbol %in% select_pathwayGenes)
dim(selected_geneNormCounts)

#filter based on samples selected
dim(metadata)
metadata$bamName <- gsub('-','.',metadata$bamName) #fix R reading error '-' is converted to .
filtered_sample_names <- intersect(names(selected_geneNormCounts),metadata$bamName)
selected_geneNormCounts <- selected_geneNormCounts[ , names(selected_geneNormCounts) %in% filtered_sample_names] 


selected_metadata <- subset(metadata, bamName %in% filtered_sample_names)
dim(selected_metadata)


annotation <- data.frame(sex = selected_metadata$donorsex.cell.lines,
                         level_1_diff_state = selected_metadata$grouplevel1differentiationstate,
                         level_3_diff_state = selected_metadata$grouplevel3differentiationstate)

#assign the sample names to row names so that the heatmap function could use them for labelling
rownames(annotation) <- selected_metadata$bamName
head(annotation)

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


pheatmap(mat.scaled,
         annotation = annotation,
         scale="none",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "average",
         fontsize_col = 5
)




#custom functions for heatmap plus
corrdist <- function(x) as.dist( 1 - cor(t(x),method="spearman"))
hclust.avl = function(x) hclust(x, method="average")


heatmap.plus(mat.scaled,
             hclustfun=hclust.avl,
             distfun=corrdist,
             scale="none",
             dendrogram="both",
             col=colorRampPalette(c("darkgreen","red","darkred"))(n = 100)
             )


source("heatmap.3.R")
heatmap.3(mat.scaled,hclustfun=hclust.avl,
          distfun=corrdist,na.rm=TRUE,
          scale="none",
          dendrogram="both",
          margins = c(4,6),
          key=TRUE,
          symbreaks=FALSE,
          symkey=FALSE,
          density.info ="none",
          trace = "none",
          labCol = FALSE,
          cexRow=1,
          col=colorRampPalette(c("darkgreen","darkred","red"))(n = 1000))

)



reg3 <- regHeatmap(mat.scaled,
                   dendrogram = list(clustfun=hclust.avl,
                                     distfun = corrdist),
                   legend=2
                   )

plot(reg3)

