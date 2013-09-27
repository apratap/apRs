library("synapseClient")
library("gdata")
library("plyr")

# #login to synapse
synapseLogin()

###
#get the MsigDB object
###
MSIGDB<-synGet("syn2227979")
load(MSIGDB@filePath) #available as MSigDB R object



###
#get the PCBC samples gene normalized counts
###
syn_geneNormCounts <- synGet('syn1968267')
#read in the file
geneNormCounts <- read.table(syn_geneNormCounts@filePath,header=T,sep='\t')


######
#get the list siginificant genes from comparative analysis in synapse
#####
sigGenes_lists <- readRDS("precomputed_sigGenes_lists.rds")


###
#get the metadata from the synapse for PCBC samples
###
METADATA_ID <- 'syn2024470'
query <- sprintf('select * from entity where parentId=="%s"', METADATA_ID)
metadata <- synQuery(query)
cols_to_be_deleted = c('entity.benefactorId', 'entity.concreteType', 'entity.createdByPrincipalId', 
                       'entity.createdOn', 'entity.createdByPrincipalId', 'entity.id', 
                       'entity.modifiedOn', 'entity.modifiedByPrincipalId', 'entity.noteType', 
                       'entity.versionLabel', 'entity.versionComment', 'entity.versionNumber', 
                       'entity.parentId', 'entity.description', 'entity.eTag')
#delete the unwanted cols
metadata <- metadata[,!names(metadata) %in% cols_to_be_deleted]

#remove the prefix 'entity.' from the df col names
names(metadata) <- gsub('entity.','',names(metadata))


#1. sex
sex <- unique(metadata$donorsex.cell.lines)
sex <- sex[sex != "None"]

#2. group level 1 diff state
level_1_diff_state <- unique(metadata$grouplevel1differentiationstate)
level_1_diff_state <- level_1_diff_state[ !level_1_diff_state %in% c('None','NA')]

#3. group level 3 diff stat
level_3_diff_state <- unique(metadata$grouplevel3differentiationstate)
level_3_diff_state <- level_3_diff_state[ !level_3_diff_state %in% c('None','NA')]

#4. 
cell_origin <- unique(metadata$celltypeoforigin)
cell_origin <- cell_origin[ !cell_origin %in% c('None','NA')]



####
#get the pathway list to show users
####
KEGG_available_pathways = list(
  "KEGG Hematopoietic Cell Lineage"  = "MSigDB$C2.CP.KEGG$KEGG_HEMATOPOIETIC_CELL_LINEAGE",
  "KEGG DNA Replication" = "MSigDB$C2.CP.KEGG$KEGG_DNA_REPLICATION",
  "KEGG Basal Transcription Factors" = "MSigDB$C2.CP.KEGG$KEGG_BASAL_TRANSCRIPTION_FACTORS",
  "custom"          = "custom"
)


###
#sample gene list of the user input area
###
sample_gene_list <- c("GP1BA","GP1BB", "EPO","IL9R", "CD33", "TNF", "GP9", "ITGAM", "CD34",      
                      "CD36", "GP5", "ITGA4", "ITGA3", "KITLG", "ITGA2B")
sample_gene_list <- paste(sample_gene_list,",")




#########
#read the precomputed enriched pathway list
########
df_precomputed_enrichedPathways_in_geneLists = readRDS("precomputed_enrichedPathways_in_geneLists.rds")

df_precomputed_enrichedPathways_in_geneLists$pathways_with_pvalue =  paste(df_precomputed_enrichedPathways_in_geneLists$pathways,
                                                                        '#p.adj_',
                                                                        format.pval(df_precomputed_enrichedPathways_in_geneLists$p.adj,digits=2),
                                                                        sep='')



df_precomputed_enrichedPathways_in_geneLists$significant_gene_list_name 
  
#creating a list of list 
precomputed_enrichedPathways_in_geneLists = split(df_precomputed_enrichedPathways_in_geneLists$pathways_with_pvalue,
                                                  df_precomputed_enrichedPathways_in_geneLists$significant_gene_list_name)




#HACK
#For each geneList add another PATHWAY TYPE "ALL" which indicates use all the pathways for the shiny SERVER/UI
# in this case genes in all the enriched pathways will be shown on the heatmap
precomputed_enrichedPathways_in_geneLists <- lapply(precomputed_enrichedPathways_in_geneLists,function(x) { x[length(x)+1] = 'ALL'; x})




#add place holder for custom gene list
precomputed_enrichedPathways_in_geneLists[['Custom gene list']] = 'NA'





