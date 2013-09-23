library("synapseClient")


# #login to synapse
synapseLogin()

#get the MsigDB object
MSIGDB<-synGet("syn2227979")
load(MSIGDB@filePath) #available as MSigDB R object


#get the PCBC samples gene normalized counts
syn_geneNormCounts <- synGet('syn1968267')
#read in the file
geneNormCounts <- read.table(syn_geneNormCounts@filePath,header=T,sep='\t')



#get the metadata from the synapse for PCBC samples
METADATA_ID <- 'syn2024470'
query <- sprintf('select * from entity where parentId=="%s"', METADATA_ID)
metadata <- synQuery(query)
cols_to_be_deleted = c('entity.benefactorId', 'entity.concreteType', 'entity.createdByPrincipalId', 
                       'entity.createdOn', 'entity.createdByPrincipalId', 'entity.id', 
                       'entity.modifiedOn', 'entity.modifiedByPrincipalId', 'entity.noteType', 
                       'entity.versionLabel', 'entity.versionComment', 'entity.versionNumber', 
                       'entity.parentId', 'entity.description')
#delete the unwanted cols
metadata <- metadata[,!names(metadata) %in% cols_to_be_deleted]

#remove the prefix 'entity.' from the df col names
names(metadata) <- gsub('entity.','',names(metadata))


apply(metadata,2,unique)

#FILTER OPTIONS FOR USERS

#1. sex
sex <- unique(metadata$donorsex.cell.lines)
sex <- sex[sex != "None"]


#2. Origin
origin <- unique(metadata$devorigin)
origin <- origin[origin != "None"]

# 
# diseases <- 
# unique(metadata$grouplevel1differentiationstate)
# metadata$disease

#   disease
# devorigin
# grouplevel1differentiationstate

#get the pathway list to show users
KEGG_available_pathways = list(
                            "KEGG Hematopoietic Cell Lineage"  = "MSigDB$C2.CP.KEGG$KEGG_HEMATOPOIETIC_CELL_LINEAGE",
                            "KEGG DNA Replication" = "MSigDB$C2.CP.KEGG$KEGG_DNA_REPLICATION",
                            "KEGG Basal Transcription Factors" = "MSigDB$C2.CP.KEGG$KEGG_BASAL_TRANSCRIPTION_FACTORS",
                            "custom pathway"          = "custom"
                          )



#sample gene list of the user input area
sample_gene_list <- c("GP1BA", "GP1BB", "EPO","IL9R", "CD33", "TNF", "GP9", "ITGAM", "CD34",      
                      "CD36", "GP5", "ITGA4", "ITGA3", "KITLG", "ITGA2B")
sample_gene_list <- paste(sample_gene_list,",")
