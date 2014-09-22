library(synapseClient)

synapseLogin()


#metadata
mRNA_metadata  <- synQuery('select * from entity where parentId == "syn2354333"')
#probe level expression summary file
mRNA_expression_summary  <- synGet('syn2345378')
#local filepath
mRNA_expression_summary@filePath


df1 <- mRNA_metadata[mRNA_metadata$entity.response_class != 'NA',]
length(df1$entity.PatientID)



#metadata
miRNA_metadata  <- synQuery('select * from entity where parentId == "syn2354336"')
#miRNA expression summary file
miRNA_expression_summary  <- synGet('syn2346232')
#local filepath
miRNA_expression_summary@filePath

miRNA_metadata[!duplicated(miRNA_metadata$entity.PatientID),]$entity.PatientID
write.table(miRNA_metadata[!duplicated(miRNA_metadata$entity.PatientID),]$entity.PatientID, file="~/Desktop/miRNA_patients_IDs.tsv", 
            col.names=F, row.names=F,quote=F, sep=",")

df2 <- miRNA_metadata[miRNA_metadata$entity.response_class != 'NA',]
dim(df2)

df2$entity.PatientID
length(df2$entity.PatientID)
length(unique(df2$entity.PatientID))

duplicated(df2$entity.PatientID)

df2$entity.PatientID[duplicated(df2$entity.PatientID)]

df2[duplicated(df2$entity.PatientID),]

df2[df2$entity.PatientID %in%  "5135",]



length(unique(df1$entity.PatientID))
length(unique(df2$entity.PatientID))
