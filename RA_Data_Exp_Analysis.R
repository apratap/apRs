library("plyr")
library("reshape")
library("ggplot2")
library("gdata")
library("plyr")
library("Rsamtools")
library("RColorBrewer")
library("gplots")
library("pheatmap")
library("GGally")

#####
#reading in the clinical data
#####
BASE_DIR = "~/Desktop/"
#first file
RA_clinical_data_file <- paste0(BASE_DIR,"AIR_Clinical_Data.xls")
patient_metadata <- read.xls(RA_clinical_data_file,sheet='demo')
#second file
RA_clinical_data_file_2 <- paste0(BASE_DIR,"AIR_clinical_data_RNA-seq.xls")
patient_metadata_extra_info <- read.xls(RA_clinical_data_file_2)
new_colnames <- c('sample_id','collab_participant_id','collab_sample_id','sample_date','collection_time',
                  'date_received','transit_days','CRP','RF','CCP','dx1_physician_confirmed',
                  'pdx1_patient_reported', 'category')
colnames(patient_metadata_extra_info) <- new_colnames
#remove the first row
patient_metadata_extra_info <- patient_metadata_extra_info[-1,]
#fix the sample id col to match with the patient_metdata qstkidid
patient_metadata_extra_info$sample_id <- sub('-','',patient_metadata_extra_info$sample_id)


#merge with the main patient meta data
cols_to_keep = c('sample_id','collab_participant_id','collab_sample_id','dx1_physician_confirmed',
                 'pdx1_patient_reported', 'category','transit_days')
patient_metadata <- merge(patient_metadata,patient_metadata_extra_info[cols_to_keep],by.x='qstkitid',by.y='sample_id')

#create a new col : disease_activity : low / high
patient_metadata[grepl('*.LOW',patient_metadata$category),'disease_activity'] = 'low'
patient_metadata[grepl('*.HI',patient_metadata$category),'disease_activity'] = 'high'

       
#########
#patient drugs
#########
patient_drugs <- read.xls(RA_clinical_data_file,sheet='drugs')
#change the datetime object
patient_drugs$date <- as.Date(as.character(patient_drugs$date),format="%Y-%m-%d")



#Histogram of C-Reactive Protein values across the patients
ggplot(data=patient_metadata,aes(x=CRP,fill=sex)) + geom_histogram(binwidth=.05) +
  xlab('C-Reactive Protein value') + ylab('#patients') + geom_vline(xintercept=c(.20), linetype="dotted",colour="red")

#Histogram of RF values
ggplot(data=patient_metadata,aes(x=RF,fill=sex)) + geom_histogram(binwidth=10) + xlab('Rheumatoid Factor(RF) clinically measured value') + ylab('#patients') 

#scatter plot RF v/s CRP vals
ggplot(data=patient_metadata,aes(x=RF,y=CRP,colour=sex)) + geom_point() + 
  xlab('Rheumatoid Factor(RF) vals') + ylab('C-Reactive Protein vals') + 
  ggtitle('Clinically measured RF v/s CRP values') + theme_bw()

#merge the drug and patient clinical meta data
df <- merge(patient_metadata,patient_drugs,by="qstkitid")
# remove the date.y column for now
# AND remove the duplicates after that
drugs_per_patient <- unique(subset(df, select=-c(date.y)))
drugs_mat <- as.matrix(subset(drugs_per_patient,select=-c(qstkitid,collectiondate.x,date.x,
                                                      age,sex,CRP,RF,collectiondate.y)))
rownames(drugs_mat) <- drugs_per_patient$qstkitid
#create a annotation data frame for heat map
annotation <- subset(drugs_per_patient,select=c(age,sex,CRP,RF))
rownames(annotation) <- drugs_per_patient$qstkitid

#render the heatmap of drugs taken by the patient
pheatmap(t(drugs_mat),
         scale="none",
         annotation = annotation,
         border_color = NA
)



###########
#COUNTING
###########
#get the annotation
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
hg19_txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
showMethods(class=class(hg19_txdb))
hg19_refGenes_by_exons <- exonsBy(hg19_txdb,by="gene")


#change the chr names to match that of BAMs to enable counting
new_labels <- sub("chr","",seqlevels(hg19_txdb))
new_labels <- toupper( sub("Un_","",new_labels) )
new_labels <- unlist(lapply(new_labels,function(x) if(startsWith(x,'GL')) paste0(x,'.1') else (x) ))
names(new_labels) <- seqlevels(hg19_txdb)
hg19_refGenes_by_exons <-  renameSeqlevels(hg19_refGenes_by_exons,new_labels)

#get a bam file list for counting
data_base_dir <- "~/projects/RA//work"
#data_base_dir <- "~/projects/RA/testing/"
bamFilesList <- dir(data_base_dir,pattern="*.bam$",recursive=T,full.names=T)
bamFilesList_obj <- BamFileList(bamFilesList,yieldSize=2000000,index=character())
#restrict the counting to primary alignment of paired reads
flag <- scanBamFlag(isNotPrimaryRead=FALSE,isProperPair=TRUE)
#selecting the params in the bam needed for counting
bam_params <- c("rname","pos","cigar")
params <- ScanBamParam(what=bam_params,flag=flag)

#do the counting
#setup the parallel environment
options(mc.cores=10)
counts <- summarizeOverlaps(features=hg19_refGenes_by_exons,
                            reads = bamFilesList_obj,
                            ignore.strand=TRUE,
                            single.end = FALSE)
save.image(file="~/projects/RA/work/expAnalysis/R.data")


#temp
###load the saved R objects
load("~/Desktop/R.data")


################
#Post counting processing
################
#change the colnames
cols <- colnames(counts)
new_cols <- sub("^.*work/Data/","",cols)
new_cols <- sub("_.*bam","",new_cols)
colnames(counts) <- new_cols


#debug
new_cols[ ! new_cols %in% patient_metadata$collab_sample_id ]
patient_metadata$collab_sample_id[ ! patient_metadata$collab_sample_id %in%  new_cols]
patient_metadata[ ! patient_metadata$collab_sample_id %in%  new_cols, ][,c('qstkitid','collab_participant_id',
                                                                           "collab_sample_id")]
#END DEBUG


#get the count mat
counts_mat <- assays(counts)$counts

#write read counts
write.table(counts_mat,"~/Desktop/RA_Data_Expression_Counts.txt",sep="\t")

#PUSH to synapse
library(synapseClient)
synapseLogin()
#get list of files used for counting : synapse ids
synBams <- synapseQuery('select id from entity WHERE parentId == "syn2280648" and fileType =="bam" ')
synBams <- as.vector(synBams$entity.id)
counts_file <- File("~/Desktop/RA_Data_Expression_Counts.txt",parentId = 'syn2290931')
counts_activity <- Activity(used=synBams,executed=list(url='https://raw.github.com/apratap/apRs/master/RA_Data_Exp_Analysis.R'))
counts_file <- synStore(counts_file,activity=counts_activity)

#total number of reads mapped in exonic region per sample
total_exonicReads_per_sample <- apply(counts_mat,2,sum)
p <- qplot(names(total_exonicReads_per_sample),total_exonicReads_per_sample,geom=c("line","point"),group="identity") 
p + theme(axis.text.x=element_text(angle=90, hjust=0))


#RPKM
returnRPKM <- function(counts, features) {
  featureLengthsInKB <- sum(width(reduce(features))) / 1000 # geneLength length of exon per gene in Kbp
  millionsMapped <- sum(counts) / 1e+06 #Factor for normalizing , converting to per million mapped reads / feature
  rpm <- counts / millionsMapped  #RPM : reads per kb of exon model
  rpkm <- rpm / featureLengthsInKB #RPKM : reads per kilobase of exon model per million mapped reads
  return(rpkm)            
}


############
# hg19 annotation conversion
############
# library("org.Hs.eg.db")
# k <- keys(org.Hs.eg.db,keytype="ENTREZID")
# keytypes(org.Hs.eg.db)
# mm10_gene_names <- select(org.Hs.eg.db, keys=k, cols=c("SYMBOL","GENENAME"), keytype="ENTREZID")
# dim(mm10_gene_names)
# res <- merge(res,mm10_gene_names,by.x='entrez_gene_id', by.y="ENTREZID")



#################
#Exploratory analysis for diff exp signal
#################
rpkm_counts_mat <- apply(counts_mat,2,function(x) returnRPKM(x,hg19_refGenes_by_exons))
rpkm_counts_mat <- log10(rpkm_counts_mat+1)  #log 2 transform the counts

#temp
bak_rpkm <- rpkm_counts_mat
rpkm_counts_mat <- bak_rpkm


#box plot of expression values across all the samples
melted_rpkm_counts <- melt(rpkm_counts_mat)
colnames(melted_rpkm_counts) <- c('gene_id','sample','value')
ggplot(data=melted_rpkm_counts, aes(x=factor(sample),y=value)) + geom_boxplot() + 
  theme(axis.text.x=element_text(angle=90, hjust=0)) + xlab('Samples') + ylab('Normalized Counts(RPKM, log10 scale)')

#filter genes with 0 exp in >= 20% samples, keep only 
genes_to_keep <- which( apply(rpkm_counts_mat, 1, function(x) {sum(x == 0.0)/length(x)} ) <= .20 )
rpkm_counts_mat <- rpkm_counts_mat[ genes_to_keep, ]
dim(rpkm_counts_mat)

#remove genes whoose expression dont vary much across ALL the samples
#cut off variance = .20
flt_rpkm_counts_mat <- rpkm_counts_mat[which(apply(rpkm_counts_mat,1,sd) > .15 ), ]
dim(flt_rpkm_counts_mat)

#keep only those samples for which we have meta-data
#mostly a sanity check step 
# for now 3 samples dont have metadata
flt_counts_mat <- counts_mat[ , colnames(counts_mat) %in% patient_metadata$collab_sample_id ]
flt_rpkm_counts_mat <- flt_rpkm_counts_mat[, colnames(flt_rpkm_counts_mat) %in% patient_metadata$collab_sample_id ]

#create a annotation data frame for heat map
annotation <- subset(patient_metadata,select=c(disease_activity))
rownames(annotation) <- patient_metadata$collab_sample_id

#scaled the data across experiments
m <- (t(scale(t(flt_rpkm_counts_mat))))
pheatmap(m,
         annotation = annotation,
         scale="none",
         border_color = NA
         )

##########
#PCA Plot
##########
pca <- prcomp(flt_rpkm_counts_mat,scale=T)
summary(pca)

pca$x

melted <- cbind(disease_type = as.vector(conditions),melt(pca$x[,1:9]))
head(melted)
ggplot(data=melted) + geom_bar(aes(x=X1,y=value,fill=disease_type),stat="identity") +
  facet_wrap(~X2)
  





############
#DESeq2 based
############
library("DESeq2")
#creating the study design
conditions <- sapply(colnames(flt_counts_mat), function(x) {
                                                    patient_metadata[patient_metadata$collab_sample_id == x,'disease_activity'] 
                                                  })

studyDesign <- data.frame(condition = as.character(conditions),
                          libType   = rep('paired-end',ncol(flt_counts_mat))
                          )
rownames(studyDesign) <- colnames(flt_counts_mat)

#main entry point
diff_exp_dataset <- DESeqDataSetFromMatrix(countData = flt_counts_mat,
                                           colData = studyDesign,
                                           design = ~ condition)

#relevelling : refer the DESeq2 vignette for more details
colData(diff_exp_dataset)$condition <- factor(colData(diff_exp_dataset)$condition,
                                              levels= c('low','high'))

diff_exp_result <- DESeq(diff_exp_dataset)
res <- as.data.frame(results(diff_exp_result))
res$entrez_gene_id <- rownames(res)



#merge annotation
library("org.Hs.eg.db")
k <- keys(org.Hs.eg.db,keytype="ENTREZID")
keytypes(org.Hs.eg.db)
hg19_gene_names <- select(org.Hs.eg.db, keys=k, columns=c("SYMBOL","GENENAME"), keytype="ENTREZID")
res <- merge(res,hg19_gene_names,by.x='entrez_gene_id', by.y="ENTREZID")

#order by p.adj
res <- res[order(res$padj),]


#get the diff exp genes with p.adj < .10
sig_genes <- subset(res, padj < .10)

#generate a table of sig genes at
head(sig_genes)
for_syn_table <- sig_genes[c('SYMBOL','GENENAME','log2FoldChange','pvalue','padj')]
x <- apply(for_syn_table,1,function(x) cat(paste(x,collapse='|'),"\n"))


#heatmap for sig genes
#getnorm matrix for heatmap
diff_exp_dataset <- estimateSizeFactors(diff_exp_dataset)
counts_norm_mat <- counts(diff_exp_dataset,normalized=TRUE)
sigGenes_counts_norm_mat <- counts_norm_mat[ sig_genes$entrez_gene_id, ]
#change the rownames of matrix to display gene names on heatmap
new_row_names <- unlist(lapply(rownames(sigGenes_counts_norm_mat), function(x) sig_genes[sig_genes$entrez_gene_id == x,]$SYMBOL))
rownames(sigGenes_counts_norm_mat) <- new_row_names

head(sigGenes_counts_norm_mat)

#heatmap col
source("gene_exp_analysis.R")
png(filename="~/apratap_bt/projects/PCBC_integrative_analysis/data/diff_exp_analysis/DESeq2_sigGenes.png",
    width=10,
    height=12,
    units="in",
    res = 600
)
draw_heatmap(log2(sigGenes_counts_norm_mat+1),scale=T,labRow=T,key=T)
dev.off()



pheatmap(t(scale(t(sigGenes_counts_norm_mat))),
         scale="none",
         annotation = annotation,
         border_color = NA
)
plotMA(diff_exp_result)
