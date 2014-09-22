library("devtools")
library("GenomicRanges")
library("synapseClient")
library("Rsamtools")


######
#test case
######

#annotation
gene1 <- GRanges("chr1", IRanges(c(1,40),c(20,60)),"+")
gene2 <- GRanges("chr2", IRanges(c(100,140),c(120,160)),"-")
annotation <- GRangesList("gene1"=gene1)
annotation


#alignments
r1 <- GAlignments("chr1",1L,"10M",strand("-"))
r2 <- GAlignments("chr1",65L,"10M",strand("+"))
pair <- GAlignmentPairs(r1,r2,isProperPair=TRUE)
pair

result <- summarizeOverlaps(annotation,pair,ignore.strand=FALSE)
assays(result)$count



?s

class(result)
?SummarizedExperiment
result@rowData
rowData(result)
?summarizeOverlaps
showMethods(class=class(result) )

hub <- AnnotationHub()
cols(hub)
keys(hub,keytype="Species")
keys(hub,keytype="Genome")

#adding filters
filters(hub) = list(Species="Homo sapiens")

refGenes <- hub$ensembl.release.70.gtf.homo_sapiens.Homo_sapiens.GRCh37.70.gtf_0.0.1.RData

refGenes

txdb

columns(txdb)
keytypes(txdb)
x <- genes(txdb)
x[1:3]

class(txdb)
genes(txdb)
showMethods(class=class(txdb),where=search())
metadata(txdb)
seqinfo(txdb)
seqnameStyle(txdb)
columns(txdb)
transcripts(txdb)
exons(txdb)
cds(txdb)


names(refGenes_by_exons)
seqnames(refGenes_by_exons[1:5])



test <- refGenes_by_exons[1:5]
length(test)
test[[1]]
elementLengths(test)
sum(width(reduce(test)))
length(refGenes_by_exons[1:100])

sapply(refGenes_by_exons[1:10],function(x) sum(width(x)))
sum(width(refGenes_by_exons[1]))
elementLengths(refGenes_by_exons)
showMethods(class=class(refGenes_by_exons))



#########
#get the annotation
#########
library("GenomicFeatures")
source("http://bioconductor.org/biocLite.R")
biocLite("TxDb.Mmusculus.UCSC.mm10.knownGene")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
mm10_txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
refGenes_by_exons <- exonsBy(mm10_txdb,by="gene")

#get a bam file list for counting
data_base_dir <- "~/apratap_bt/projects/PCBC_integrative_analysis/data/"
bam_files_list <- dir(data_base_dir,pattern="accepted_hits.bam",recursive=T)
bam_files_list <- BamFileList(bam_files_list,yieldSize=2000000)




?BamFileList

#restrict the counting to primary alignment of paired reads
flag <- scanBamFlag(isNotPrimaryRead=FALSE,isProperPair=TRUE)
#selecting the params in the bam needed for counting
bam_params <- c("rname","pos","cigar")
params <- ScanBamParam(what=bam_params,flag=flag)

sum

?summarizeOverlaps
summarizeOverlaps(refGenes_by_exons)


?system.file
synapseLogin()
syn_sampleBam <- synGet('syn2246875')
sampleBam <- syn_sampleBam@filePath 

?summarizeOverlaps
?findOverlaps

summarizeOverlaps(features = ,
                  reads = BamFileList(sampleBam),
                  ignore.strand = TRUE,
                  
                  )

?summarizeOverlaps,GRanges,BamFileList-method

?mcols


library("GenomicRanges")
library("Rsamtools")
GAlignments("chr1",1L,"10M",strand("+"))

sessionInfo()
