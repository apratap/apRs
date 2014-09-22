library(plyr)
library(reshape)
library(ggplot2)
library("Rsamtools")
library("rtracklayer")

#setup the parallel environment
options(mc.cores=8)

base_analysis_dir <- "~/projects/PCBC_integrative_analysis/"

#get a list of all the bam files
all_bamFiles <- dir(base_analysis_dir,pattern="accepted_hits.bam$",recursive=T, full.names=T)

##########
#COUNTING
##########
#get the annotation
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
mm10_txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#showMethods(class=class(mm10_txdb))
mm10_refGenes_by_exons <- exonsBy(mm10_txdb,by="gene")


###############
#Using read 1
###############
idx <- grep('tophat_mapped_read1',all_bamFiles)
to_use_bamfiles <- all_bamFiles[idx]
bamFilesList_obj <- BamFileList(to_use_bamfiles,yieldSize=2000000,index=character())

#restrict the counting to primary alignment of paired reads
#selecting the params in the bam needed for counting
bam_params <- c("rname","pos","cigar","strand")
params <- ScanBamParam(what=bam_params)


#do the counting
read1_counts <- summarizeOverlaps(features=mm10_refGenes_by_exons,
                                  reads = bamFilesList_obj,
                                  mode = "Union",
                                  ignore.strand=TRUE,
                                  singleEnd = TRUE,
                                  params = params)
read1_counts <- assays(read1_counts)$counts
new_colnames <- gsub("^.*(TRA.+?)/.+$", '\\1',colnames(read1_counts),perl=T)
colnames(read1_counts) <- new_colnames
head(read1_counts)

#save read 1 counts
outfile = "~/projects/PCBC_integrative_analysis/analysis/read1/read1_based_counts.tsv"
write.table(read1_counts, file=outfile, sep="\t",quote=FALSE)


###############
#Using read 2
###############
idx <- grep('tophat_mapped_read2',all_bamFiles)
read2_bamfiles <- all_bamFiles[idx]
read_2_bamFilesList_obj <- BamFileList(read2_bamfiles,yieldSize=2000000,index=character())

#restrict the counting to primary alignment of paired reads
#selecting the params in the bam needed for counting
bam_params <- c("rname","pos","cigar","strand")
params <- ScanBamParam(what=bam_params)

#do the counting
read2_counts <- summarizeOverlaps(features=mm10_refGenes_by_exons,
                                  reads = read_2_bamFilesList_obj,
                                  mode = "Union",
                                  ignore.strand=TRUE,
                                  singleEnd = TRUE,
                                  params = params)

read2_counts <- assays(read2_counts)$counts
new_colnames <- gsub("^.*(TRA.+?)/.+$", '\\1',colnames(read2_counts),perl=T)
colnames(read2_counts) <- new_colnames

#save read 2 counts
outfile = "~/projects/PCBC_integrative_analysis/analysis/read2/read2_based_counts.tsv"
write.table(read2_counts, file=outfile, sep="\t",quote=FALSE)


###############
#Using read 1 & 2 (Paired reads)
###############
idx <- grep('tophat_mapped_paired_reads',all_bamFiles)
paired_bamfiles <- all_bamFiles[idx]
paired_bamFilesList_obj <- BamFileList(paired_bamfiles,yieldSize=1000000,index=character())

#restrict the counting to primary alignment of paired reads
#selecting the params in the bam needed for counting
bam_params <- c("rname","pos","cigar","strand")
params <- ScanBamParam(what=bam_params)

#do the counting
#setup the parallel environment
options(mc.cores=4)
pairedread_counts <- summarizeOverlaps(features=mm10_refGenes_by_exons,
                                      reads = paired_bamFilesList_obj,
                                      mode = "Union",
                                      ignore.strand=TRUE,
                                      singleEnd = FALSE,
                                      fragments = TRUE,
                                      params = params)

pairedread_counts <- assays(pairedread_counts)$counts
new_colnames <- gsub("^.*(TRA.+?)/.+$", '\\1',colnames(pairedread_counts),perl=T)
colnames(pairedread_counts) <- new_colnames

#save read 1 & 2 counts
outfile = "~/projects/PCBC_integrative_analysis/analysis/paired/pairedread_based_counts.tsv"
write.table(pairedread_counts, file=outfile, sep="\t",quote=FALSE)

