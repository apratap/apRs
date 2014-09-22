#!/usr/bin/env Rscript
library(plyr)
library(reshape)
library(ggplot2)
library("Rsamtools")

#setup the parallel environment
options(mc.cores=6)

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
#Using read 1 & 2 (Paired reads)
###############
idx <- grep('tophat_mapped_paired_reads',all_bamFiles)
paired_bamfiles <- all_bamFiles[idx]
#paired_bamfiles
#paired_bamfiles <- paired_bamfiles[6]
paired_bamFilesList_obj <- BamFileList(paired_bamfiles,yieldSize=1000000,index=character())

#restrict the counting to primary alignment of paired reads
#selecting the params in the bam needed for counting
bam_params <- c("rname","pos","cigar","strand")
params <- ScanBamParam(what=bam_params)

#do the counting
pairedread_counts <- summarizeOverlaps(features=mm10_refGenes_by_exons,
                                       reads = paired_bamFilesList_obj,
                                       mode = "Union",
                                       ignore.strand=TRUE,
                                       singleEnd = FALSE,
                                       fragments = TRUE,
                                       params = params)


pairedread_counts <- assays(pairedread_counts)$counts
#sum(pairedread_counts)
new_colnames <- gsub("^.*(TRA.+?)/.+$", '\\1',colnames(pairedread_counts),perl=T)
colnames(pairedread_counts) <- new_colnames

#rownames(pairedread_counts) == '100039795'

#pairedread_counts['23793']


#save read 1 & 2 counts
outfile = "~/projects/PCBC_integrative_analysis/analysis/paired/pairedread_Rbased_counts.tsv"
write.table(pairedread_counts, file=outfile, sep="\t",quote=FALSE)

