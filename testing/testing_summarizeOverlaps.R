reads <- GAlignments(
  names = c("a","b","c","d","e","f","g"),
  seqnames = Rle(c(rep(c("chr1", "chr2"), 3), "chr1")),
  pos = as.integer(c(1400, 2700, 3400, 7100, 4000, 3100, 5200)),
  cigar = c("500M", "100M", "300M", "500M", "300M", 
            "50M200N50M", "50M150N50M"),
  strand = strand(rep("+", 7)))

gr <- GRanges(
  seqnames = c(rep("chr1", 7), rep("chr2", 4)), strand = "+", 
  ranges = IRanges(c(1000, 3000, 3600, 4000, 4000, 5000, 5400, 
                     2000, 3000, 7000, 7500), 
                   width = c(500, 500, 300, 500, 900, 500, 500, 
                             900, 500, 600, 300),
                   names=c("A", "B", "C1", "C2", "D1", "D2", "E", "F",
                           "G", "H1", "H2"))) 
groups <- factor(c(1,2,3,3,4,4,5,6,7,8,8))
grl <- splitAsList(gr, groups)
names(grl) <- LETTERS[seq_along(grl)]

## ---------------------------------------------------------------------
## Counting modes. 
## ---------------------------------------------------------------------

## First we count with a GRanges as the 'features'. Note that
## 'Union' is the most conservative counting mode followed by 
## 'IntersectionStrict' then 'IntersectionNotEmpty'.
counts1 <- 
  data.frame(union=assays(summarizeOverlaps(gr, reads))$counts, 
             intStrict=assays(summarizeOverlaps(gr, reads, 
                                                mode="IntersectionStrict"))$counts,
             intNotEmpty=assays(summarizeOverlaps(gr, reads,
                                                  mode="IntersectionNotEmpty"))$counts)

colSums(counts1)

## Split the 'features' into a GRangesList and count again.
counts2 <- 
  data.frame(union=assays(summarizeOverlaps(grl, reads))$counts, 
             intStrict=assays(summarizeOverlaps(grl, reads, 
                                                mode="IntersectionStrict"))$counts,
             intNotEmpty=assays(summarizeOverlaps(grl, reads,
                                                  mode="IntersectionNotEmpty"))$counts)
colSums(counts2)

## The GRangesList ('grl' object) has 8 features whereas the GRanges 
## ('gr' object) has 11. The affect on counting can be seen by looking
## at feature 'H' with mode 'Union'. In the GRanges this feature is 
## represented by ranges 'H1' and 'H2',
gr[c("H1", "H2")]

## and by list element 'H' in the GRangesList, 
grl["H"]

## Read "d" hits both 'H1' and 'H2'. This is considered a multi-hit when
## using a GRanges (each range is a separate feature) so the read was 
## dropped and not counted.
counts1[c("H1", "H2"), ]

## When using a GRangesList, each list element is considered a feature.
## The read hits multiple ranges within list element 'H' but only one 
## list element. This is not considered a multi-hit so the read is counted.
counts2["H", ]

## ---------------------------------------------------------------------
## Counting multi-hit reads.
## ---------------------------------------------------------------------

## The goal of the counting modes is to provide a set of rules that
## resolve reads hitting multiple features so each read is counted
## a maximum of once. However, sometimes it may be desirable to count 
## a read for each feature it overlaps. This can be accomplished by 
## setting 'inter.feature' to FALSE.

## When 'inter.feature=FALSE', modes 'Union' and 'IntersectionStrict'
## essentially reduce to countOverlaps() with type="any" and 
## type="within", respectively.

## When 'inter.feature=TRUE' only features "A", "F" and "G" have counts.
se1 <- summarizeOverlaps(gr, reads, mode="Union", inter.feature=TRUE)
assays(se1)$counts

## When 'inter.feature=FALSE' all 11 features have a count. There are 
## 7 total reads so one or more reads were counted more than once.
se2 <- summarizeOverlaps(gr, reads, mode="Union", inter.feature=FALSE)
assays(se2)$counts

## ---------------------------------------------------------------------
## Counting Bam files.
## ---------------------------------------------------------------------

library(Rsamtools)
library(pasillaBamSubset)
library("TxDb.Dmelanogaster.UCSC.dm3.ensGene")
exbygene <- exonsBy(TxDb.Dmelanogaster.UCSC.dm3.ensGene, "gene")

## (i) Single-end reads:

## 'yieldSize' can be set to iterate over large files in chunks.
## When counting, 'singleEnd' should be TRUE (default).
bf_s <- BamFile(untreated1_chr4(), yieldSize=50000)
se_s <- summarizeOverlaps(exbygene, bf_s, singleEnd=TRUE)
table(assays(se_s)$counts > 0)

## (ii) Paired-end reads:

## A paired-end file may contain singletons, reads with unmapped
## pairs or reads with more than two fragments. When 'fragments=FALSE'
## only reads paired by the algorithm are included in the counting. 
bf <- BamFile(untreated3_chr4())
se_nofrag <- summarizeOverlaps(exbygene, bf, singleEnd=FALSE, fragments=FALSE)
table(assays(se_nofrag)$counts > 0)

## When 'fragments=TRUE' (default) all singletons, reads with
## unmapped pairs and other fragments will be included in the
## counting.
bf <- BamFile(untreated3_chr4(), asMates=TRUE)
se_frag <- summarizeOverlaps(exbygene, bf, singleEnd=FALSE)
table(assays(se_frag)$counts > 0)

## As expected, using 'fragments=TRUE' results in a larger number 
## of total counts because singletons, unmapped pairs etc. are 
## included in the counting.

## Total reads in the file:
countBam(untreated3_chr4())

## Reads counted with 'fragments=FALSE':
sum(assays(se_nofrag)$counts)

## Reads counted with 'fragments=TRUE':
sum(assays(se_frag)$counts)

## ---------------------------------------------------------------------
## Count tables for DESeq or edgeR.
## ---------------------------------------------------------------------

fls <- list.files(system.file("extdata",package="GenomicRanges"),
                  recursive=TRUE, pattern="*bam$", full=TRUE)
names(fls) <- basename(fls)
bf <- BamFileList(fls, index=character(), yieldSize=1000)
genes <- GRanges(
  seqnames = c(rep("chr2L", 4), rep("chr2R", 5), rep("chr3L", 2)),
  ranges = IRanges(c(1000, 3000, 4000, 7000, 2000, 3000, 3600, 
                     4000, 7500, 5000, 5400), 
                   width=c(rep(500, 3), 600, 900, 500, 300, 900, 
                           300, 500, 500))) 
se <- summarizeOverlaps(genes, bf)

## When the reads are Bam files, the 'colData' contains summary 
## information from a call to countBam().
colData(se)

## Create count tables.
library(DESeq)
deseq <- newCountDataSet(assays(se)$counts, rownames(colData(se)))
library(edgeR)
edger <- DGEList(assays(se)$counts, group=rownames(colData(se)))

## ---------------------------------------------------------------------
## User supplied 'mode'. 
## ---------------------------------------------------------------------

## A user defined count function must have the same arguments as 
## the current counting modes.
## Not run: 
counter <- function(x, y,  ignore.strand, inter.feature) {
  ## count ...
}

se <- summarizeOverlaps(gr, reads, mode=counter) 

## End(Not run)