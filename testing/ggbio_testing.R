library(ggbio)
p.ideo <- Ideogram(genome = "hg19")
p.ideo


library(GenomicRanges)
## special highlights instead of zoomin!
p.ideo + xlim(GRanges("chr2", IRanges(1e8, 1e8+10000000)))

library("Homo.sapiens")
class(Homo.sapiens)
## [1] "OrganismDb"
## attr(,"package")
## [1] "OrganismDbi"
##
data(genesymbol, package = "biovizBase")
wh <- genesymbol[c("BRCA1", "NBR1")]
wh <- range(wh, ignore.strand = TRUE)
p.txdb <- autoplot(Homo.sapiens, which = wh)
p.txdb
