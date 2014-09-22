hgid_to_keggID <- read.table("~/apratap_bt/reference/id_converts/hsa_hgnc.list",sep="\t")

#delete the third col
hgid_to_keggID$V3 <- NULL
colnames(hgid_to_keggID) <- c("kegg_gene_id","hugo_gene_id")

hgid_to_keggID$kegg_gene_id <- sub('hsa:', '',hgid_to_keggID$kegg_gene_id)
hgid_to_keggID$hugo_gene_id <- sub('hgnc:', '',hgid_to_keggID$hugo_gene_id)

x <- sample(hgid_to_keggID$kegg_gene_id,100,replace=T)
paste(x,collapse=',')


http://pathways.embl.de/mapExport.cgi?batch=1&amp;jobID=G0rHoPd9VzeaIo9ABLdrpA&amp;include_metabolic=true&amp;include_regulatory=true&amp;&amp;export_dpi=72&amp;export_format=svg



http://pathways.embl.de/mapExport.cgi?batch=1&amp;jobID=2zkmRUiHOuGpzMvH4BORcQ&amp;include_metabolic=true&amp;include_regulatory=true&amp;&amp;export_dpi=72&amp;export_format=png