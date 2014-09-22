counts_file <- star_read1_counts_file

counts = read_counts(counts_file)
counts$TRA00010814 <- NULL
studyDesign <- data.frame(row.names = colnames(counts),
                          condition = c("ablated","ablated","control","ablated",
                                        "control","control"))
#changing the level so that log2 ratio makes more sense and are w.r.t control
studyDesign$condition <- factor(studyDesign$condition, level=c("control", "ablated"))

dds <- DESeqDataSetFromMatrix(countData=counts,
                              colData = studyDesign,
                              design = ~ condition)
cat('Running DESeq2\n')
#DESeq2 based differential expression
dds <- DESeq(dds)
res <- as.data.frame(results(dds))

#filter out genes with padj == NA
res <- res[complete.cases(res),]
diffExp_genes <- res[res$padj < .10, ]-
write.table(res,file=paste0(counts_file,'.difExp_results_TRA00010814_removed.txt'),quote=F, sep="\t")
write.table(diffExp_genes,file=paste0(counts_file,'.difExp_genes_TRA00010814_removed_pvalue_lt.10.txt'),quote=F, sep="\t")

#Volcano plot
cat('Volcano plot \n')
res['significance'] = 'not-sig'
res$significance[res$padj < min_padj] =  paste('padj <',min_padj)
ggplot(data=res,aes(x=log2FoldChange, y=-log(padj), color=significance)) + geom_jitter(size=.8)
ggsave(paste0(counts_file,'.volcano_plot_TRA00010814_removed.png'),width=8,height=8,units="in")

#MA plot
cat('plotting the MA plot\n')
png(paste0(counts_file,'.DE_MA_plot_TRA00010814_removed.png'))
plotMA(dds)
dev.off()

#regularized log transformation
cat('creating the regularized log2 counts\n')
rld <- rlogTransformation(dds, blind=TRUE)
write.table(assay(rld),file=paste0(counts_file,'.regularizedLog_counts_TRA00010814_removed.txt'), quote=F, sep="\t")
#variance stabilized transformation
cat('creating the variance stabilized counts\n')
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
write.table(assay(rld),file=paste0(counts_file,'.varstabilized_counts_TRA00010814_removed.txt'), quote=F, sep="\t")

#PCA plot
cat('PCA plot \n')
#save the lattice plot
png(paste0(counts_file,'.PCA_plot_TRA00010814_removed.png'),width=8,height=6,units="in",res=150)
names <- rownames(colData(vsd))
pcaplot <- plotPCA(vsd,ntop=1000)
pcaplot = update(pcaplot, panel = function(x, y, ...) {
  lattice::panel.xyplot(x, y, ...);
  lattice::ltext(x=x, y=y, labels=names, pos=1, offset=1, cex=0.6)
})
print(pcaplot)
dev.off()





###############
#Non paramteric analysis with SAMR (Significance Analysis of Microarrays)
###############



samples = c('TRA00010811' = 2,
            'TRA00010812' = 2,
            'TRA00010813' = 2,
            'TRA00010814' = 2,
            'TRA00010815' = 1,
            'TRA00010816' = 1,
            'TRA00010817' = 1,
            'TRA00010818' = 1)


sample_names = sapply(colnames(counts), function(x){samples[[x]] })

samr_out <- SAMseq(counts, sample_names, 
                   geneid=rownames(star_read1_counts),
                   genenames=rownames(star_read1_counts),
                   resp.type='Two class unpaired',
                   fdr.output=.10)
#names(samr_out)
#names(samr_out[['siggenes.table']])
#head(samr_out[['samr.obj']])

genes_up <- samr_out[['siggenes.table']][['genes.up']]
genes_down <- samr_out[['siggenes.table']][['genes.lo']]

write.table(genes_up,file=paste0(star_read1_counts_file,'nonparametric_diffExp_genes_padj_lt.10.tsv'), quote=F, sep="\t", row.names=F)

head(genes_up,n=100)
cat(genes_up[,1])