install.packages("devtools") # you don't need to run this command if you already have the devtools package installed.

library(devtools)
install_github("clickme", "nachocab")

library(clickme)
clickme(points, rnorm(100))

data(microarray)
clickme(points, x = microarray$significance, y = microarray$logFC,
        color_groups = ifelse(microarray$adj.P.Val < 1e-4, "Significant", "Noise"),
        names = microarray$gene_name,
        xlab = "Significance (-log10)", ylab = "Fold-change (log2)",
        extra = list(Probe = microarray$probe_name))