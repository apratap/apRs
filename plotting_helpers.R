library(gridExtra)
library(ggplot2)

tableToPlot <- function(xyTable){
  qplot(0:1, 0:1, geom = "blank") + 
    theme_bw() +
    theme(line = element_blank(),
          text = element_blank()) +
    annotation_custom(grob = tableGrob(xyTable,
                                       # change font sizes:
                                       gpar.coltext = gpar(cex = 1),
                                       gpar.rowtext = gpar(cex = 1)),
                      xmin = 0, xmax = 1, ymin = 0, ymax = 1) 
}



dvipng.dvi <-
  function (object, file, ...) 
  {
    cmd <- if (missing(file)) 
      paste("dvipng -T tight", shQuote(object$file))
    else paste("dvipng -T tight", "-o", file, shQuote(object$file))
    invisible(system(cmd))
  }



panel.correlation <- function(x, y, digits=2, prefix="", cex.cor, col="black",...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y,method="spearman")
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  if(abs(r) > .8 ) cex.cor = cex.cor + .5
  if(abs(r) > .9) col = "red"
  text(0.5, 0.5, txt, cex = cex.cor, col = col )
}

panel.smooth <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
                          cex = 1, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
          col = col.smooth, ...)
}



corPlot <- function(data_matrix,...){
  m <- as.matrix(data_matrix)
  #remove any rows with NA
  to_keep <- !apply(m,1, function(x) any(is.na(x)))
  m <- m[to_keep,]
  pairs(m, upper.panel=panel.correlation, lower.panel=panel.smooth, ...)
}



# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}