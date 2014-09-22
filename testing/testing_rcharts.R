require(devtools)
install_github('rCharts', 'ramnathv')
library('rCharts')

names(iris) = gsub("\\.", "", names(iris))
rPlot(SepalLength ~ SepalWidth | Species, data = iris, color = 'Species', type = 'point')


## Example 2 Facetted Barplot
hair_eye = as.data.frame(HairEyeColor)
rPlot(Freq ~ Hair | Eye, color = 'Eye', data = hair_eye, type = 'bar')


read.url <- function(url, ...){
  tmpFile <- tempfile()
  download.file(url, destfile = tmpFile, method = "curl")
  url.data <- read.csv(tmpFile, ...)
  return(url.data)
}

url <- ''
tmpFile <- tempfile()
download.file(url, destfile = tmpFile, method = "curl")
tmpF
url.data <- read.csv(tmpFile, ...)




library(rCharts)
library(plyr)
library(reshape2)
library(scales)
library("RCurl")

x <- getURL("https://raw.github.com/patilv/rChartsTutorials/master/findata.csv",.opts=list(followlocation=T))
findata <- read.csv(text=x)
head(findata)
head(ddply(melt(findata),.(variable),transform,rescale=rescale(value)))
# These are data regarding NCAA athletic department expenses at public universities. Please see the blog post where these charts were originally used 
# regarding more details on the origins of these data.: http://analyticsandvisualization.blogspot.com/2013/10/subsidies-revenues-and-expenses-of-ncaa.html
findata=findata[,-c(1:2)] # removing first dummy column - the csv quirk - and second column on Rank.
findatamelt=ddply(melt(findata),.(variable),transform,rescale=rescale(value))
hmap <- rPlot(variable ~ School, color = 'rescale', data = findatamelt, type = 'tile')
hmap$addParams(height = 400, width=1000)
hmap$guides(reduceXTicks = FALSE)
hmap$guides("{color: {scale: {type: gradient, lower: white, upper: red}}}")
hmap$guides(y = list(numticks = length(unique(findatamelt$value))))
hmap$guides(x = list(numticks = 5))
hmap

?read.csv

sessionInfo()
