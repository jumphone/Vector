
setwd('F:/Vector/data/HumanHepDiff_GSE81252')

D=read.csv('GSE81252_data.cast.log2.lineage.csv',header=T,row.names=1)

DATA=t(D[,2:ncol(D)])
TYPE=D[,1]


source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

############################################
# Analyze all cells
pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)
pbmc@meta.data$type=TYPE
pbmc@meta.data$type[which(pbmc@meta.data$type %in% c('he1','he2'))]='HE'
pbmc@meta.data$type[which(pbmc@meta.data$type %in% c('de'))]='DE'
pbmc@meta.data$type[which(pbmc@meta.data$type %in% c('ih1'))]='IH'
pbmc@meta.data$type[which(pbmc@meta.data$type %in% c('ipsc'))]='IPSC'
pbmc@meta.data$type[which(pbmc@meta.data$type %in% c('mh1','mh2'))]='MH'
