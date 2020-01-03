

##########################################
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseDentateGyrus_GSE104323/')
pbmc=readRDS(file='pbmc.RDS')

VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings
PCA=vector.rankPCA(PCA)


OUT=vector.buildGrid(VEC, N=40,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)


##############
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
##############


OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL,AL=40)

FeaturePlot(pbmc,features='Nes')



############################
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/LUOZAILI')
pbmc=readRDS('12W_pbmc3k_final.rds')

FeaturePlot(pbmc,features=c('nFeature_RNA'))





VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings

PCA=vector.rankPCA(PCA)


OUT=vector.buildGrid(VEC, N=20,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)

OUT=vector.getValue(OUT, PCA, SHOW=TRUE)

#OUT$VALUE[OUT$USED_INDEX]=vector.removeOut(OUT$VALUE[OUT$USED_INDEX])
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
##############

OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL,AL=40)







