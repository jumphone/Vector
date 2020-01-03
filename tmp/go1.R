

##########################################
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseDentateGyrus_GSE104323/')
pbmc=readRDS(file='pbmc.RDS')

VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings

OUT=vector.buildGrid(VEC, N=40,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)


##############
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)


plot(VEC,col='grey90')
i=2
while(i<=length(OUT$VALUE)){
    var
    #if(i %% 100==1){
    #    points(VEC[order(-OUT$VALUE)[1:i],])
    #    Sys.sleep(0.2)
    #}
    i=i+1}
      
      





YES=cbind(VEC,OUT$VALUE)[OUT$USED_INDEX,]
S.YES=apply(YES,2,scale)

library(gmodels)
PCA.OUT=fast.prcomp(S.YES, retx = TRUE, center = FALSE, scale. = FALSE, tol = NULL)
YES.PCA=PCA.OUT$x



source('https://raw.githubusercontent.com/jumphone/VISA/master/VISA.R')

visa.plot3d(S.YES,COL='black')



OUT$VALUE[OUT$USED_INDEX]=YES.PCA[,1]#vector.removeOut(OUT$VALUE[OUT$USED_INDEX])
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
##############


OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL,AL=40)





############################
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/LUOZAILI')
pbmc=readRDS('12W_pbmc3k_final.rds')

FeaturePlot(pbmc,features=c('nFeature_RNA'))





VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings


OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)

OUT=vector.getValue(OUT, PCA, SHOW=TRUE)

#OUT$VALUE[OUT$USED_INDEX]=vector.removeOut(OUT$VALUE[OUT$USED_INDEX])
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
##############







