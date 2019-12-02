
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/AIM/MouseDentateGyrus_GSE104323/')


pbmc=readRDS(file='pbmc.RDS')


#umap=BEER.bbknn(pbmc, PCUSE, NB=2, NT=10)
#pbmc@reductions$umap@cell.embeddings=umap
#DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1,label=F)

DimPlot(pbmc, group.by='type',label = TRUE)+NoLegend()
#FeaturePlot(pbmc,features=c('Nes','Olig2','Pdgfra','Gfap','Gdnf'))
#######################################


VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings
#######################

OUT=vector.buildGrid(VEC, N=40,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
VALUE=vector.getValue(PCA)
OUT=vector.gridValue(OUT,VALUE, SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
#OUT=vector.selectCenter(OUT)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)

#######################
OUT=vector.reDrawArrow(OUT, COL=OUT$COL)
OUT=vector.selectRegion(OUT)


