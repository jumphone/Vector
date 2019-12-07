source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

setwd('F:/Vector/data/HumanHepDiff_GSE81252')

D=read.csv('GSE81252_data.cast.log2.lineage.csv',header=T,row.names=1)

DATA=t(D[,2:ncol(D)])
TYPE=as.character(D[,1])

TYPE=c(as.character(TYPE))
TYPE[which(TYPE %in% c('ipsc'))]='IPS'
TYPE[which(TYPE %in% c('de'))]='DE'
TYPE[which(TYPE %in% c('he1','he2'))]='HE'
TYPE[which(TYPE %in% c('ih1'))]='IH'
TYPE[which(TYPE %in% c('mh1','mh2'))]='MH'


pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
#pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs = 150)

#pbmc <- RunUMAP(pbmc, dims = 1:100, n.neighbors=10, min.dist=0.5,spread=1)

VARG=VariableFeatures(object = pbmc)
library(sva)
library(limma)
pheno = data.frame(batch=as.matrix(TYPE))
orig.data=pbmc@assays$RNA@data
used.gene.index=which(rownames(orig.data) %in% VARG)
edata = as.matrix(orig.data)[used.gene.index,]
batch = pheno$batch
modcombat = model.matrix(~1, data=pheno)
combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

EXP=combat_edata
EXP[which(EXP<0)]=0

pbmc <- CreateSeuratObject(counts = EXP, project = "pbmc3k", min.cells = 0, min.features = 0)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs = 150)
#pbmc <- RunUMAP(pbmc, dims = 1:20, n.neighbors=10, min.dist=0.5,spread=1)
pbmc <- RunUMAP(pbmc, dims = 1:70)

pbmc@meta.data$type=TYPE
DimPlot(pbmc, reduction = "umap",group.by='type',label=TRUE)

VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings[,1:5]


OUT=vector.buildGrid(VEC, N=10,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)

OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)

OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)






























table(TYPE)
############################################

USED.1=c(which(TYPE=='IPS')[1:15], which(TYPE=='DE')[1:30], which(TYPE=='HE')[1:45], 
         which(TYPE=='IH')[1:60], which(TYPE=='MH')[1:75])

D.1=DATA[,USED.1]
T.1=TYPE[USED.1]

# Analyze all cells
pbmc <- CreateSeuratObject(counts = D.1, project = "pbmc3k", min.cells = 0, min.features = 0)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs = 150)

pbmc <- RunUMAP(pbmc, dims = 1:100, n.neighbors=10, min.dist=0.5,spread=1)

pbmc@meta.data$type=T.1
DimPlot(pbmc, reduction = "umap",group.by='type',label=TRUE)



EXP=as.matrix(pbmc@assays$RNA@scale.data)
EXP=apply(EXP,2,rank)
library(destiny, quietly = TRUE)
dm <- DiffusionMap(t(EXP))
plot(dm$DC1,dm$DC2)

UMAP=pbmc@reductions$umap@cell.embeddings
DM=cbind(dm$DC1, dm$DC2)
DM=apply(DM,2,scale)
colnames(DM)=colnames(UMAP)
rownames(DM)=rownames(UMAP)


pbmc@reductions$umap@cell.embeddings=DM

DimPlot(pbmc, reduction = "umap",group.by='type',label=TRUE)



VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings[,1:5]


OUT=vector.buildGrid(VEC, N=10,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)

OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)

OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)




































PCA= pbmc@reductions$pca@cell.embeddings

EXP=as.matrix(pbmc@assays$RNA@scale.data)
EXP=t(apply(t(EXP),2,rank))
#DimPlot(pbmc, reduction = "pca",dims=c(1,2),group.by='type',label=TRUE)


library(destiny, quietly = TRUE)
#dm <- DiffusionMap(PCA[,1:150])
dm <- DiffusionMap(t(EXP))
plot(dm$DC1,dm$DC2)

UMAP=pbmc@reductions$umap@cell.embeddings
DM=cbind(dm$DC1, dm$DC2)
DM=apply(DM,2,scale)
colnames(DM)=colnames(UMAP)
rownames(DM)=rownames(UMAP)

pbmc@reductions$umap@cell.embeddings=DM



VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings


OUT=vector.buildGrid(VEC, N=20,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)


DimPlot(pbmc, reduction = "umap",group.by='type',label=TRUE)










PCA= pbmc@reductions$pca@cell.embeddings


EXP=as.matrix(pbmc@assays$RNA@scale.data)
EXP=t(apply(t(EXP),2,rank))


DimPlot(pbmc, reduction = "pca",dims=c(1,2),group.by='type',label=TRUE)




library(destiny, quietly = TRUE)
#dm <- DiffusionMap(PCA[,1:150])
dm <- DiffusionMap(t(EXP))
plot(dm$DC1,dm$DC2)

UMAP=pbmc@reductions$umap@cell.embeddings
DM=cbind(dm$DC1, dm$DC2)
DM=apply(DM,2,scale)
colnames(DM)=colnames(UMAP)
rownames(DM)=rownames(UMAP)


pbmc@reductions$umap@cell.embeddings=DM
#PCA= pbmc@reductions$pca@cell.embeddings
#R.PCA=apply(PCA,2,rank)
#pbmc@reductions$pca@cell.embeddings=R.PCA

DimPlot(pbmc, reduction = "umap",group.by='type',label=TRUE)

saveRDS(pbmc,file='pbmc.RDS')


#########################################################




VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings


OUT=vector.buildGrid(VEC, N=20,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)


####################################################################











source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseOligo_GSE75330')

DATA=as.matrix(readRDS('DATA.RDS'))
TYPE=readRDS('GSE75330.LABEL.RDS')

DATA=.simple_combine(DATA,RD)$combine

RD=matrix( sample(as.numeric(DATA), 100*nrow(DATA), replace = TRUE), ncol=100, nrow=nrow(DATA))
rownames(RD)=rownames(DATA)
colnames(RD)=paste0('RD_',1:100)

pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)
pbmc@meta.data$type=c(as.character(TYPE),rep('RD',100))



pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 3000)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs = 150)

pbmc <- RunUMAP(pbmc, dims = 1:100, n.neighbors=10, min.dist=0.5,spread=1)

DimPlot(pbmc, reduction = "umap",group.by='type',label=TRUE)





VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings


OUT=vector.buildGrid(VEC, N=20,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)












