source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

setwd('F:/Vector/data/kBET/')

#/home/zy/single_cell/mouse_EED
DATA=read.csv('Counts.csv',header=TRUE,row.names=1,sep=',')

DATA=t(DATA)


META=read.csv('pData.csv',header=TRUE,row.names=1,sep=',')
TYPE=META$sample
BATCH=META$batch

saveRDS(DATA,file='DATA.RDS')
saveRDS(META,file='META.RDS')

####################################################################
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/kBET/')

DATA=readRDS(file='DATA.RDS')
META=readRDS(file='META.RDS')

TYPE=META$sample
BATCH=META$batch

        library(sva)
        library(limma)
        pheno = data.frame(batch=as.matrix(BATCH))
        orig.data=DATA
        used.gene.index=1:nrow(DATA)#which(rownames(orig.data) %in% VARG)
        edata = as.matrix(orig.data)[used.gene.index,]
        batch = pheno$batch
        modcombat = model.matrix(~1, data=pheno)
        combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
        rownames(combat_edata)=rownames(edata)
        colnames(combat_edata)=colnames(edata)
        combat_edata=as.matrix(combat_edata)
        combat_edata[which(combat_edata<0)]=0



pbmc <- CreateSeuratObject(counts = combat_edata, project = "pbmc3k", min.cells = 0, min.features = 0)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs = 150)




pbmc@meta.data$type=TYPE


pbmc <- RunUMAP(pbmc, dims = 1:150)
DimPlot(pbmc, reduction = "umap",group.by='type',label=TRUE)



saveRDS(pbmc,file='pbmc.RDS')


PCA=pbmc@reductions$pca@cell.embeddings[,c(1:150)]
VEC=pbmc@reductions$umap@cell.embeddings#pbmc@reductions$pca@cell.embeddings[,c(1, 2)]



OUT=vector.buildGrid(VEC, N=15,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)














