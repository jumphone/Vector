
########################################
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseIntestine_GSE92332/')




DATA=readRDS('REF.RDS')
saveRDS(DATA,file='DATA.RDS')

########################################
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseIntestine_GSE92332/')








DATA=readRDS('DATA.RDS')
TYPE=readRDS('LABEL.RDS')
pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)
BATCH=pbmc@meta.data$orig.ident



used_cell_index=which((BATCH %in% c('B3')) & TYPE %in% c('Enterocyte.Immature.Distal',
                                 'Enterocyte.Immature.Proximal',
                                 'Enterocyte.Mature.Distal',
                                  'Enterocyte.Mature.Proximal',
                                  'Enterocyte.Progenitor',
                                  'Enterocyte.Progenitor.Early',
                                  'Enterocyte.Progenitor.Late',
                                  'Stem',
                                  'TA.Early',
                                  'TA.G1','TA.G1'
                                 ))

###################################
DATA=DATA[,used_cell_index]
TYPE=TYPE[used_cell_index]
BATCH=BATCH[used_cell_index]
####################################################

pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs=50)
pbmc <- RunUMAP(pbmc, dims = 1:50)
DimPlot(pbmc, reduction = "umap")
######################
pbmc@meta.data$type=TYPE
DimPlot(pbmc, reduction = "umap",group.by='type',label=TRUE)


VEC=pbmc@reductions$umap@cell.embeddings
#VEC=pbmc@reductions$pca@cell.embeddings[,1:2]
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings#[,10:50]

#FeaturePlot(pbmc,features = c('Lgr5','Itgb1','Smoc2','Sox9','Bmi1','Lrig1','Ascl2'),ncol=3)

FeaturePlot(pbmc,features = c('Lgr5'))
OUT=VECTOR(VEC, PCA, N=20)

LGR5=pbmc@assays$RNA@data[which(rownames(pbmc)=='Lgr5'),]


#########################




mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=100, ROUND=1, GN=5000, SEED=1, COMBAT=TRUE, RMG=NULL) 


#mybeer=readRDS(file='mybeer.RDS')
PCUSE=mybeer$select
COL=rep('black',length(mybeer$cor))
COL[PCUSE]='red'
plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,
    xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))
#####################
pbmc <- mybeer$seurat
PCUSE=mybeer$select   
pbmc=BEER.combat(pbmc) #Adjust PCs using ComBat
umap=BEER.bbknn(pbmc, PCUSE, NB=3, NT=10)

pbmc@reductions$umap@cell.embeddings=umap
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1,label=F)
pbmc@meta.data$type=TYPE
DimPlot(pbmc, reduction = "umap",group.by='type')
DimPlot(pbmc, reduction = "pca",group.by='type')




VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings


SCORE=apply(apply(abs(PCA),2,rank),1,mean)


#OUT=VECTOR(VEC, PCA, N=20)

plot(SCORE,LGR5)



SCORE=apply((PCA.N-PCA.C.MAT)^2,1,mean)
SCORE=(SCORE-min(SCORE))/ (max(SCORE)-min(SCORE))
SCORE=1-SCORE

plot(SCORE,LGR5)

OUT=vector.buildGrid(VEC, N=40,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
VALUE=vector.getScore(PCA)
OUT=vector.gridValue(OUT,VALUE, SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
#OUT=vector.selectCenter(OUT)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)

#######################
OUT=vector.reDrawArrow(OUT, COL=OUT$COL)
OUT=vector.selectRegion(OUT)


















ppp=DimPlot(pbmc, reduction.use='umap', pt.size=0.5)
used.cells <- CellSelector(plot = ppp)

##############
used_cell_index=which(colnames(pbmc) %in% used.cells)
DATA=DATA[,used_cell_index]
TYPE=TYPE[used_cell_index]
BATCH=pbmc@meta.data$orig.ident[used_cell_index]


pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs=50)
pbmc <- RunUMAP(pbmc, dims = 1:50)
DimPlot(pbmc, reduction = "umap")
pbmc@meta.data$type=TYPE
DimPlot(pbmc, reduction = "umap",group.by='type')





VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings




VECTOR(VEC, PCA, N=40)















used_cell_index=which(TYPE %in% c('Enterocyte.Immature.Distal',
                                 'Enterocyte.Immature.Proximal',
                                 'Enterocyte.Mature.Distal',
                                  'Enterocyte.Mature.Proximal',
                                  'Enterocyte.Progenitor',
                                  'Enterocyte.Progenitor.Early',
                                  'Enterocyte.Progenitor.Late',
                                  'Stem',
                                  'TA.Early',
                                  'TA.G1','TA.G1'
                                 ))

###################################
DATA=DATA[,used_cell_index]
TYPE=TYPE[used_cell_index]
BATCH=pbmc@meta.data$orig.ident[used_cell_index]
####################################################



mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=20, ROUND=1, GN=2000, SEED=1, COMBAT=TRUE, RMG=NULL) 

#mybeer=readRDS(file='mybeer.RDS')
PCUSE=mybeer$select
COL=rep('black',length(mybeer$cor))
COL[PCUSE]='red'
plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,
    xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))
#####################
pbmc <- mybeer$seurat
PCUSE=mybeer$select   
pbmc=BEER.combat(pbmc) #Adjust PCs using ComBat
umap=BEER.bbknn(pbmc, PCUSE, NB=3, NT=10)

pbmc@reductions$umap@cell.embeddings=umap
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1,label=F)
pbmc@meta.data$type=TYPE
DimPlot(pbmc, reduction = "umap",group.by='type')




ppp=DimPlot(pbmc, reduction.use='umap', pt.size=0.5)
used.cells <- CellSelector(plot = ppp)















pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs=20)
pbmc <- RunUMAP(pbmc, dims = 1:20)
DimPlot(pbmc, reduction = "umap")

pbmc@meta.data$type=TYPE
DimPlot(pbmc, reduction = "umap",group.by='type')




VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings[,mybeer$select]



VECTOR(VEC, PCA, N=15)


















ppp=DimPlot(pbmc, reduction.use='umap', pt.size=0.5)
used.cells <- CellSelector(plot = ppp)


used_cell_index=which(colnames(pbmc) %in% used.cells & TYPE %in% c('Enterocyte.Immature.Distal',
                                 'Enterocyte.Immature.Proximal',
                                 'Enterocyte.Mature.Distal',
                                  'Enterocyte.Mature.Proximal',
                                  'Enterocyte.Progenitor',
                                  'Enterocyte.Progenitor.Early',
                                  'Enterocyte.Progenitor.Late',
                                  'Stem',
                                  'TA.Early',
                                  'TA.G1','TA.G1'
                                 ))

###################################
DATA=DATA[,used_cell_index]
TYPE=TYPE[used_cell_index]
BATCH=pbmc@meta.data$orig.ident[used_cell_index]
####################################################


pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs=50)
pbmc <- RunUMAP(pbmc, dims = 1:50)
DimPlot(pbmc, reduction = "umap")

pbmc@meta.data$type=TYPE
DimPlot(pbmc, reduction = "umap",group.by='type')










mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=50, ROUND=1, GN=5000, SEED=1, COMBAT=TRUE, RMG=NULL)  
saveRDS(mybeer, file='mybeer.RDS')




#mybeer=readRDS(file='mybeer.RDS')
PCUSE=mybeer$select
COL=rep('black',length(mybeer$cor))
COL[PCUSE]='red'
plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,
    xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))
#####################
pbmc <- mybeer$seurat
PCUSE=mybeer$select   
pbmc=BEER.combat(pbmc) #Adjust PCs using ComBat
umap=BEER.bbknn(pbmc, PCUSE, NB=3, NT=10)
pbmc@reductions$umap@cell.embeddings=umap
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1,label=F)
pbmc@meta.data$type=TYPE
DimPlot(pbmc, reduction = "umap",group.by='type')
#####################
saveRDS(pbmc,file='pbmc.RDS')
#################################################



##########################################
pbmc@meta.data$type[which(pbmc@meta.data$type %in% c('Enterocyte.Immature.Distal',
                                                     'Enterocyte.Immature.Proximal'))]='Enterocyte.Immature'
pbmc@meta.data$type[which(pbmc@meta.data$type %in% c('Enterocyte.Mature.Distal',
                                                     'Enterocyte.Mature.Proximal'))]='Enterocyte.Mature'

pbmc@meta.data$type[which(pbmc@meta.data$type %in% c('Enterocyte.Progenitor',
                                  'Enterocyte.Progenitor.Early',
                                  'Enterocyte.Progenitor.Late'))]='Enterocyte.Progenitor'


pbmc@meta.data$type[which(pbmc@meta.data$type %in% c('TA.Early',
                                  'TA.G1','TA.G1'))]='TA'
##########################################

DimPlot(pbmc, reduction = "umap",group.by='type')


###################################

VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings

SCORE=apply(apply(abs(PCA),2,rank),1,mean)

pbmc@meta.data$score=apply(apply(abs(PCA),2,rank),1,mean)


LCOL=vector.lcol(pbmc@meta.data$type)
pdf('f1.pdf',width=3.5,height=4)
plot(VEC,col=LCOL,pch=16,cex=0.2)
dev.off()



VCOL=vector.vcol(pbmc@meta.data$score, c(min(pbmc@meta.data$score),quantile(pbmc@meta.data$score,0.1),
                 median(pbmc@meta.data$score), quantile(pbmc@meta.data$score,0.9), max(pbmc@meta.data$score)),
                 c('blue3','blue3','grey95','red3','red3'))
pdf('f2.pdf',width=3.5,height=4)
plot(VEC,col=VCOL,pch=16,cex=0.2)
dev.off()




VECTOR(VEC, PCA, N=20)


