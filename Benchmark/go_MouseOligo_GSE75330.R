
setwd('F:/Vector/data/MouseOligo_GSE75330')
DATA=read.table(file='GSE75330_Marques_et_al_mol_counts2.tab',sep='\t',header=TRUE,row.names=1)
saveRDS(DATA,file='DATA.RDS')

L=readRDS('LABEL.RDS')
L=L[(length(L)-5069+1):length(L)]
library(stringr)
L=str_replace(L, '_batch2', '')
saveRDS(L,file='GSE75330.LABEL.RDS')


############################################
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseOligo_GSE75330')
DATA=readRDS('DATA.RDS')
pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)
pbmc@meta.data$type=readRDS('GSE75330.LABEL.RDS')
pbmc@meta.data$type[which(pbmc@meta.data$type=='Differentiation-committed oligodendrocyte precursors')]='Differentiation-committed OPC'


############################################
used_cell_index=which( pbmc@meta.data$type != 'PPR')

DATA=DATA[,used_cell_index]
TYPE=pbmc@meta.data$type[used_cell_index]
pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)
pbmc@meta.data$type=TYPE
##########################

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
##########################

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs = 50)
pbmc <- RunUMAP(pbmc, dims = 1:50)
DimPlot(pbmc, reduction = "umap")
saveRDS(pbmc,file='pbmc.RDS')

###################



#############################
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

setwd('F:/Vector/data/MouseOligo_GSE75330')
pbmc=readRDS('pbmc.RDS')

FeaturePlot(pbmc,features=c('Pdgfra','Bmp4','Sema4f','Mog'))



DimPlot(pbmc, group.by='type',label = TRUE,
        label.size=5,pt.size=0.5)+NoLegend()

###################################

VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings

DATA=as.matrix(pbmc@assays$RNA@data)


OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)

SCORE=vector.calScore(OUT,PCA,SHOW=TRUE)
#VALUE=SCORE#vector.getScore(PCA)

OUT=vector.gridValue(OUT,SCORE, SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
#OUT=vector.autoCenterNew(OUT,SHOW=TRUE)
#OUT=vector.selectCenter(OUT)
#OUT=vector.nonCenter(OUT)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)

OUT=vector.nonCenter(OUT)



OUT=vector.drawArrow(OUT,P=1,SHOW=TRUE, COL=OUT$COL)










OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)

SCORE=vector.getValue(PCA)
VALUE=SCORE#vector.getScore(PCA)

OUT=vector.gridValue(OUT,VALUE, SHOW=TRUE)

OUT=vector.autoCenterNew(OUT,SHOW=TRUE)
#OUT=vector.selectCenter(OUT)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)














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


###########################################################















