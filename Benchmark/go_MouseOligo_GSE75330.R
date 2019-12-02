
setwd('F:/Vector/data/MouseOligo_GSE75330')
DATA=read.table(file='GSE75330_Marques_et_al_mol_counts2.tab',sep='\t',header=TRUE,row.names=1)
saveRDS(DATA,file='DATA.RDS')


############################################

setwd('F:/Vector/data/MouseOligo_GSE75330')
DATA=readRDS('DATA.RDS')
pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs = 50)
pbmc <- RunUMAP(pbmc, dims = 1:50)
DimPlot(pbmc, reduction = "umap")
saveRDS(pbmc,file='pbmc.RDS')


############################################
L=readRDS('LABEL.RDS')
L=L[(length(L)-5069+1):length(L)]
library(stringr)
L=str_replace(L, '_batch2', '')
saveRDS(L,file='GSE75330.LABEL.RDS')


#############################
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

setwd('F:/Vector/data/MouseOligo_GSE75330')
pbmc=readRDS('pbmc.RDS')

pbmc@meta.data$type=readRDS('GSE75330.LABEL.RDS')
pbmc@meta.data$type[which(pbmc@meta.data$type=='Differentiation-committed oligodendrocyte precursors')]='Differentiation-committed OPC'


DimPlot(pbmc, group.by='type',label = TRUE,
        label.size=5,pt.size=0.5)+NoLegend()
FeaturePlot(pbmc,features=c('Pdgfra','Bmp4','Sema4f','Mog'))

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



VCOL=vector.vcol(pbmc@meta.data$score, c(min(pbmc@meta.data$score),
                                         median(pbmc@meta.data$score),
                                         max(pbmc@meta.data$score)),c('blue3','grey80','red3'))
pdf('f2.pdf',width=3.5,height=4)
plot(VEC,col=VCOL,pch=16,cex=0.2)
dev.off()


###########################################################





























VlnPlot(pbmc,features = 'score',group.by='celltype')
Idents(pbmc)=pbmc@meta.data$celltype
VlnPlot(pbmc,features = 'score',#group.by='celltype',
       sort=TRUE, 
       pt.size=0.5,
       idents=c('Differentiation-committed OPC',
                         'OPC',
                         'Mature Oligodendrocytes',
                         'Myelin-forming Oligodendrocytes',
                         'Newly-formed Oligodendrocytes'))+NoLegend()





DimPlot(pbmc,label = TRUE,
        label.size=5,pt.size=0.5)+NoLegend()







