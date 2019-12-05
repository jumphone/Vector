
setwd('F:/Vector/data/MouseOligo_GSE75330')

##########################
# Read raw data
DATA=read.table(file='GSE75330_Marques_et_al_mol_counts2.tab',sep='\t',header=TRUE,row.names=1)
saveRDS(DATA,file='DATA.RDS')

###############################
# Reload data and cell labels
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseOligo_GSE75330')
LABEL=readRDS('GSE75330.LABEL.RDS')
DATA=readRDS('DATA.RDS')

############################################
# Analyze all cells
pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)
pbmc@meta.data$type=readRDS('GSE75330.LABEL.RDS')
pbmc@meta.data$type[which(pbmc@meta.data$type=='Differentiation-committed oligodendrocyte precursors')]='Differentiation-committed OPC'


###############################################
# Remove non-oligodendrocyte lineage cells
used_cell_index=which( pbmc@meta.data$type != 'PPR')
DATA=DATA[,used_cell_index]
TYPE=pbmc@meta.data$type[used_cell_index]

############################################
# Analyze oligodendrocyte lineage cells
pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)
pbmc@meta.data$type=TYPE
##########################

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
##########################

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs = 300)
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




OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)

RPPCA=vector.RPPCA(PCA)

OUT=vector.getValue(OUT, RPPCA[,1:50], SHOW=TRUE)
OUT=vector.gridValue(OUT,SHOW=TRUE)

OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)


OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)

eigs <- pbmc@reductions$pca@stdev^2



N=200
X=1:N
fit=lm(eigs[1:N]~X+I(X^2))




sum(eigs[1:N])/sum(VAR)
(sum(eigs[1:N])+(eigs[N]+0)*(length(VariableFeatures(pbmc))-N)/2)






OUT=vector.selectCenter(OUT)
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






GTEX=read.table('F:/Vector/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct',header=TRUE,row.names=1,sep='\t')
GTEX=GTEX[GTEX[,1] %in% names(which(table(GTEX[,1])==1)),]
RN=GTEX[,1]
GTEX=GTEX[,c(2:ncol(GTEX))]
rownames(GTEX)=RN
colnames(GTEX)=paste0('GTEX_',1:ncol(GTEX))


saveRDS(GTEX,file='F:/Vector/data/GTEX.RDS')



vector.scoreGTEX <- function(pbmc){
    pbmc=pbmc
    D1=as.matrix(pbmc@assays$RNA@data)
    rownames(D1)=toupper(rownames(D1))
    D1=D1[which(rownames(D1) %in% names(which(table(rownames(D1))==1))),]
    colnames(D1)=paste0('D1_',colnames(D1))
    DGTEX=.simple_combine(D1, GTEX)$combine
    pbmc <- CreateSeuratObject(counts = DGTEX, project = "pbmc3k", min.cells = 0, min.features = 0)
    pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
    pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
    pbmc <- ScaleData(pbmc, features = VariableFeatures(pbmc))
    pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs = 50)
    D1_INDEX=which(pbmc@meta.data$orig.ident=='D1')
    D2_INDEX=which(pbmc@meta.data$orig.ident=='GTEX')
    PCA=pbmc@reductions$pca@cell.embeddings
    PCA.N=apply(PCA,2,rank)
    MEAN=apply(PCA[D2_INDEX,],2,mean)
    MEAN.MAT=matrix(rep(MEAN,each=nrow(PCA)),nrow=nrow(PCA))
    SCORE=apply(apply(abs(PCA-MEAN.MAT),2,rank),1,mean)
    SCORE=SCORE[D1_INDEX]
    return(SCORE)
    }






SD=apply(PCA[D2_INDEX,],2,sd)

SD.MAT=matrix(rep(SD,each=nrow(PCA)),nrow=nrow(PCA))

PCA.N=(PCA-MEAN.MAT)/SD.MAT



PCA.N=()

SCORE=apply(apply(abs(PCA.N),2,rank),1,mean)


SCORE1=apply(apply(abs(PCA),2,rank),1,mean)


pbmc1=readRDS('pbmc.RDS')
pbmc@meta.data$orig.ident=as.character(pbmc@meta.data$orig.ident)
pbmc@meta.data$orig.ident[D1_INDEX]=as.character(pbmc1@meta.data$type)




DimPlot(pbmc,reductions='pca',dims=c(1,2),group.by='orig.ident')

PCA=pbmc@reductions$pca@cell.embeddings[1:ncol(D1),]











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















