
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

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs = 150)
pbmc <- RunUMAP(pbmc, dims = 1:50)
DimPlot(pbmc, reduction = "umap")
saveRDS(pbmc,file='pbmc.RDS')
######################################

#############################
# Reload oligodendrocyte lineage cells

source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseOligo_GSE75330')
pbmc=readRDS('pbmc.RDS')
#FeaturePlot(pbmc,features=c('Pdgfra','Bmp4','Sema4f','Mog'))

VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings


###########################################
#Draw heatmap

DrawHeatMap<-function(TAG){
    TAG=TAG
    R.PCA=apply(PCA,2,rank)
    THIS.INDEX=which(pbmc@meta.data$type==TAG)
    THIS.R.PCA=R.PCA[THIS.INDEX,]


    MAT=matrix(0,nrow=ncol(PCA),ncol=nrow(PCA))
    rownames(MAT)=colnames(PCA)
    colnames(MAT)=paste0('R_',1:nrow(PCA))
    print(MAT[1:3,1:3])
	
    i=1
    while(i<=ncol(PCA)){
        MAT[i,THIS.R.PCA[,i]]=1
        i=i+1}
    #######################
    #MAT=t(apply(t(MAT),2,smooth))
    ########################
	
    library('ComplexHeatmap')
    library('circlize')
    library('seriation')
    
    mat=MAT
    o.mat=mat
    col_fun =colorRamp2(c(0,1), c('white','#000080'))
     
    LLL=apply(mat,2,mean)
    ha = HeatmapAnnotation(
	M = anno_lines(LLL, add_points = FALSE,smooth=FALSE,
		       gp = gpar(col = 'grey50',lwd=0.5),
		      axis=FALSE),
	name=c(''),
	show_annotation_name=FALSE
        )	
	
    tiff(paste0("IMG/",TAG,".HEAT.tiff"),width=4,height=1,units='in',res=600)
    draw(Heatmap(o.mat,row_title='',name="",cluster_rows=FALSE,
        cluster_columns=FALSE,show_heatmap_legend=FALSE，
	show_column_dend = FALSE, show_row_dend = FALSE, 
	show_column_names=FALSE, show_row_names=FALSE,
	col=col_fun, border = TRUE,
        top_annotation = ha
	))
    dev.off()
    }

#########################################
length(TYPE)
#4993

table(TYPE)
#Differentiation-committed OPC         Mature Oligodendrocytes 
#                            140                            2748 
#Myelin-forming Oligodendrocytes   Newly-formed Oligodendrocytes 
#                           1283                             512 
#                            OPC 
#                            310 

#########################################
TAG='OPC'
DrawHeatMap(TAG)
###################################################
TAG='Differentiation-committed OPC'
DrawHeatMap(TAG)
###################################################
TAG='Newly-formed Oligodendrocytes'
DrawHeatMap(TAG)
###################################################
TAG='Myelin-forming Oligodendrocytes'
DrawHeatMap(TAG)
###################################################
TAG='Mature Oligodendrocytes'
DrawHeatMap(TAG)


########################################################
# Draw rank curve

source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseOligo_GSE75330')
pbmc=readRDS('pbmc.RDS')
#FeaturePlot(pbmc,features=c('Pdgfra','Bmp4','Sema4f','Mog'))

VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings



TAG='OPC'
R.PCA=apply(PCA,2,rank)
THIS.INDEX=which(pbmc@meta.data$type==TAG)
#THIS.R.PCA=R.PCA[THIS.INDEX,]

MEAN=c()
UP=c()
LW=c()

N=1
PCA.RC=.normX(rank(PCA[,1]))
PCA.RC=abs(PCA.RC-0.5)   
VALUE=PCA.RC
R.VALUE=rank(VALUE)/length(VALUE)
this_mean=median(R.VALUE[THIS.INDEX])
this_up=quantile(R.VALUE[THIS.INDEX],0.975)
this_lw=quantile(R.VALUE[THIS.INDEX],0.025)
MEAN=c(MEAN, this_mean)
UP=c(UP,this_up)
LW=c(LW,this_lw)


N=2
while(N<=150){
    PCA.RC=apply(apply(PCA[,1:N],2,rank), 2, .normX)
    PCA.RC=abs(PCA.RC-0.5)   
    VALUE=apply(PCA.RC,1,mean)
    R.VALUE=rank(VALUE)/length(VALUE)
    this_mean=median(R.VALUE[THIS.INDEX])
    this_up=quantile(R.VALUE[THIS.INDEX],0.975)
    this_lw=quantile(R.VALUE[THIS.INDEX],0.025)
    MEAN=c(MEAN, this_mean)
    UP=c(UP,this_up)
    LW=c(LW,this_lw)
    print(N)
    N=N+1
    }

tiff(paste0("IMG/OPC.SCORE.tiff"),width=3,height=3,units='in',res=600)

par(mar=c(2,2,2,2))
plot(MEAN,type='l',lwd=1,ylim=c(0,1))
points(LW,type='l',pch=16,cex=1,col='grey70')
points(UP,type='l',pch=16,cex=1, col='grey70')

segments(x0=c(1:length(MEAN)),
	 y0=LW,
	 x1=c(1:length(MEAN)),
	 y1=UP,
	 col='grey70',lwd=0.5)
abline(h=0.5,lty=2, lwd=2)
points(MEAN,type='l',lwd=2)

dev.off()






############################################################################
# PC500

source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseOligo_GSE75330')
pbmc=readRDS('pbmc.RDS')
#pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs = 500)
#saveRDS(pbmc,file='pbmc_500.RDS')

library(gmodels)
D=as.matrix(pbmc@assays$RNA@scale.data)
D=D[which(rownames(D) %in% VariableFeatures(pbmc)),]
PCA.OUT=fast.prcomp(t(D), retx = TRUE, center = FALSE, scale. = FALSE, tol = NULL)
PCA=PCA.OUT$x

saveRDS(PCA.OUT,file='PCA.OUT_500.RDS')


#####

####

source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseOligo_GSE75330')
PCA.OUT=readRDS(file='PCA.OUT_500.RDS')
pbmc=readRDS('pbmc.RDS')

PCA=PCA.OUT$x

VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA=PCA #pbmc@reductions$pca@cell.embeddings



TAG='OPC'
R.PCA=apply(PCA,2,rank)
THIS.INDEX=which(pbmc@meta.data$type==TAG)
#THIS.R.PCA=R.PCA[THIS.INDEX,]

MEAN=c()
UP=c()
LW=c()

N=1
PCA.RC=.normX(rank(PCA[,1]))
PCA.RC=abs(PCA.RC-0.5)   
VALUE=PCA.RC
R.VALUE=rank(VALUE)/length(VALUE)
this_mean=median(R.VALUE[THIS.INDEX])
this_up=quantile(R.VALUE[THIS.INDEX],0.975)
this_lw=quantile(R.VALUE[THIS.INDEX],0.025)
MEAN=c(MEAN, this_mean)
UP=c(UP,this_up)
LW=c(LW,this_lw)


N=2
while(N<=nrow(PCA)){
    PCA.RC=.normX(rank(PCA[,N]))
    PCA.RC=abs(PCA.RC-0.5)      
    VALUE = VALUE+PCA.RC
    R.VALUE=rank(VALUE)/length(VALUE)
    this_mean=median(R.VALUE[THIS.INDEX])
    this_up=quantile(R.VALUE[THIS.INDEX],0.975)
    this_lw=quantile(R.VALUE[THIS.INDEX],0.025)
    MEAN=c(MEAN, this_mean)
    UP=c(UP,this_up)
    LW=c(LW,this_lw)
    print(N)
    N=N+1
    }


tiff(paste0("IMG/OPC.SCORE.500.tiff"),width=3,height=3,units='in',res=600)


INDEX=c(1:ncol(PCA))
INDEX=log(INDEX,2)
par(mar=c(2,2,2,2))
plot( INDEX,MEAN,type='l',lwd=1,ylim=c(0,1))
points( INDEX, LW,type='p',pch=16,cex=0.5,col='grey70')
points( INDEX, UP,type='p',pch=16,cex=0.5, col='grey70')

segments(x0=INDEX,
	 y0=LW,
	 x1=INDEX,
	 y1=UP,
	 col='grey70',lwd=0.5)
abline(h=0.5,lty=2, lwd=2)
abline(v=INDEX[150],lty=2, lwd=2,col='black')
points(INDEX,MEAN,type='l',lwd=2)

EXP.VAR=cumsum(PCA.OUT$sdev^2)/sum(PCA.OUT$sdev^2)
points( INDEX, EXP.VAR,type='l',pch=16, lwd=2, col='red3')


dev.off()

########################################







###################################

source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseOligo_GSE75330')
pbmc=readRDS('pbmc.RDS')




VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings

OUT=vector.buildGrid(VEC, N=20,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)


tiff(paste0("IMG/VECTOR.1.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.buildGrid(VEC, N=20,SHOW=TRUE)
dev.off()


tiff(paste0("IMG/VECTOR.2.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
dev.off()

tiff(paste0("IMG/VECTOR.3.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
dev.off()

tiff(paste0("IMG/VECTOR.4.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.gridValue(OUT,SHOW=TRUE)
dev.off()

tiff(paste0("IMG/VECTOR.5.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
dev.off()


tiff(paste0("IMG/VECTOR.6.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)
dev.off()



tiff(paste0("IMG/NEW_VECTOR.6.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)
dev.off()




tiff(paste0("IMG/Pdgfra.tiff"),width=3.5,height=3.2,units='in',res=600)
par(mar=c(0,0,0,0))
FeaturePlot(pbmc,features='Pdgfra',order = TRUE)
dev.off()




tiff(paste0("IMG/Marker.tiff"),width=5,height=5,units='in',res=600)
p <- FeaturePlot(pbmc, features=c('Pdgfra','Bmp4','Erbb3','Mog'),order=TRUE, combine = FALSE)

for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}

cowplot::plot_grid(plotlist = p)
dev.off()


