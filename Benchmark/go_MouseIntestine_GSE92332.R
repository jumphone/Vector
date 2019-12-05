
########################################
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseIntestine_GSE92332/')

#######################################
# Load data of GSE92332
DATA=readRDS('DATA.RDS')
TYPE=readRDS('LABEL.RDS')
pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)
BATCH=pbmc@meta.data$orig.ident


used_cell_index=which( TYPE %in% c('Enterocyte.Immature.Distal',
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
# Select enterocyte lineage cells

DATA=DATA[,used_cell_index]
TYPE=TYPE[used_cell_index]
BATCH=BATCH[used_cell_index]
####################################################


mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=150, ROUND=1, GN=2000, SEED=1, COMBAT=TRUE, RMG=NULL) 
saveRDS(mybeer, file='mybeer.RDS')

#mybeer=readRDS(file='mybeer.RDS')
pbmc <- mybeer$seurat
PCUSE=mybeer$select   
pbmc=BEER.combat(pbmc)
umap=BEER.bbknn(pbmc, PCUSE, NB=3, NT=10)

pbmc@reductions$umap@cell.embeddings=umap
#DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1,label=F)
pbmc@meta.data$type=TYPE
DimPlot(pbmc, reduction = "umap",group.by='type')
saveRDS(pbmc, file='pbmc.RDS')
saveRDS(PCUSE, file='PCUSE.RDS')

#####################################################

VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings

pbmc@meta.data$type[which(pbmc@meta.data$type %in% c('Enterocyte.Immature.Distal',
                                 'Enterocyte.Immature.Proximal') )]='EIM'
pbmc@meta.data$type[which(pbmc@meta.data$type %in% c('Enterocyte.Mature.Distal',
                                 'Enterocyte.Mature.Proximal') )]='EM'
pbmc@meta.data$type[which(pbmc@meta.data$type %in% c('Enterocyte.Progenitor',
                                  'Enterocyte.Progenitor.Early',
                                  'Enterocyte.Progenitor.Late') )]='EP'
pbmc@meta.data$type[which(pbmc@meta.data$type %in% c('TA.Early',
                                  'TA.G1','TA.G1') )]='TA'
pbmc@meta.data$type[which(pbmc@meta.data$type %in% c('Stem') )]='STEM'


length(pbmc@meta.data$type)
#[1] 5560

table(pbmc@meta.data$type)

# EIM   EM   EP STEM   TA 
# 809  822 1589 1267 1073 

###########################################
#Draw heatmap
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseIntestine_GSE92332/')

pbmc=readRDS(file='pbmc.RDS')


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

table(pbmc@meta.data$type)

# EIM   EM   EP STEM   TA 
# 809  822 1589 1267 1073 

VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings

pbmc@meta.data$type[which(pbmc@meta.data$type %in% c('Enterocyte.Immature.Distal',
                                 'Enterocyte.Immature.Proximal') )]='EIM'
pbmc@meta.data$type[which(pbmc@meta.data$type %in% c('Enterocyte.Mature.Distal',
                                 'Enterocyte.Mature.Proximal') )]='EM'
pbmc@meta.data$type[which(pbmc@meta.data$type %in% c('Enterocyte.Progenitor',
                                  'Enterocyte.Progenitor.Early',
                                  'Enterocyte.Progenitor.Late') )]='EP'
pbmc@meta.data$type[which(pbmc@meta.data$type %in% c('TA.Early',
                                  'TA.G1','TA.G1') )]='TA'
pbmc@meta.data$type[which(pbmc@meta.data$type %in% c('Stem') )]='STEM'


#########################################
TAG='STEM'
DrawHeatMap(TAG)
###################################################
TAG='TA'
DrawHeatMap(TAG)
###################################################
TAG='EP'
DrawHeatMap(TAG)
###################################################
TAG='EIM'
DrawHeatMap(TAG)
###################################################
TAG='EM'
DrawHeatMap(TAG)



########################################################
# Draw rank curve


TAG='STEM'
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

tiff(paste0("IMG/STEM.SCORE.tiff"),width=3,height=3,units='in',res=600)

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


#######################################################################################
#Draw heatmap
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseIntestine_GSE92332/')

pbmc=readRDS(file='pbmc.RDS')




VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings
#######################



tiff(paste0("IMG/VECTOR.1.tiff"),width=1,height=1,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.buildGrid(VEC, N=8,SHOW=TRUE)
dev.off()


tiff(paste0("IMG/VECTOR.2.tiff"),width=1,height=1,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
dev.off()

tiff(paste0("IMG/VECTOR.3.tiff"),width=1,height=1,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
dev.off()

tiff(paste0("IMG/VECTOR.4.tiff"),width=1,height=1,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.gridValue(OUT,SHOW=TRUE)
dev.off()

tiff(paste0("IMG/VECTOR.5.tiff"),width=1,height=1,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.autoCenter(OUT,UP=0.65,SHOW=TRUE)
dev.off()

tiff(paste0("IMG/VECTOR.6.tiff"),width=1,height=1,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.drawArrow(OUT,P=1,SHOW=TRUE, COL=OUT$COL)
dev.off()




pbmc@meta.data$type[which(pbmc@meta.data$type %in% c('Enterocyte.Immature.Distal',
                                 'Enterocyte.Immature.Proximal') )]='EIM'
pbmc@meta.data$type[which(pbmc@meta.data$type %in% c('Enterocyte.Mature.Distal',
                                 'Enterocyte.Mature.Proximal') )]='EM'
pbmc@meta.data$type[which(pbmc@meta.data$type %in% c('Enterocyte.Progenitor',
                                  'Enterocyte.Progenitor.Early',
                                  'Enterocyte.Progenitor.Late') )]='EP'
pbmc@meta.data$type[which(pbmc@meta.data$type %in% c('TA.Early',
                                  'TA.G1','TA.G1') )]='TA'
pbmc@meta.data$type[which(pbmc@meta.data$type %in% c('Stem') )]='STEM'



tiff(paste0("IMG/TYPE.tiff"),width=3,height=3,units='in',res=600)
par(mar=c(0,0,0,0))
DimPlot(pbmc,group.by='type',label=TRUE,pt.size=0.01)+NoLegend()
dev.off()

tiff(paste0("IMG/LGR5.tiff"),width=3.5,height=3.2,units='in',res=600)
par(mar=c(0,0,0,0))
FeaturePlot(pbmc,features='Lgr5',order = TRUE)
dev.off()




##########################################################################################
# Smart-seq

source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseIntestine_GSE92332/')
DATA=read.table('GSE92332_AtlasFullLength_TPM.txt',header=TRUE,row.names=1)

TMP=strsplit(colnames(DATA),'_')
LABEL=c()
i=1
while(i<=length(TMP)){
    LABEL=c(LABEL,TMP[[i]][5])
    i=i+1}

TYPE=LABEL


pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs=150)
pbmc <- RunUMAP(pbmc, dims = 1:150)
DimPlot(pbmc, reduction = "umap")
saveRDS(pbmc,file='pbmc_smart.RDS')

#################################

source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseIntestine_GSE92332/')
pbmc=readRDS(file='pbmc_smart.RDS')


VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings




tiff(paste0("IMG/SMART.VECTOR.1.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
dev.off()


tiff(paste0("IMG/SMART.VECTOR.2.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
dev.off()

tiff(paste0("IMG/SMART.VECTOR.3.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
dev.off()

tiff(paste0("IMG/SMART.VECTOR.4.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.gridValue(OUT,SHOW=TRUE)
dev.off()

tiff(paste0("IMG/SMART.VECTOR.5.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
dev.off()

tiff(paste0("IMG/SMART.VECTOR.6.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)
dev.off()




tiff(paste0("IMG/SMART.Lgr5.tiff"),width=3.5,height=3.2,units='in',res=600)
par(mar=c(0,0,0,0))
FeaturePlot(pbmc,features='Lgr5',order = TRUE)
dev.off()





























































DATA=readRDS('REF.RDS')
saveRDS(DATA,file='DATA.RDS')





########################################
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseIntestine_GSE92332/')

#
DATA=read.table('GSE92332_AtlasFullLength_TPM.txt',header=TRUE,row.names=1)

TMP=strsplit(colnames(DATA),'_')
LABEL=c()
i=1
while(i<=length(TMP)){
    LABEL=c(LABEL,TMP[[i]][5])
    i=i+1}

TYPE=LABEL


pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs=100)
pbmc <- RunUMAP(pbmc, dims = 1:100)
DimPlot(pbmc, reduction = "umap")




#PCA.OUT=vector.SeuratRandomPCA(pbmc)


#PCA=PCA.OUT$PRED.PCA[,1:PCA.OUT$N]


FeaturePlot(pbmc,features='Lgr5')


VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings#[,1:50]#[,1:100]#[,50:100]#[,30:50]



RPPCA=vector.RPPCA(PCA)
OUT=vector.buildGrid(VEC, N=40,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, RPPCA[,1:20], SHOW=TRUE)




OUT=vector.gridValue(OUT,SHOW=TRUE)
#OUT=vector.gridValueSmooth(OUT,CUT=0.98,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)




####################################################################




########################################
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseIntestine_GSE92332/')


DATA=readRDS('DATA.RDS')
TYPE=readRDS('LABEL.RDS')
pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)
BATCH=pbmc@meta.data$orig.ident


#pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)



VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

used_cell_index=which( TYPE %in% c('Enterocyte.Immature.Distal',
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

mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=100, ROUND=1, GN=2000, SEED=1, COMBAT=TRUE, RMG=NULL) 

#mybeer=readRDS(file='mybeer.RDS')
PCUSE=mybeer$select
COL=rep('black',length(mybeer$cor))
COL[PCUSE]='red'
plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,
    xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))
#####################
pbmc <- mybeer$seurat

pbmc=RunTSNE(pbmc,)



PCUSE=mybeer$select   

#pbmc=BEER.combat(pbmc) #Adjust PCs using ComBat
umap=BEER.bbknn(pbmc, PCUSE, NB=2, NT=10)



pbmc@reductions$umap@cell.embeddings=umap
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1,label=F)
pbmc@meta.data$type=TYPE
DimPlot(pbmc, reduction = "umap",group.by='type')





pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs =200)


VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings


OUT=vector.buildGrid(VEC, N=20,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)


OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.8,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)

































OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)

SCORE=vector.getValue(PCA)


OUT=vector.gridValue(OUT,VALUE, SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
#OUT=vector.selectCenter(OUT)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)
























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
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings




OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)

SCORE=vector.getValue(PCA)
VALUE=SCORE#vector.getScore(PCA)
VALUE=LGR5
OUT=vector.gridValue(OUT,VALUE, SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
#OUT=vector.selectCenter(OUT)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)


















DimPlot(pbmc,reduction = 'pca',group.by = 'type',dims = c(3,2))

VEC=pbmc@reductions$umap@cell.embeddings
#VEC=pbmc@reductions$pca@cell.embeddings[,1:2]
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings#[,10:50]

#FeaturePlot(pbmc,features = c('Lgr5','Itgb1','Smoc2','Sox9','Bmi1','Lrig1','Ascl2'),ncol=3)

FeaturePlot(pbmc,features = c('Lgr5'))
#OUT=VECTOR(VEC, PCA, N=20)

LGR5=pbmc@assays$RNA@data[which(rownames(pbmc)=='Lgr5'),]


#########################


OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)

SCORE=vector.getValue(PCA)
VALUE=SCORE#vector.getScore(PCA)

OUT=vector.gridValue(OUT,VALUE, SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.7,SHOW=TRUE)
#OUT=vector.selectCenter(OUT)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)






ppp=DimPlot(pbmc, reduction.use='umap', pt.size=0.5)
used.cells <- CellSelector(plot = ppp)
used_cell_index=which(colnames(pbmc) %in% used.cells)
DATA=DATA[,used_cell_index]
TYPE=TYPE[used_cell_index]
BATCH=BATCH[used_cell_index]

pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs=50)
pbmc <- RunUMAP(pbmc, dims = 1:50)
DimPlot(pbmc, reduction = "umap")




VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings




OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)

SCORE=vector.getValue(PCA)
VALUE=SCORE#vector.getScore(PCA)
OUT=vector.gridValue(OUT,VALUE, SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
#OUT=vector.selectCenter(OUT)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)

#######################
#OUT=vector.reDrawArrow(OUT, COL=OUT$COL)
#OUT=vector.selectRegion(OUT)



TMP=strsplit(rownames(VEC),'_')
LABEL=c()
i=1
while(i<=length(TMP)){
    LABEL=c(LABEL,TMP[[i]][3])
    i=i+1}


pbmc@meta.data$type=LABEL
DimPlot(pbmc, reduction = "umap",group.by='orig.ident',label=TRUE)

DimPlot(pbmc, reduction = "umap",group.by='type',label=TRUE)

FeaturePlot(pbmc,features=c('Bub1','Cdk1','Nusap1','Hist1h2al','Cdt1','Orc1'),ncol=3,order=TRUE)


DimPlot(pbmc, reduction = "pca",dims=c(1,2),group.by='type')






























.normX <- function(x){
    y=(x-min(x))/(max(x)-min(x))
    return(y)
    }

.densityCenter <- function(x){
    x=x
    dx=density(x) 
    cx=dx$x[which(dx$y==max(dx$y))]
    return(cx)
    }



vector.getValueNew <- function(PCA){
    PCA=PCA
    PCA.N=apply(PCA,2,.normX)
    PCA.C=apply(PCA.N,2,.densityCenter)
    PCA.NP=(PCA.C-0.5)/abs(PCA.C-0.5)
    PCA.NP.MAT=matrix(rep(PCA.NP,each=nrow(PCA)),nrow=nrow(PCA))
    SCORE=apply(apply((PCA*PCA.NP.MAT),2,rank),1,mean)
    return(SCORE)
    }



vector.getScore <- function(PCA, CUT=0.5){
    PCA=PCA
    CUT=CUT
    ################
    PCA.N=apply(PCA,2,.normX)
    PCA.C=apply(PCA.N,2,.densityCenter)
    PCA.C.MAT=matrix(rep(PCA.C,each=nrow(PCA)),nrow=nrow(PCA))
    PCA.N.D=(PCA.N-PCA.C.MAT)^2
    COR=cor(PCA.N.D,method='spearman')
    diag(COR)=0
    COR[which(COR<CUT)]=0
    SSS=apply(COR,2,sum)
    USED.PCA=which(SSS>0)
    #print(length(USED.PCA))
    print(USED.PCA)
    SCORE=apply(PCA.N.D[,USED.PCA],1,mean)
    SCORE=1-.normX(SCORE)
    ###############
    return(SCORE)
    }





mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=100, ROUND=1, GN=5000, SEED=1, COMBAT=TRUE, RMG=NULL) 


#mybeer=readRDS(file='mybeer.RDS')
PCUSE=mybeer$select
COL=rep('black',length(mybeer$cor))
COL[PCUSE]='red'
plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,
    xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))
DimPlot(pbmc, reduction = "umap",group.by='type',label=TRUE)

DimPlot(pbmc,reduction = 'pca',group.by = 'type',dims = c(1,2))

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


#PCA.NP=-(PCA.C-0.5)/abs(PCA.C-0.5)

#PCA.NP.MAT=matrix(rep(PCA.NP,each=nrow(PCA)),nrow=nrow(PCA))


#USED=which(abs(PCA.C-0.5)>0.1)

SCORE=vector.getValueNew(PCA)        #apply(apply((PCA*PCA.NP.MAT),2,rank),1,mean)

#SCORE=rank(PCA[,1])
#PCA.N=apply(PCA,2,.normX)
#    PCA.C=apply(PCA.N,2,.densityCenter)



#OUT=VECTOR(VEC, PCA, N=20)

plot(SCORE,LGR5)



OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)

SCORE=vector.getValue(PCA)
VALUE=SCORE#vector.getScore(PCA)
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


