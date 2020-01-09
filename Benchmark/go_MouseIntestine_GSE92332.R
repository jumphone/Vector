
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
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseIntestine_GSE92332/')
pbmc=readRDS( file='pbmc.RDS')



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




#######################################################################
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





############################################################################
# All PCs
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseIntestine_GSE92332/')

pbmc=readRDS('pbmc.RDS')
#pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs = 500)
#saveRDS(pbmc,file='pbmc_500.RDS')

library(gmodels)
D=as.matrix(pbmc@assays$RNA@scale.data)
D=D[which(rownames(D) %in% VariableFeatures(pbmc)),]
PCA.OUT=fast.prcomp(t(D), retx = TRUE, center = FALSE, scale. = FALSE, tol = NULL)
PCA=PCA.OUT$x

saveRDS(PCA.OUT,file='PCA.OUT_500.RDS')


####
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseIntestine_GSE92332/')
PCA.OUT=readRDS(file='PCA.OUT_500.RDS')
pbmc=readRDS('pbmc.RDS')

PCA=PCA.OUT$x

#################################

tiff(paste0("IMG/VAR_EXP.tiff"),width=3,height=2.5,units='in',res=600)

par(mar=c(2,2,2,2))
EXP.VAR=cumsum(PCA.OUT$sdev^2)/sum(PCA.OUT$sdev^2)
INDEX=c(1:ncol(PCA))
INDEX=log(INDEX,2)

plot( INDEX, EXP.VAR,type='l',pch=16, lwd=5, col='grey70')

THIS=150
segments(x0=INDEX[THIS],
	 y0=0,
	 x1=INDEX[THIS],
	 y1=EXP.VAR[THIS],
	 col='black',lwd=2,lty=5)
segments(x0=-10,
	 y0=EXP.VAR[THIS],
	 x1=INDEX[THIS],
	 y1=EXP.VAR[THIS],
	 col='black',lwd=2,lty=5)
dev.off()



########################

VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA=PCA #pbmc@reductions$pca@cell.embeddings



TAG='Stem'
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

















#####################
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseIntestine_GSE92332/')

pbmc=readRDS(file='pbmc.RDS')


VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings
#######################


OUT=vector.buildGrid(VEC, N=15,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.75,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=1,SHOW=TRUE,  COL=OUT$COL,AL=40)


tiff(paste0("IMG/NEW_VECTOR_TRY.6.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.drawArrow(OUT,P=1,SHOW=TRUE, COL=OUT$COL,AL=60,OL=1.5,AW=2,AC='black',BD=FALSE)
dev.off()



tiff(paste0("IMG/Marker.big.tiff"),width=4,height=2,units='in',res=600)
p <- FeaturePlot(pbmc, features=c('Lgr5','Pcna'),order=TRUE, combine = FALSE)

for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}

cowplot::plot_grid(plotlist = p)
dev.off()



tiff(paste0("IMG/Marker.big.new.tiff"),width=4,height=2,units='in',res=600)
p <- FeaturePlot(pbmc, features=c('Lgr5','Fabp1'),order=TRUE, combine = FALSE)

for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}

cowplot::plot_grid(plotlist = p)
dev.off()






tiff(paste0("IMG/NEW_VECTOR.6.tiff"),width=1,height=1,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.drawArrow(OUT,P=1,SHOW=TRUE, COL=OUT$COL,AL=20,CEX=0.2,SHOW.SUMMIT=TRUE)
dev.off()


#################

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


TTT=as.factor(pbmc@meta.data$type)
TTT=factor(TTT , levels=c("STEM",'TA','EP','EIM','EM'))

tiff(paste0("IMG/CHANGE.tiff"),width=3,height=1,units='in',res=600)
par(mar=c(0,0,0,0))
boxplot(OUT$VALUE~TTT,outline=FALSE,xlab='',ylab='',las=2)
dev.off()


tiff(paste0("IMG/CHANGE1.tiff"),width=3,height=1,units='in',res=600)
#par(mar=c(0,0,0,0))
par(mar=c(0,3,0,0))
boxplot(OUT$VALUE~TTT,outline=FALSE,xlab='',ylab='',las=2)
dev.off()


tiff(paste0("IMG/CHANGE2.tiff"),width=4.2,height=1.4,units='in',res=600)
par(mar=c(1,1,0.2,0.2))
boxplot(OUT$VALUE~TTT,outline=FALSE,xlab='',ylab='',las=2)
dev.off()


###################



OUT=vector.buildGrid(VEC, N=15,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.75,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=1,SHOW=TRUE,  COL=OUT$COL,AL=40)

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



tiff(paste0("IMG/VECTOR.6.big.tiff"),width=2,height=2,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.drawArrow(OUT,P=1,SHOW=TRUE, COL=OUT$COL,AL=20)
dev.off()




tiff(paste0("IMG/Marker.big.tiff"),width=2,height=2,units='in',res=600)
p <- FeaturePlot(pbmc, features=c('Lgr5'),order=TRUE, combine = FALSE)

for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}

cowplot::plot_grid(plotlist = p)
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




OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)


tiff(paste0("IMG/NEW_SMART.VECTOR_TRY.6.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL,AL=60,AW=1,AC='black',BD=FALSE)
dev.off()



tiff(paste0("IMG/NEW_SMART.VECTOR_TRY_BLUE.6.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL,AL=60,AW=1,AC='black',BD=FALSE)
dev.off()


tiff(paste0("IMG/NEW_SMART.VECTOR.6.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL,AL=40)
dev.off()




tiff(paste0("IMG/SMART.Marker.tiff"),width=5,height=5,units='in',res=600)
p <- FeaturePlot(pbmc, features=c('Lgr5','Pcna','Muc2','Dclk1'),order=TRUE, combine = FALSE,pt.size=0.5)

for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}

cowplot::plot_grid(plotlist = p)
dev.off()



tiff(paste0("IMG/SMART.Marker.new.tiff"),width=5,height=5,units='in',res=600)
p <- FeaturePlot(pbmc, features=c('Lgr5','Fabp1','Muc2','Dclk1'),order=TRUE, combine = FALSE,pt.size=0.5)

for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}

cowplot::plot_grid(plotlist = p)
dev.off()




tiff(paste0("IMG/SMART.Marker.new2.tiff"),width=5,height=5,units='in',res=600)
p <- FeaturePlot(pbmc, features=c('Lgr5','Fabp6','Muc2','Dclk1'),order=TRUE, combine = FALSE,pt.size=0.5)

for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}

cowplot::plot_grid(plotlist = p)
dev.off()




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




tiff(paste0("IMG/SMART.Marker.tiff"),width=5,height=5,units='in',res=600)
p <- FeaturePlot(pbmc, features=c('Lgr5','Pcna','Muc2','Dclk1'),order=TRUE, combine = FALSE)

for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}

cowplot::plot_grid(plotlist = p)
dev.off()






##########################3
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseIntestine_GSE92332/')
pbmc=readRDS(file='pbmc_smart.RDS')


VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings

OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)



###########################################
VALUE=list()
#######################
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=FALSE)
VALUE$msPCA=OUT$VALUE
############################


pdf('./IMG/TEST_nFeature_RNA.pdf')
OUT$VALUE=pbmc@meta.data$nFeature_RNA
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=TRUE)
dev.off()
VALUE$nGene=OUT$VALUE



pdf('./IMG/TEST_MS_VGENE.pdf')
VGENE=t(as.matrix(pbmc@assays$RNA@data[which(rownames(pbmc) %in% VariableFeatures(pbmc)),]))
OUT=vector.getValue(OUT, VGENE, SHOW=TRUE)
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=TRUE)
dev.off()
VALUE$msGene=OUT$VALUE




pdf('./IMG/TEST_VAR_PCA.pdf')
OUT$VALUE=apply(PCA,1,var)
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=TRUE)
dev.off()
VALUE$varPCA=OUT$VALUE



saveRDS(VALUE,'IMG/TEST_VALUE.RDS')

###############
MAT <- matrix(unlist(VALUE), ncol = 4, byrow = FALSE)
colnames(MAT)=names(VALUE)
rownames(MAT)=colnames(pbmc)
saveRDS(MAT,'IMG/TEST_MAT.RDS')
####################

COR=cor(MAT, method='spearman')

COR
#            msPCA      nGene     msGene     varPCA
#msPCA   1.0000000 -0.4431614  0.6708386  0.7994030
#nGene  -0.4431614  1.0000000 -0.9102853 -0.3923117
#msGene  0.6708386 -0.9102853  1.0000000  0.5461314
#varPCA  0.7994030 -0.3923117  0.5461314  1.0000000






######################################################






###
# Test pathway feature

##########################
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseIntestine_GSE92332/')
pbmc=readRDS(file='pbmc_smart.RDS')


VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings

OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)



###########################################
VALUE=list()
#######################
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=FALSE)
VALUE$msPCA=OUT$VALUE
############################


###########
#Prepare EXP

EXP=as.matrix(pbmc@assays$RNA@scale.data[which(rownames(pbmc@assays$RNA@scale.data) %in% VariableFeatures(pbmc)),])
#EXP=t(apply(t(EXP),2,scale))
#EXP[which(is.na(EXP))]=0
#colnames(EXP)=colnames(pbmc)
UP=toupper(rownames(EXP))
USED=which(UP %in% names(which(table(UP)==1)))
EXP=EXP[USED,]
rownames(EXP)=UP[USED]


###########
#KEGG
library(qusage)
BIO=read.gmt('../KEGG_2019_Mouse.gmt')
BIO.MAT=matrix(0,nrow=length(BIO),ncol=ncol(EXP))
colnames(BIO.MAT)=colnames(EXP)
rownames(BIO.MAT)=names(BIO)

i=1
while(i<=length(BIO)){
    used=which(rownames(EXP) %in% BIO[[i]])
    if(length(used)>1){
        BIO.MAT[i,]=apply(EXP[used,],2,mean) 
    }else if(length(used)==1){
	BIO.MAT[i,]=EXP[used,]
    }else{
        BIO.MAT[i,]=0}
    
    i=i+1}

VAR=apply(BIO.MAT,1,var)
BIO.MAT=BIO.MAT[which(VAR>0 & (!is.na(VAR))),]
this_value=vector.calValue(t(BIO.MAT))$VALUE
####################################

pdf('./IMG/TEST_msKEGG.pdf')
OUT$VALUE=this_value
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=TRUE)
dev.off()
VALUE$msKEGG=OUT$VALUE

##########################################################


###########
#MSigDB
library(qusage)
BIO=read.gmt('../MSigDB_Computational.gmt')
BIO.MAT=matrix(0,nrow=length(BIO),ncol=ncol(EXP))
colnames(BIO.MAT)=colnames(EXP)
rownames(BIO.MAT)=names(BIO)

i=1
while(i<=length(BIO)){
    used=which(rownames(EXP) %in% BIO[[i]])
    if(length(used)>1){
        BIO.MAT[i,]=apply(EXP[used,],2,mean) 
    }else if(length(used)==1){
	BIO.MAT[i,]=EXP[used,]
    }else{
        BIO.MAT[i,]=0}
    
    i=i+1}

VAR=apply(BIO.MAT,1,var)
BIO.MAT=BIO.MAT[which(VAR>0 & (!is.na(VAR))),]
this_value=vector.calValue(t(BIO.MAT))$VALUE
####################################

pdf('./IMG/TEST_msMSigDB.pdf')
OUT$VALUE=this_value
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=TRUE)
dev.off()
VALUE$msMSigDB=OUT$VALUE

##########################################################



MAT <- matrix(unlist(VALUE), ncol = 3, byrow = FALSE)
colnames(MAT)=names(VALUE)
rownames(MAT)=colnames(pbmc)
saveRDS(MAT,'IMG/TEST_BIO_MAT.RDS')


########################################






######################
# Test random feature


##########################
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseIntestine_GSE92332/')
pbmc=readRDS(file='pbmc_smart.RDS')


VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings

OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)



###########################################
VALUE=list()
#######################
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=FALSE)
VALUE$msPCA=OUT$VALUE
############################


###########

set.seed(123)

pdf('./IMG/TEST_random.pdf')
OUT$VALUE=rnorm(length(OUT$VALUE))
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=TRUE)
dev.off()
VALUE$random=OUT$VALUE

saveRDS(MAT,'IMG/TEST_RANDOM_MAT.RDS')







































