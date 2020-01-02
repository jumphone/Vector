source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/LUOZAILI')

pbmc=readRDS('12W_pbmc3k_final.rds')


pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs = 150)



VEC = pbmc@reductions$umap@cell.embeddings
rownames(VEC) = colnames(pbmc)
PCA = pbmc@reductions$pca@cell.embeddings

# Define pixel
OUT=vector.buildGrid(VEC, N=40,SHOW=TRUE)

# Build network
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)

# Calculate Margin Score (MS)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)

# Get pixel's MS
OUT=vector.gridValue(OUT,SHOW=TRUE)

# Find summit
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)

# Infer vector
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=TRUE)
####################################


FeaturePlot(pbmc,features='NES')




DATA=readRDS('5_all_cells.rds')

D=DATA@assays$RNA@counts
BATCH=DATA@meta.data$batch
DATA=D
UB=table(BATCH)
FOLD=rep(5,length(UB))
names(FOLD)=names(UB)

AGG=BEER.AGG(DATA, BATCH, FOLD, PCNUM=50, GN=2000, CPU=4, print_step=10, SEED=123, N=3, RMG=NULL)
saveRDS(AGG, 'AGG.RDS')


########################
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/LUOZAILI')
AGG=readRDS('AGG.RDS')

DATA=AGG$data.agg
BATCH=AGG$data.agg.batch

mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=150, ROUND=1, GN=5000, SEED=1, COMBAT=TRUE, RMG=NULL,N=3)   

saveRDS(mybeer, 'mybeer.RDS')






############
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/LUOZAILI')

mybeer=readRDS('mybeer.RDS')


PCUSE=mybeer$select
COL=rep('black',length(mybeer$cor))
COL[PCUSE]='red'
plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,
    xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))



pbmc <- mybeer$seurat

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs = 300)



saveRDS(pbmc,file='pbmc300.RDS')






############
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/LUOZAILI')

pbmc=readRDS(file='pbmc300.RDS')

PCA = pbmc@reductions$pca@cell.embeddings


###########################
CCG=as.character(read.table('CellCycle.txt')[,1])
SD=as.matrix(pbmc@assays$RNA@scale.data)
CCG.EXP=apply(SD[which(rownames(SD) %in% CCG),],2,mean)
pbmc@meta.data$ccg=CCG.EXP

#########################
N.PCA=PCA
i=1
while(i<=ncol(PCA)){
    this_fit=loess(PCA[,i]~CCG.EXP)
    N.PCA[,i]= PCA[,i] -  predict(this_fit)	
    #this_order=order(CCG.EXP)
    #N.PCA[this_order,i]= PCA[this_order,i] -  smooth(PCA[this_order,i]) 
    print(i)
    i=i+1}


pbmc@reductions$pca@cell.embeddings=N.PCA



#####################


PCUSE=1:80



pbmc <- RunUMAP(object = pbmc, reduction='pca',dims = PCUSE, n.components=2, 
		#n.neighbors=70, min.dist=0.3,seed.use=123,
		check_duplicates=FALSE)

DimPlot(pbmc, reduction='umap', group.by='batch', pt.size=0.1) 






nUMAP=round(length(PCUSE))
pbmc <- RunUMAP(object = pbmc, reduction='pca',dims = PCUSE,n.components=nUMAP, check_duplicates=FALSE)
pbmc <- RunUMAP(object = pbmc, reduction='umap',dims = 1:nUMAP,n.components=2, check_duplicates=FALSE)

DimPlot(pbmc, reduction='umap', group.by='batch', pt.size=0.1) 


FeaturePlot(pbmc,reduction='umap',features = c('PCNA','NES'),ncol=2,cols=c('blue','gold','red'))






FeaturePlot(pbmc,reduction='umap',features = c('PCNA','NES'),ncol=2,cols=c('blue','gold','red'))

FeaturePlot(pbmc,reduction='umap',features = CCG,ncol=2,cols=c('blue','gold','red'))

FeaturePlot(pbmc,reduction='umap',features = 'ccg',ncol=2,cols=c('blue','gold','red'))



VEC = pbmc@reductions$umap@cell.embeddings
rownames(VEC) = colnames(pbmc)
PCA = pbmc@reductions$pca@cell.embeddings

# Define pixel
OUT=vector.buildGrid(VEC, N=50,SHOW=TRUE)

# Build network
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)

# Calculate Margin Score (MS)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)

# Get pixel's MS
OUT=vector.gridValue(OUT,SHOW=TRUE)

# Find summit
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)

# Infer vector
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=TRUE)



















###########
UMAP=pbmc@reductions$umap@cell.embeddings
library(destiny, quietly = TRUE)
dm <- DiffusionMap(UMAP)

##########

DimPlot(pbmc, reduction='umap', group.by='batch', pt.size=0.1) 





nUMAP=round(length(PCUSE))
pbmc <- RunUMAP(object = pbmc, reduction='pca',dims = PCUSE,n.components=nUMAP, check_duplicates=FALSE)
pbmc <- RunUMAP(object = pbmc, reduction='umap',dims = 1:nUMAP,n.components=2, check_duplicates=FALSE)



DimPlot(pbmc, reduction='umap', group.by='batch', pt.size=0.1) 


FeaturePlot(pbmc,reduction='umap',features = c('NES','PCNA','SOX11','LHX1','PTPRZ1','PAX6'),ncol=3,cols=c('blue','gold','red'))







library(gmodels)
D=as.matrix(pbmc@assays$RNA@scale.data)
D=D[which(rownames(D) %in% VariableFeatures(pbmc)),]
PCA.OUT=fast.prcomp(t(D), retx = TRUE, center = FALSE, scale. = FALSE, tol = NULL)
PCA=PCA.OUT$x




PCUSE <- mybeer$select




nUMAP=round(length(PCUSE)/2)
pbmc <- RunUMAP(object = pbmc, reduction='pca',dims = PCUSE,n.components=nUMAP, check_duplicates=FALSE)

pbmc <- RunUMAP(object = pbmc, reduction='umap',dims = 1:nUMAP,n.components=2, check_duplicates=FALSE)


DimPlot(pbmc, reduction='umap', group.by='batch', pt.size=0.1) 


FeaturePlot(pbmc,reduction='umap',features = c('NES','PCNA','SOX11','LHX1','PTPRZ1','PAX6'),ncol=3,cols=c('blue','gold','red'))



















pbmc <- RunUMAP(object = pbmc, reduction.use='tsne',dims = PCUSE, check_duplicates=FALSE)




#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)


PCUSE <- mybeer$select[1:10]
pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims = PCUSE, check_duplicates=FALSE)
DimPlot(pbmc, reduction='umap', group.by='batch', pt.size=0.1) 
FeaturePlot(pbmc,features = 'PCNA', reduction.use='umap')



PCUSE <- mybeer$select[1:50]
pbmc <- RunTSNE(object = pbmc, reduction.use='pca',dims = PCUSE, check_duplicates=FALSE)



TSNEPlot(pbmc, group.by='batch', pt.size=0.1) 
FeaturePlot(pbmc,reduction='tsne',features = c('PCNA','SOX11','NES'))


#VEC=pbmc@reductions$umap@cell.embeddings
VEC=pbmc@reductions$tsne@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings

# Define pixel
OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)

# Build network
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)

# Calculate Margin Score (MS)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)

# Get pixel's MS
OUT=vector.gridValue(OUT,SHOW=TRUE)

# Find summit
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)

# Infer vector
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=TRUE)


FeaturePlot(pbmc,reduction='tsne',features = c('NES','PCNA','SOX11','LHX1','PTPRZ1','PAX6'),ncol=3,cols=c('blue','gold','red'))

FeaturePlot(pbmc,reduction='tsne',features = c('NES','PCNA','SOX11'))















source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseNeuralCrest_GSE129114')


D1=read.table('GSE129114_E8.5_whole_embryo_Wnt1_counts.txt',sep=' ',header=T,row.names=1)
D2=read.table('GSE129114_E9.5_anterior_Sox10_counts.txt',sep=' ',header=T,row.names=1)
D3=read.table('GSE129114_E9.5_posterior_Sox10_counts.txt',sep=' ',header=T,row.names=1)
D4=read.table('GSE129114_E9.5_trunk_Wnt1_counts.txt',sep=' ',header=T,row.names=1)
D5=read.table('GSE129114_E10.5_post-otic_Sox10_counts.txt',sep=' ',header=T,row.names=1)
D6=read.table('GSE129114_E10.5_tail_Wnt1_counts.txt',sep=' ',header=T,row.names=1)
D7=read.table('GSE129114_E10.5_trunk_Wnt1_counts.txt',sep=' ',header=T,row.names=1)

BATCH=c(rep('E8.5.WH',ncol(D1)),
	rep('E9.5.AN.SOX10',ncol(D2)),
	rep('E9.5.PO.SOX10',ncol(D3)),
	rep('E9.5.TR.WNT1',ncol(D4)),
	rep('E10.5.PP.SOX10',ncol(D5)),
	rep('E10.5.TA.WNT1',ncol(D6)),
	rep('E10.5.TR.WNT1',ncol(D7))
	   )


       
D12=.simple_combine(D1,D2)$combine
D34=.simple_combine(D3,D4)$combine
D56=.simple_combine(D5,D6)$combine
D1234=.simple_combine(D12,D34)$combine
D123456=.simple_combine(D1234,D56)$combine
DATA=.simple_combine(D123456,D7)$combine

saveRDS(DATA,file='DATA.RDS')
saveRDS(BATCH,file='BATCH.RDS')

##################################



mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=150, ROUND=1, GN=5000, SEED=1, COMBAT=TRUE, RMG=NULL)   
saveRDS(mybeer,file='mybeer.RDS')




################
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseNeuralCrest_GSE129114')
mybeer=readRDS(file='mybeer.RDS')
#####################


##########
#E9.5.TR.WNT1
pbmc <- CreateSeuratObject(counts = D4, project = "pbmc3k", min.cells = 0, min.features = 0)
pbmc@meta.data$type=TYPE

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs = 150)
pbmc <- RunUMAP(pbmc, dims = 1:50)
DimPlot(pbmc, reduction = "umap")


FeaturePlot(pbmc,features = c('Sox9','Mafb','Olig3','Dlx5','Ret'))

saveRDS(pbmc,file='pbmc_D4.RDS')








VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings

OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)

































###############

#Biological Pathway

source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseOligo_GSE75330')
PCA.OUT=readRDS(file='PCA.OUT_500.RDS')
pbmc=readRDS('pbmc.RDS')

PCA=PCA.OUT$x

VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA=PCA #pbmc@reductions$pca@cell.embeddings



EXP=as.matrix(pbmc@assays$RNA@data)
UP=toupper(rownames(EXP))
USED=which(UP %in% names(which(table(UP)==1)))
EXP=EXP[USED,]
rownames(EXP)=UP[USED]




library(qusage)
BIO=read.gmt('MSigDB_Computational.gmt')

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

SUM=apply(BIO.MAT,1,sum)
BIO.MAT=BIO.MAT[which(SUM>0),]

TYPE=pbmc@meta.data$type


DrawHeatMap<-function(TAG){
    PCA=t(BIO.MAT)
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
	
    tiff(paste0("IMG/",TAG,".BIO.HEAT.tiff"),width=4,height=1,units='in',res=600)
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




#####################################################################











###############

#Biological Pathway

source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseIntestine_GSE92332/')

PCA.OUT=readRDS(file='PCA.OUT_500.RDS')
pbmc=readRDS('pbmc.RDS')

PCA=PCA.OUT$x

VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA=PCA #pbmc@reductions$pca@cell.embeddings

EXP=as.matrix(pbmc@assays$RNA@data)
EXP=t(apply(t(EXP),2,scale))
EXP[which(is.na(EXP))]=0
colnames(EXP)=colnames(pbmc)


UP=toupper(rownames(EXP))
USED=which(UP %in% names(which(table(UP)==1)))
EXP=EXP[USED,]
rownames(EXP)=UP[USED]




library(qusage)
BIO=read.gmt('MSigDB_Computational.gmt')

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

#BIO.MAT.ABS=vector.calValue(t(BIO.MAT))$PCA.RC




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
TYPE=pbmc@meta.data$type


DrawHeatMap<-function(TAG){
    PCA=t(BIO.MAT)#.ABS
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
	
    tiff(paste0("IMG/",TAG,".BIO.HEAT.tiff"),width=4,height=1,units='in',res=600)
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



#####################################################################















vector.gridValueSmooth <- function(OUT, SHOW=TRUE){
    OUT=OUT
    SHOW=SHOW
    ####################################
    #INDEX_LIST=OUT$INDEX_LIST
    #VALUE=OUT$VALUE
    USED=OUT$USED
    CENTER_VALUE=OUT$CENTER_VALUE
    CENTER_VEC=OUT$CENTER_VEC
    p1=OUT$p1
    p2=OUT$p2
    p1_index=as.numeric(str_replace(p1,'P',''))
    p2_index=as.numeric(str_replace(p2,'P',''))
    library(igraph)
    #################
    DEG=degree(OUT$GRAPH,v = V(OUT$GRAPH))
    names(DEG)=as_ids(V(OUT$GRAPH))
    #W=DEG/max(DEG)
    #W=rank(DEG)/length(DEG)
    DEG=DEG[order(as.numeric(str_replace(names(DEG),'P','')))]
    
    #################
    NEW_CENTER_VALUE=CENTER_VALUE
    NB_VALUE_MEAN=c()
    NB_VALUE_MIN=c()
    NB_VALUE_MAX=c()
    NB_VALUE_SD=c()
    NB_DEG_MEAN=c()
    ###########################################
    
    i=1
    while(i<=length(CENTER_VALUE)){
        this_value=CENTER_VALUE[i]
        neighbor_p1_index=p1_index[which(p2_index==i)]
        neighbor_p2_index=p2_index[which(p1_index==i)]
        neighbor_index=unique(sort(c(neighbor_p1_index,neighbor_p2_index)))
        nb_value=CENTER_VALUE[neighbor_index]
        nb_value_sd=sd(nb_value)
        nb_value_mean=mean(c(nb_value,this_value))
        nb_value_min=min(c(nb_value,this_value))
        nb_value_max=max(c(nb_value,this_value))
        #this_value=CENTER_VALUE[i]
        #this_value= nb_value_mean + (this_value-nb_value_mean)
        NB_VALUE_MEAN=c(NB_VALUE_MEAN,nb_value_mean)
        NB_VALUE_MIN=c(NB_VALUE_MIN,nb_value_min)
        NB_VALUE_MAX=c(NB_VALUE_MAX,nb_value_max)
        NB_VALUE_SD=c(NB_VALUE_SD,nb_value_sd) 
        NB_DEG_MEAN=c(NB_DEG_MEAN, mean(DEG[c(neighbor_index,i)]) )
        i=i+1}
    
    #RATIO = NB_VALUE_SD / max(NB_VALUE_SD)
    #W = NB_DEG_MEAN/max(NB_DEG_MEAN)
    #V =  NB_VALUE_MEAN / max(NB_VALUE_MEAN)
    
    RATIO=1#rank(NB_VALUE_SD)
    W=rank(NB_DEG_MEAN)
    V=rank(CENTER_VALUE)
    
    
    NEW_CENTER_VALUE =  V *  W  * RATIO# + NB_VALUE_MIN  #* W
    
    
    #####################################
    N.VALUE=.normX(NEW_CENTER_VALUE)#(VALUE-min(VALUE))/(max(VALUE)-min(VALUE))
    ORIG.CENTER.COL=vector.vcol(N.VALUE, c(0,0.5,1),c('#009FFF','#FFF200','#ec2F4B')) 
    #######################
    if(SHOW==TRUE){
            plot(OUT$VEC[,1],OUT$VEC[,2],col='grey80',cex=0.5, pch=16)
            points(CENTER_VEC[USED,1],CENTER_VEC[USED,2], col=ORIG.CENTER.COL[USED], pch=15, cex=1.5)
            }
    #################
    OUT$CENTER_VALUE=NEW_CENTER_VALUE       
    OUT$ORIG.CENTER.COL=ORIG.CENTER.COL
    return(OUT)
    }



vector.gridValueSmooth <- function(OUT,CUT=0.95, SHOW=TRUE, MAX=2000){
    OUT=OUT
    SHOW=SHOW
    CUT=CUT
    MAX=MAX
    INDEX_LIST=OUT$INDEX_LIST
    VALUE=OUT$VALUE
    USED=OUT$USED
    ################
    CENTER_VALUE=OUT$CENTER_VALUE
    CENTER_VEC=OUT$CENTER_VEC
    p1=OUT$p1
    p2=OUT$p2
    
    
    vector.smooth <- function(CUT, MAX, CENTER_VALUE, p1, p2){
    
        NEW_CENTER_VALUE=CENTER_VALUE
        p1_value=c()
        p2_value=c()
        DIFF=c()
        i=1
        while(i<=length(p1)){
            this_p1=p1[i]
            this_p2=p2[i]
            this_p1_index=as.numeric(str_replace(this_p1,'P',''))
            this_p2_index=as.numeric(str_replace(this_p2,'P',''))
            this_p1_value=NEW_CENTER_VALUE[this_p1_index]
            this_p2_value=NEW_CENTER_VALUE[this_p2_index]
            p1_value=c(p1_value, this_p1_value)
            p2_value=c(p2_value, this_p2_value)
            this_diff=this_p1_value-this_p2_value
            DIFF=c(DIFF, this_diff)
            i=i+1}    
        ABS_DIFF=abs(DIFF)
        POS_NUM=length(which(ABS_DIFF>0))
        ##############################################
        ABS_DIFF_COR=cor(1:length(ABS_DIFF),sort(ABS_DIFF))
        ############################
        COR_HIST=c()
        TIME=1
        while(ABS_DIFF_COR < CUT & TIME <= MAX){
            ################
            target_abs_diff= median(ABS_DIFF)
        
            ###############################
        
            this_max_index=which(ABS_DIFF==max(ABS_DIFF))[1]

            max_p1_index=as.numeric(str_replace( p1[this_max_index] ,'P',''))
            max_p2_index=as.numeric(str_replace( p2[this_max_index] ,'P',''))
            
            ###########################
            #NEW_CENTER_VALUE=NEW_CENTER_VALUE-target_abs_diff/2
            #NEW_CENTER_VALUE[which(NEW_CENTER_VALUE<0)]=0
            #################################
            if(NEW_CENTER_VALUE[max_p1_index] < NEW_CENTER_VALUE[max_p2_index]){
            
                NEW_CENTER_VALUE[max_p2_index]=NEW_CENTER_VALUE[max_p1_index] + target_abs_diff #/2          
                #NEW_CENTER_VALUE[max_p1_index] = max(0, NEW_CENTER_VALUE[max_p1_index]- target_abs_diff/2)
                ##########################################
                p1_changed_index=which( p1==p2[this_max_index] )
                ABS_DIFF[p1_changed_index]=abs( p1_value[p1_changed_index]- target_abs_diff - p2_value[p1_changed_index] )
                p2_changed_index=which( p2==p2[this_max_index] )
                ABS_DIFF[p2_changed_index]=abs( p1_value[p2_changed_index]+ target_abs_diff - p2_value[p2_changed_index] )
                ######################################
                }else{
            
                NEW_CENTER_VALUE[max_p1_index]=NEW_CENTER_VALUE[max_p2_index] + target_abs_diff #/2
                #NEW_CENTER_VALUE[max_p2_index]= max(0, NEW_CENTER_VALUE[max_p2_index]- target_abs_diff/2)
                ##########################################
                p1_changed_index=which( p1==p1[this_max_index] )
                ABS_DIFF[p1_changed_index]=abs( p1_value[p1_changed_index]+ target_abs_diff - p2_value[p1_changed_index] )
                
                p2_changed_index=which( p2==p1[this_max_index] )
                ABS_DIFF[p2_changed_index]=abs( p1_value[p2_changed_index]- target_abs_diff - p2_value[p2_changed_index] )
                ######################################
                }
        
            #ABS_DIFF[this_max_index]= target_abs_diff
            
            ##############################################
         
            ###############################
            ABS_DIFF_COR=cor(1:length(ABS_DIFF),sort(ABS_DIFF))
        
            ############################
            COR_HIST=c(COR_HIST, ABS_DIFF_COR)
            POS_NUM=length(which(ABS_DIFF>0))
            TIME=TIME+1
        }
        ###################################################
    
        NEW_CENTER_VALUE[which(!c(1:nrow(OUT$CENTER_VEC)) %in% USED)]=0
        NEW_CENTER_VALUE[USED]=.normX(NEW_CENTER_VALUE[USED])
        
        N.VALUE=.normX(NEW_CENTER_VALUE)#(VALUE-min(VALUE))/(max(VALUE)-min(VALUE))
        ORIG.CENTER.COL=vector.vcol(N.VALUE, c(0,0.5,1),c('#009FFF','#FFF200','#ec2F4B')) 
        #####################
        OUT=list()
        OUT$NEW_CENTER_VALUE=NEW_CENTER_VALUE
        OUT$TIME=TIME
        OUT$COR_HIST=COR_HIST
        OUT$ABS_DIFF_COR=ABS_DIFF_COR
        OUT$ORIG.CENTER.COL=ORIG.CENTER.COL
        return(OUT)
        }        
        
    SMOOTH.OUT= vector.smooth(CUT, MAX, CENTER_VALUE, p1, p2)
    ABS_DIFF_COR=SMOOTH.OUT$ABS_DIFF_COR
    NEW_CENTER_VALUE=SMOOTH.OUT$NEW_CENTER_VALUE
    ORIG.CENTER.COL=SMOOTH.OUT$ORIG.CENTER.COL
    COR_HIST=SMOOTH.OUT$COR_HIST
    TIME=SMOOTH.OUT$TIME
    ###############
        
    ####################################
    if(ABS_DIFF_COR<CUT){
        
        print('CUT is too high!!!')
        print(paste0('Max CUT should be less than: ', max(COR_HIST) ) )
        print(paste0('CUT is changed to: ',  max(COR_HIST) ) )
        ##########################
        SMOOTH.OUT= vector.smooth(max(COR_HIST), MAX, CENTER_VALUE, p1, p2)
        ABS_DIFF_COR=SMOOTH.OUT$ABS_DIFF_COR
        NEW_CENTER_VALUE=SMOOTH.OUT$NEW_CENTER_VALUE
        ORIG.CENTER.COL=SMOOTH.OUT$ORIG.CENTER.COL
        COR_HIST=SMOOTH.OUT$COR_HIST
        TIME=SMOOTH.OUT$TIME
        
        }
    
    if(SHOW==TRUE){
            plot(OUT$VEC[,1],OUT$VEC[,2],col='grey80',cex=0.5, pch=16)
            points(CENTER_VEC[USED,1],CENTER_VEC[USED,2], col=ORIG.CENTER.COL[USED], pch=15, cex=1.5)
            }   
    OUT$CENTER_VALUE=NEW_CENTER_VALUE
    OUT$ORIG.CENTER.COL=ORIG.CENTER.COL
    OUT$COR_HIST=COR_HIST
    return(OUT)
    }



vector.getValue_old <-function(PCA){
    PCA=PCA
    r.pca=apply(abs(PCA), 2, rank)
    mean.r.pca=apply(r.pca,1,mean)
    return(mean.r.pca)
    }


vector.autoCenterNew <- function(OUT,  SHOW=TRUE){
    #####################
    OUT=OUT
    SHOW=TRUE
    #############################
    DIST=OUT$DIST
    USED=OUT$USED
    USED_NAME=OUT$USED_NAME
    CENTER_VALUE=OUT$CENTER_VALUE
    CENTER_VEC=OUT$CENTER_VEC
    CENTER_INDEX=OUT$CENTER_INDEX
    INDEX_LIST=OUT$INDEX_LIST
    ############
    .tmp<-function(x){
        return(cor(x,CENTER_VALUE[USED],method='spearman'))
        }
    DIST_USED.COR=apply(DIST[USED,USED],2,.tmp)
    SUMMIT=USED[which(DIST_USED.COR==min(DIST_USED.COR))]
    SUMMIT_NAME=paste0('P',SUMMIT)
    SCORE=c()
    i=1
    while(i<=length(USED_NAME)){
        this_name=USED_NAME[i]
        this_dist=DIST[which(colnames(DIST)==this_name),which(rownames(DIST) %in% paste0('P',SUMMIT))]
        this_dist=this_dist
        #this_value=CENTER_VALUE[USED]
        this_score=min(this_dist)#sum(rank(-this_dist) * rank(this_value))
        #this_cor=cor(-this_value, this_dist)#,method='spearman')
        SCORE=c(SCORE, this_score)
        i=i+1}
    SCORE=max(SCORE)-SCORE
    VALUE=SCORE
    #plot(OUT$VEC, col='grey70',pch=16)
    N.VALUE=(VALUE-min(VALUE))/(max(VALUE)-min(VALUE))
    COL=vector.vcol(N.VALUE, c(0,0.5,1),c('#009FFF','#FFF200','#ec2F4B'))
    
    
    ###################################################
    OUT$COL=rep('grey70',nrow(OUT$VEC))
    OUT$ORIG.COL=rep('grey70',nrow(OUT$VEC))
    OUT$P.SCORE=rep(0,nrow(OUT$VEC))
    i=1
    while(i<=length(USED)){
        this_index=INDEX_LIST[[USED[i]]]
        OUT$COL[this_index]=COL[i]
        OUT$ORIG.COL[this_index]=OUT$ORIG.CENTER.COL[USED][i]
        OUT$P.SCORE[this_index]=SCORE[i]
        i=i+1
        }
    ###########
    ######################
    ############
    OUT$SCORE=SCORE
    OUT$SUMMIT=SUMMIT
    ################################
    if(SHOW==TRUE){
        #plot(OUT$VEC, col=OUT$COL, pch=16, cex=0.5 )
        plot(OUT$VEC, col=OUT$ORIG.COL, pch=16, cex=0.5)
        #text(CENTER_VEC[HIGH,],labels=PCH,cex=1,pos=2)
        #points(CENTER_VEC[HIGH,], col='black',pch=16,cex=1)
        points(CENTER_VEC[SUMMIT,1], CENTER_VEC[SUMMIT,2], col='black',pch=16,cex=1.5)
        points(CENTER_VEC[SUMMIT,1], CENTER_VEC[SUMMIT,2], col='red',pch=16,cex=1)  
        }
    #######################
    return(OUT)


     }





vector.calScore <-function(OUT,PCA,SHOW=TRUE){
    OUT=OUT
    INDEX_LIST=OUT$INDEX_LIST
    #PCA=t(DATA[which(rownames(DATA) %in% VariableFeatures(pbmc)),])
    PCA=PCA
    SHOW=SHOW
    DIST=OUT$DIST
    CENTER_VEC=OUT$CENTER_VEC
    #CENTER_PCA=c()
    USED=OUT$USED
    USED_NAME=OUT$USED_NAME
    ################
    CENTER_PCA=c()
    i=1
    while(i<=length(INDEX_LIST)){
        this_index=INDEX_LIST[[i]]
        if(length(this_index)==1){
            CENTER_PCA=cbind(CENTER_PCA, PCA[this_index,])
            }else{
            CENTER_PCA=cbind(CENTER_PCA, apply(PCA[this_index,],2,mean))
            }
        #print(i)
        i=i+1}
    ######################################
    CENTER_PCA=t(CENTER_PCA)
    rownames(CENTER_PCA)=colnames(DIST)
    colnames(CENTER_PCA)=colnames(PCA)
    ####################################

    used_p_index=which(OUT$p1 %in% USED_NAME & OUT$p2 %in% USED_NAME)
    used_p1=OUT$p1[used_p_index]
    used_p2=OUT$p2[used_p_index]
    
    CENTER_PCA.N=apply(CENTER_PCA,2,.normX)
    PCA.N=apply(PCA,2,rank)
    
    LENGTH=c()
    SSS=c()
    i=1
    while(i<=length(USED)){
        #this_pca_dist=#CENTER_PCA_SCORE[USED]-CENTER_PCA_SCORE[USED[i]]#CENTER_PCA_DIST[USED[i],USED]
        #CENTER_PCA[USED[i],]
        
        this_net_dist=DIST[USED[i],USED]
        
        this_net_dist_cluster=round(.normX(rank(this_net_dist)),1)
        step_list=sort(unique(this_net_dist_cluster))
        
        step_score=c()
        used_list=c()
        j=1
        while(j<=length(step_list)){
            this_step=step_list[j]
            this_index=which(this_net_dist_cluster==this_step)
            this_orig_index=c()
            t=1
            while(t<=length(this_index)){
                this_orig_index=c(this_orig_index, INDEX_LIST[[USED[this_index[t]]]])
                t=t+1}
               
            if(length(this_orig_index)>1){
                this_tmp=abs(cor(PCA[this_orig_index,],method='spearman'))
                this_step_score = mean(unique(this_tmp))
                step_score=c(step_score,this_step_score)
                used_list=c(used_list, this_step)
                }
            j=j+1}
        
        this_score=cor(used_list, step_score,method='spearman')#mean(this_mut_list)
        SSS=c(SSS,this_score)
        
        if(i %%10==1){print(paste0(i,'/',length(USED)))}
        i=i+1}
        
        
    #plot(SSS)
  
    #################################
    VALUE=SSS#TMP
    N.VALUE=(VALUE-min(VALUE))/(max(VALUE)-min(VALUE))
    COL=vector.vcol(N.VALUE, c(0,0.5,1),c('#009FFF','#FFF200','#ec2F4B'))
    
    
    SSS.N=.normX(SSS)
    
    SCORE=rep(0,nrow(OUT$VEC))
    i=1
    while(i<=length(USED)){
        this_index=INDEX_LIST[[USED[i]]]
        SCORE[this_index]=SSS.N[i]
        if(is.na(SSS.N[USED[i]])){print(i)}
        #SUMMIT_NEG_DIST[USED[i]]
        i=i+1
        }
    
    names(SCORE)=rownames(OUT$VEC)
    
    VALUE=SCORE
    N.VALUE=(VALUE-min(VALUE))/(max(VALUE)-min(VALUE))
    COL=vector.vcol(N.VALUE, c(0,0.5,1),c('#009FFF','#FFF200','#ec2F4B'))
    ################
    if(SHOW==TRUE){
        
        plot(OUT$VEC,col=COL,pch=16,cex=0.5)

        }
   
    return(SCORE)

    }




