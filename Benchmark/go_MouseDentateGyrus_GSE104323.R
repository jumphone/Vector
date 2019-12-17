source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

########################################
setwd('F:/Vector/data/MouseDentateGyrus_GSE104323/')

DATA=read.table(file='10X_expression_data.tab',sep='\t',header=TRUE,row.names=1)
saveRDS(DATA,file='DATA.RDS')


#################################
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseDentateGyrus_GSE104323/')


DATA=readRDS('DATA.RDS')
pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)
#pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
BATCH=pbmc@meta.data$orig.ident
saveRDS(BATCH,file='BATCH.RDS')
###########################



mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=150, ROUND=1, GN=5000, SEED=1, COMBAT=TRUE, RMG=NULL)   
saveRDS(mybeer,file='mybeer.RDS')



########################
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseDentateGyrus_GSE104323/')


mybeer=readRDS(file='mybeer.RDS')
PCUSE=mybeer$select
COL=rep('black',length(mybeer$cor))
COL[PCUSE]='red'
plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,
    xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))


pbmc <- mybeer$seurat
PCUSE=mybeer$select   
pbmc=BEER.combat(pbmc) #Adjust PCs using ComBat
umap=BEER.bbknn(pbmc, PCUSE, NB=3, NT=10)
pbmc@reductions$umap@cell.embeddings=umap
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1,label=F)


############################

#colnames(pbmc)[1:3]
library(stringr)

META=read.table(file='GSE104323_metadata_barcodes_24185cells.txt',sep='\t',header=TRUE,row.names=NULL)
META=META[1:ncol(pbmc),]

META.NAME=str_replace(META[,1],'-','.')
META.NAME=str_replace(META.NAME,'10X','X10X')
META=META[match(colnames(pbmc), META.NAME),]
pbmc@meta.data$type=META$characteristics..cell.cluster

saveRDS(pbmc,file='pbmc.RDS')

#########################







##########################################
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseDentateGyrus_GSE104323/')
pbmc=readRDS(file='pbmc.RDS')




tiff(paste0("IMG/ALL.tiff"),width=5.5,height=5.5,units='in',res=600)
DimPlot(pbmc, group.by='type',label = TRUE)+NoLegend()+NoAxes()
dev.off()



VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings





OUT=vector.buildGrid(VEC, N=40,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL,AL=40)


tiff(paste0("IMG/NEW_VECTOR.6.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL,OL=2,AL=50)
dev.off()





tiff(paste0("IMG/VECTOR.1.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.buildGrid(VEC, N=40,SHOW=TRUE)
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










####
#Change MS to NES
NES.EXP = pbmc@assays$RNA@data[which(rownames(pbmc) =='Nes'),]


OUT=vector.buildGrid(VEC, N=40,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)



OUT$VALUE=NES.EXP

tiff(paste0("IMG/NES.1.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.showValue(OUT)
dev.off()

tiff(paste0("IMG/NES.2.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.gridValue(OUT, SHOW=TRUE)
dev.off()

tiff(paste0("IMG/NES.3.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.autoCenter(OUT,UP=0.95,SHOW=TRUE)
dev.off()

tiff(paste0("IMG/NES.4.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL,OL=2,AL=50)
dev.off()



###########
#Select starting point



OUT=vector.buildGrid(VEC, N=40,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
OUT=vector.gridValue(OUT,SHOW=TRUE)


OUT=vector.selectCenter(OUT)

tiff(paste0("IMG/MSP.2.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)
dev.off()





###########
#Select arrows


OUT=vector.buildGrid(VEC, N=40,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)

OUT=vector.reDrawArrow(OUT, COL=OUT$COL)
OUT=vector.selectRegion(OUT)


saveRDS(OUT,'Select.OUT.RDS')


SELECT_SCORE=OUT$SELECT_SCORE
SELECT_INDEX=OUT$SELECT_INDEX

EXP=as.matrix(pbmc@assays$RNA@scale.data)[,SELECT_INDEX]


print(nrow(EXP))
#10646
COR=c()
i=1
while(i<=nrow(EXP)){
    this_cor=cor(SELECT_SCORE, EXP[i,],method='spearman')
    COR=c(COR,this_cor)
    if(i %%100==1){print(i)}
    i=i+1}


names(COR)=rownames(EXP)

head(sort(COR),n=10)
#Tubb4a      Aplp1      Stmn4      Cryab       Plp1     Slain1      Cntn2        App 
#-0.8320211 -0.7922493 -0.7217838 -0.6819902 -0.6730290 -0.6663356 -0.6648373 -0.6537134 
#      Ank3      Pacs2 
#-0.6534026 -0.6464254


tail(sort(COR),n=10)
#Rpl22l1     Hmgb3       Ran     Fabp7      Sox9     Hmgb2     H2afv    Slc1a3     Hmgn1 
#0.5718373 0.5811396 0.5818859 0.5865809 0.5908710 0.6133952 0.6231035 0.6544590 0.6549967 
#    H2afz 
#0.7163760 


ORDER=order(-SELECT_SCORE)


SOX9=EXP[which(rownames(EXP)=='Sox9'),ORDER]
INDEX=1:length(SOX9)
SOX9.USED=which(SOX9>0)
SOX9.INDEX=INDEX[SOX9.USED]
SOX9=SOX9[SOX9.USED]
SOX9.fit=lm(SOX9~SOX9.INDEX)



TUBB=EXP[which(rownames(EXP)=='Tubb4a'),ORDER]
TUBB.INDEX=1:length(TUBB)
TUBB.USED=which(TUBB>0)
TUBB.INDEX=INDEX[TUBB.USED]
TUBB=TUBB[TUBB.USED]
TUBB.fit=lm(TUBB~TUBB.INDEX)





tiff(paste0("IMG/SOX9_TUBB.tiff"),width=2,height=2,units='in',res=600)

par(mar=c(0,0,0,0))
par(mfrow=c(2,1))
plot(SOX9.INDEX, SOX9, pch=16,col=OUT$COL[OUT$SELECT_INDEX[ORDER][SOX9.INDEX]],cex=0.5,axes=FALSE,frame.plot=FALSE)
abline(SOX9.fit,col='black',lwd=2)
plot(TUBB.INDEX, TUBB, pch=16,col=OUT$COL[OUT$SELECT_INDEX[ORDER][TUBB.INDEX]],cex=0.5,axes=FALSE,frame.plot=FALSE)
abline(TUBB.fit,col='black',lwd=2)

dev.off()








##################################3



###########################

source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseDentateGyrus_GSE104323/')
pbmc=readRDS(file='pbmc.RDS')


VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings


OUT=vector.buildGrid(VEC, N=40,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)


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
 #       msPCA       nGene       msGene      varPCA
#msPCA   1.0000000 -0.03457670  0.122639236 0.855785020
#nGene  -0.0345767  1.00000000 -0.988302236 0.085748360
#msGene  0.1226392 -0.98830224  1.000000000 0.007726971
#varPCA  0.8557850  0.08574836  0.007726971 1.000000000







######################
# Test pathway feature


###########################

source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseDentateGyrus_GSE104323/')
pbmc=readRDS(file='pbmc.RDS')


VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings


OUT=vector.buildGrid(VEC, N=40,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)


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

EXP=as.matrix(pbmc@assays$RNA@data)
EXP=t(apply(t(EXP),2,scale))
EXP[which(is.na(EXP))]=0
colnames(EXP)=colnames(pbmc)
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








