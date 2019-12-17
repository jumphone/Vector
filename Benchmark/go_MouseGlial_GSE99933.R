source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

setwd('F:/Vector/data/MouseGlial_GSE99933')



D1=read.csv('E12.5_counts.txt',header=T,row.names=1,sep='\t')
D2=read.csv('E13.5_counts.txt',header=T,row.names=1,sep='\t')

DATA=.simple_combine(D1,D2)$combine
TYPE=c(rep('E12.5',ncol(D1)),rep('E13.5',ncol(D2)))


pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)


pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs = 150)
#pbmc <- RunUMAP(pbmc, dims = 1:20, n.neighbors=10, min.dist=0.5,spread=1)
pbmc <- RunUMAP(pbmc, dims = 1:50)

pbmc@meta.data$type=TYPE
DimPlot(pbmc, reduction = "umap",group.by='type',label=TRUE)
saveRDS(pbmc,file='pbmc.RDS')

##################

source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseGlial_GSE99933')

pbmc=readRDS(file='pbmc.RDS')


VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings




OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL,AL=40)

tiff(paste0("IMG/NEW_VECTOR.6.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL,AL=40)
dev.off()







OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)

OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)

OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)


FeaturePlot(pbmc,features=c('Sox10','Htr3a','Th','Cartpt'))





tiff(paste0("IMG/VECTOR.1.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
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
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL,AL=90)
dev.off()




tiff(paste0("IMG/Marker.tiff"),width=5,height=5,units='in',res=600)
p <- FeaturePlot(pbmc, features=c('Sox10','Htr3a','Th','Cartpt'),order=TRUE, combine = FALSE)

for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}

cowplot::plot_grid(plotlist = p)
dev.off()


















vector.drawArrow <- function(OUT, P=0.9, SHOW=TRUE, COL='grey70',AL=70){
    ################
    OUT=OUT
    SHOW=SHOW
    P=P
    COL=COL
    USED=OUT$USED
    DIST=OUT$DIST
    ALL_VEC=OUT$VEC
    USED_NAME=OUT$USED_NAME
    USED_CENTER_VEC=OUT$CENTER_VEC[USED,]
    USED_DIST=OUT$DIST[which(rownames(DIST) %in% USED_NAME),which(rownames(DIST) %in% USED_NAME)]
    SCORE=OUT$SCORE
    AL=AL

    ###################
    .norm_one <-function(x,one=1){
        one=one
        if(var(x)!=0){
            x=x/sqrt(sum(x^2)) * one }
        return(x)
        }
    
    DIV=1/P
    if(SHOW==TRUE){
        plot(ALL_VEC,col=COL,pch=16,cex=0.2)
        }
    
    N.SCORE=.normX(SCORE)
    SCORE.COL=vector.vcol(N.SCORE, c(0,0.5,1),c('#009FFF','#FFF200','#ec2F4B'))
    
    A1_VEC=c()
    A2_VEC=c()
    A_LENGTH=c()
    one=1
    #one=min(dist(USED_CENTER_VEC))
    i=1
    while(i<=length(USED)){
        this_p1_loc=USED_CENTER_VEC[i,]
        
        vector_list=cbind(USED_CENTER_VEC[,1]-this_p1_loc[1],USED_CENTER_VEC[,2]-this_p1_loc[2])
        vector_list_norm=t(apply(vector_list,1,.norm_one, one))
        
        vector_weight_1= DIV^-(rank(USED_DIST[i,])-1)  
        vector_weight_2= SCORE[i]-SCORE
        
        vector_weight = vector_weight_1 * vector_weight_2        
        vector_weight = vector_weight/sum(abs(vector_weight))
        
        final_vec=t(vector_list_norm) %*% vector_weight
        
        this_p2_loc=c(this_p1_loc[1]+final_vec[1],this_p1_loc[2]+final_vec[2])
                
        #plot(ALL_VEC,col=COL,pch=16,cex=0.2)
        this_arrow_length=0.1*sqrt(sum(final_vec^2))
        #this_arrow_length=sqrt(sum(final_vec^2))#0.1 #* (1+sqrt(sum(final_vec^2)))
        #this_arrow_length=dev.size()[1]/AL *  sqrt(sum(final_vec^2))/one   # * sqrt(sum(final_vec^2)) #0.25
        if(SHOW==TRUE){
            arrows(x0=this_p1_loc[1],y0=this_p1_loc[2],
                   x1=this_p2_loc[1],y1=this_p2_loc[2],
                   lwd=2, length=this_arrow_length,
                   col='black'
                   #col=SCORE.COL[i]
                   )
            }
        A1_VEC=cbind(A1_VEC,this_p1_loc)
        A2_VEC=cbind(A2_VEC,this_p2_loc)
        A_LENGTH=c(A_LENGTH,this_arrow_length)
        i=i+1}
    #################################
    A1_VEC=t(A1_VEC)
    A2_VEC=t(A2_VEC)
   ######
    #OUT=list()
    OUT$A1_VEC=A1_VEC
    OUT$A2_VEC=A2_VEC
    OUT$A_LENGTH=A_LENGTH
    OUT$A_COL=SCORE.COL
    ###########
    return(OUT)
    }




##########################




###########################
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseGlial_GSE99933')

pbmc=readRDS(file='pbmc.RDS')

VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings


OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
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
#            msPCA      nGene     msGene     varPCA
#msPCA   1.0000000 -0.4092968  0.5214183  0.8720072
#nGene  -0.4092968  1.0000000 -0.9510347 -0.3504877
#msGene  0.5214183 -0.9510347  1.0000000  0.4741413
#varPCA  0.8720072 -0.3504877  0.4741413  1.0000000




######################
# Test pathway feature



###########################
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseGlial_GSE99933')

pbmc=readRDS(file='pbmc.RDS')

VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings


OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
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


###########################
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseGlial_GSE99933')

pbmc=readRDS(file='pbmc.RDS')

VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings


OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
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




set.seed(123)

pdf('./IMG/TEST_random.pdf')
OUT$VALUE=rnorm(length(OUT$VALUE))
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=TRUE)
dev.off()
VALUE$random=OUT$VALUE

saveRDS(MAT,'IMG/TEST_RANDOM_MAT.RDS')






