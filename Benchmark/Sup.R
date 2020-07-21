#########################################################
## Olig

setwd('F:/Vector/data/MouseOligo_GSE75330')
library(Seurat)

pbmc=readRDS('pbmc.RDS')

VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings

R.PCA=apply(PCA,2,rank)
RR.PCA=R.PCA/nrow(R.PCA)

CUT=10

MAT=matrix(0,ncol=ncol(PCA),nrow=nrow(PCA))
rownames(MAT)=rownames(PCA)
colnames(MAT)=colnames(PCA)
MAT[which(R.PCA>nrow(PCA)-CUT |  R.PCA <= CUT)]=1


#########################################
TAG=as.character(pbmc@meta.data$type)
TAG[which(TAG=='OPC')]='OPC'
TAG[which(TAG=='Differentiation-committed OPC')]='COP'
TAG[which(TAG=='Newly-formed Oligodendrocytes')]='NFOL'
TAG[which(TAG=='Myelin-forming Oligodendrocytes')]='MFOL'
TAG[which(TAG=='Mature Oligodendrocytes')]='MOL'
################################################

TMP=c()
this_mat=MAT[which(TAG=='OPC'),]
ratio_out=apply(this_mat, 2, sum)/nrow(this_mat)
TMP=c(TMP,length(which(ratio_out>0))/ncol(PCA))
#plot(ratio_out,ylim=c(0,1),type='l')

this_mat=MAT[which(TAG=='COP'),]
ratio_out=apply(this_mat, 2, sum)/nrow(this_mat)
TMP=c(TMP,length(which(ratio_out>0))/ncol(PCA))
#lines(ratio_out,type = "l",col='red')

this_mat=MAT[which(TAG=='NFOL'),]
ratio_out=apply(this_mat, 2, sum)/nrow(this_mat)
TMP=c(TMP,length(which(ratio_out>0))/ncol(PCA))
#lines(ratio_out,type = "l",col='green')

this_mat=MAT[which(TAG=='MFOL'),]
ratio_out=apply(this_mat, 2, sum)/nrow(this_mat)
TMP=c(TMP,length(which(ratio_out>0))/ncol(PCA))
#lines(ratio_out,type = "l",col='blue')

this_mat=MAT[which(TAG=='MOL'),]
ratio_out=apply(this_mat, 2, sum)/nrow(this_mat)
TMP=c(TMP,length(which(ratio_out>0))/ncol(PCA))
#lines(ratio_out,type = "l",col='gold')

names(TMP)=c('OPC','COP','NFOL','MFOL','MOL')
barplot(TMP,ylim=c(0,1))










N.TMP=TMP
N.TMP[1]=length(which(TAG=='OPC'))
N.TMP[2]=length(which(TAG=='COP'))
N.TMP[3]=length(which(TAG=='NFOL'))
N.TMP[4]=length(which(TAG=='MFOL'))
N.TMP[5]=length(which(TAG=='MOL'))
barplot(N.TMP)











##########################################################

# Intestine

setwd('F:/Vector/data/MouseIntestine_GSE92332/')
library(Seurat)

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


TAG=pbmc@meta.data$type


R.PCA=apply(PCA,2,rank)
RR.PCA=R.PCA/nrow(R.PCA)

CUT=10

MAT=matrix(0,ncol=ncol(PCA),nrow=nrow(PCA))
rownames(MAT)=rownames(PCA)
colnames(MAT)=colnames(PCA)
MAT[which(R.PCA> nrow(PCA)-CUT |  R.PCA <= CUT)]=1



TMP=c()
this_mat=MAT[which(TAG=='STEM'),]
ratio_out=apply(this_mat, 2, sum)/nrow(this_mat)
TMP=c(TMP,length(which(ratio_out>0))/ncol(PCA))
#plot(ratio_out,ylim=c(0,1),type='l')

this_mat=MAT[which(TAG=='TA'),]
ratio_out=apply(this_mat, 2, sum)/nrow(this_mat)
TMP=c(TMP,length(which(ratio_out>0))/ncol(PCA))
#lines(ratio_out,type = "l",col='red')

this_mat=MAT[which(TAG=='EP'),]
ratio_out=apply(this_mat, 2, sum)/nrow(this_mat)
TMP=c(TMP,length(which(ratio_out>0))/ncol(PCA))
#lines(ratio_out,type = "l",col='green')

this_mat=MAT[which(TAG=='EIM'),]
ratio_out=apply(this_mat, 2, sum)/nrow(this_mat)
TMP=c(TMP,length(which(ratio_out>0))/ncol(PCA))
#lines(ratio_out,type = "l",col='blue')

this_mat=MAT[which(TAG=='EM'),]
ratio_out=apply(this_mat, 2, sum)/nrow(this_mat)
TMP=c(TMP,length(which(ratio_out>0))/ncol(PCA))
#lines(ratio_out,type = "l",col='gold')

names(TMP)=c('STEM','TA','EP','EIM','EM')
barplot(TMP,ylim=c(0,1))












N.TMP=TMP
N.TMP[1]=length(which(TAG=='STEM'))
N.TMP[2]=length(which(TAG=='TA'))
N.TMP[3]=length(which(TAG=='EP'))
N.TMP[4]=length(which(TAG=='EIM'))
N.TMP[5]=length(which(TAG=='EM'))
barplot(N.TMP)
barplot(TMP/N.TMP)









