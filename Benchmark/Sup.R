# Olig

setwd('F:/Vector/data/MouseOligo_GSE75330')
library(Seurat)

pbmc=readRDS('pbmc.RDS')

VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings

R.PCA=apply(PCA,2,rank)
RR.PCA=R.PCA/nrow(R.PCA)


CUT=0.05

MAT=matrix(0,ncol=ncol(PCA),nrow=nrow(PCA))
rownames(MAT)=rownames(PCA)
colnames(MAT)=colnames(PCA)
MAT[which(RR.PCA>1-CUT |  RR.PCA <CUT)]=1


#########################################
TAG=as.character(pbmc@meta.data$type)
TAG[which(TAG=='OPC')]='OPC'
TAG[which(TAG=='Differentiation-committed OPC')]='COP'
TAG[which(TAG=='Newly-formed Oligodendrocytes')]='NFOL'
TAG[which(TAG=='Myelin-forming Oligodendrocytes')]='MFOL'
TAG[which(TAG=='Mature Oligodendrocytes')]='MOL'
################################################



