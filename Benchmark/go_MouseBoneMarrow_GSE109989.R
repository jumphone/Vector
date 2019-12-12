source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

########################################
setwd('F:/Vector/data/MouseBoneMarrow_GSE109989/')

DATA=read.csv(file='count_matrix.csv',sep=',',header=TRUE,row.names=1)
saveRDS(DATA,file='DATA.RDS')



########################################
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseBoneMarrow_GSE109989/')

DATA=readRDS('DATA.RDS')

pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)
#pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs = 150)
pbmc <- RunUMAP(pbmc, dims = 1:150)
DimPlot(pbmc, reduction = "umap")
saveRDS(pbmc,file='pbmc.RDS')




###########################
setwd('F:/Vector/data/MouseBoneMarrow_GSE109989/')
pbmc=readRDS(file='pbmc.RDS')

FeaturePlot(pbmc, features=c('Gfi1','Cenpa','Cd79a','Mpeg1','Cd47'))

#######################
VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings




OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
tiff(paste0("IMG/NEW_VECTOR.6.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL,AL=40)
dev.off()






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
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)
dev.off()



tiff(paste0("IMG/Marker.tiff"),width=5.5,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
FeaturePlot(pbmc, features=c('Cenpa','Il1b','Cd79a','Mpeg1'),order=TRUE)
dev.off()




tiff(paste0("IMG/Marker.tiff"),width=5,height=5,units='in',res=600)
p <- FeaturePlot(pbmc, features=c('Cenpa','Il1b','Cd79a','Mpeg1'),order=TRUE, combine = FALSE)

for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}

cowplot::plot_grid(plotlist = p)
dev.off()








###########################
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseBoneMarrow_GSE109989/')
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


pdf('./IMG/TEST_VAR_rPCA.pdf')
rPCA=apply(PCA,2,rank)
OUT$VALUE=apply(rPCA,1,var)
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=TRUE)
dev.off()
VALUE$varRPCA=OUT$VALUE

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
#msPCA      nGene     msGene     varPCA
#msPCA   1.0000000  0.6178928 -0.6173950  0.9674628
#nGene   0.6178928  1.0000000 -0.9819733  0.6700861
#msGene -0.6173950 -0.9819733  1.0000000 -0.6523595
#varPCA  0.9674628  0.6700861 -0.6523595  1.0000000




