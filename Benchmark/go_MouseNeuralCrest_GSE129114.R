


source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseNeuralCrest_GSE129114')


D4=read.table('GSE129114_E9.5_trunk_Wnt1_counts.txt',sep=' ',header=T,row.names=1)


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


FeaturePlot(pbmc,features = c('Mafb','Olig3','Dlx5','Ret'))

saveRDS(pbmc,file='pbmc_D4.RDS')





