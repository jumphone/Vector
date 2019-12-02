source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

########################################
setwd('F:/AIM/MouseBoneMarrow_GSE109989/')

DATA=read.table(file='count_matrix.csv',sep='\t',header=TRUE,row.names=1)
saveRDS(DATA,file='DATA.RDS')




pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)
#pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs = 50)
pbmc <- RunUMAP(pbmc, dims = 1:50)
DimPlot(pbmc, reduction = "umap")

