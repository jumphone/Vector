#Citation: 
# https://data.humancellatlas.org/explore/projects/cc95ff89-2e68-4a08-a234-480eca21ce79

#Python:

fi=open('F:/Vector/data/HumanBoneMarrow/5635f5a8-532e-483b-81a2-9c53495dbd90.csv/expression.csv')
fo=open('F:/Vector/data/HumanBoneMarrow/5635f5a8-532e-483b-81a2-9c53495dbd90.csv/expression_10000.csv','w')
l1=fi.readline()
fo.write(l1)
i=1
while i<=10000:
    fo.write(fi.readline())
    print(i)
    i=i+1
fi.close()
fo.close()
#################################

#R:
setwd('F:/Vector/data/HumanBoneMarrow/5635f5a8-532e-483b-81a2-9c53495dbd90.csv/')
DATA=read.csv(file='expression_10000.csv',sep=',',header=TRUE,row.names=1)
DATA=t(DATA)
saveRDS(DATA,file='DATA.RDS')



library(Seurat)

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

###############################

###########################
setwd('F:/Vector/data/MouseBoneMarrow_GSE109989/')
pbmc=readRDS(file='pbmc.RDS')
# https://www.rndsystems.com/cn/pathways/hematopoietic-stem-cell-differentiation-pathways-lineage-specific-markers

# HSC
# GATA2:ENSG00000179348
# AK2: ENSG00000004455


FeaturePlot(pbmc, features=c('ENSG00000179348',
                            'ENSG00000102145'
                            ),
           sort.cell=TRUE, pt.size=1)


VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings

OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)







