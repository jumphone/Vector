<img src="https://github.com/jumphone/BEER/blob/master/DATA/Vector_LOGO.png" width="266">

# Developing Vector Inference

#### Environment: R (3.6.1)

#### Please install following R packages before using VECTOR:

    install.packages('circlize')
    install.packages('gatepoints')
    install.packages('stringr')
    install.packages('igraph')

## Usage:

### Please prepare a Seurat object with UMAP and 150 PCs.
Users can follow https://satijalab.org/seurat/ to generate Seurat object.

    library(Seurat)
    # DATA: Expression matrix. Rownames are gene names. Colnames are cell names.
    pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)
    pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

    pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
    all.genes <- rownames(pbmc)
    pbmc <- ScaleData(pbmc, features = all.genes)
    pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs = 150)
    pbmc <- RunUMAP(pbmc, dims = 1:50)
    DimPlot(pbmc, reduction = "umap")
    saveRDS(pbmc,file='pbmc.RDS')

### Get UMAP and PCs from Seurat3, pbmc: a Seurat object:

    VEC = pbmc@reductions$umap@cell.embeddings
    rownames(VEC) = colnames(pbmc)
    PCA = pbmc@reductions$pca@cell.embeddings

### (Optional) Get UMAP and PCs from Monocle3, cds: a Monocle object:
   
    VEC = cds@reducedDims$UMAP
    colnames(VEC) = c('UMAP_1','UMAP_2')
    PCA = cds@reducedDims$PCA

###  Use VECTOR:


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




#### Change MS to NES gene expression:

    NES.EXP = pbmc@assays$RNA@data[which(rownames(pbmc) =='Nes'),]
    OUT=vector.buildGrid(VEC, N=40,SHOW=TRUE)
    OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
    OUT=vector.getValue(OUT, PCA, SHOW=TRUE)

    OUT$VALUE=NES.EXP

    OUT=vector.showValue(OUT)
    OUT=vector.gridValue(OUT, SHOW=TRUE)
    OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
    OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)

    
#### Select starting point:

    OUT=vector.buildGrid(VEC, N=40,SHOW=TRUE)
    OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
    OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
    OUT=vector.gridValue(OUT,SHOW=TRUE)

    OUT=vector.selectCenter(OUT)

    OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)

#### Select region of interest:

    OUT=vector.buildGrid(VEC, N=40,SHOW=TRUE)
    OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
    OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
    OUT=vector.gridValue(OUT,SHOW=TRUE)
    OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
    OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)

    OUT=vector.reDrawArrow(OUT, COL=OUT$COL)
    OUT=vector.selectRegion(OUT)

    SELECT_PS=OUT$SELECT_PS
    SELECT_INDEX=OUT$SELECT_INDEX


