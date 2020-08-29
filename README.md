<img src="https://github.com/jumphone/BEER/raw/master/DATA/Vector_LOGO.png" width="266">

# Single-cell Developing Vector Inference

#### Environment: R (3.6.1)

#### Please install following R packages before using VECTOR:

    install.packages('circlize')
    install.packages('gatepoints')
    install.packages('stringr')
    install.packages('igraph')
    install.packages('gmodels')

## Usage:

### Step 1. Please prepare a Seurat object with UMAP and 150 PCs.
Users can follow https://satijalab.org/seurat/ to generate Seurat object (V3.0.0).

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

### Step 2. Get UMAP and PCs from Seurat3 object. (pbmc: a Seurat object):

    VEC = pbmc@reductions$umap@cell.embeddings
    rownames(VEC) = colnames(pbmc)
    PCA = pbmc@reductions$pca@cell.embeddings
    
    source('https://raw.githubusercontent.com/jumphone/Vector/master/Vector.R')
    
    # Remove quantile-based colinearity among PCs (new feature in VECTOR 0.0.3):   
    PCA=vector.rankPCA(PCA)

### Step 3. Use VECTOR:
<img src="https://raw.githubusercontent.com/jumphone/BEER/master/DATA/TMP/WF.jpg" width="600">

    source('https://raw.githubusercontent.com/jumphone/Vector/master/Vector.R')

    # Define pixel
    OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
    
    # Build network
    OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
    
    # Calculate Quantile Polarization (QP) score
    OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
    
    # Get pixel's QP score
    OUT=vector.gridValue(OUT,SHOW=TRUE)
    
    # Find starting point
    OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
    
    # Infer vector
    OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=TRUE)


## Additional function 1: Change QP score to a given gene's expression value (e.g. Nes):
<img src="https://raw.githubusercontent.com/jumphone/BEER/master/DATA/TMP/WF1.jpg" width="400">

    NES.EXP = pbmc@assays$RNA@data[which(rownames(pbmc) =='Nes'),]
    OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
    OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
    OUT=vector.getValue(OUT, PCA, SHOW=TRUE)

    OUT$VALUE=NES.EXP

    OUT=vector.showValue(OUT)
    OUT=vector.gridValue(OUT, SHOW=TRUE)
    OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
    OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)

    
## Additional function 2: Manually select starting point:
<img src="https://raw.githubusercontent.com/jumphone/BEER/master/DATA/TMP/WF2.jpg" width="200">

    OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
    OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
    OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
    OUT=vector.gridValue(OUT,SHOW=TRUE)

    OUT=vector.selectCenter(OUT)

    OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)

## Additional function 3: Manually select region of interest:
<img src="https://raw.githubusercontent.com/jumphone/BEER/master/DATA/TMP/WF3.jpg" width="200">

    OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
    OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
    OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
    OUT=vector.gridValue(OUT,SHOW=TRUE)
    OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
    OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)

    #######################
    OUT=vector.reDrawArrow(OUT, COL=OUT$COL)
    OUT=vector.selectRegion(OUT)
    
    #######################
    SELECT_PS=OUT$SELECT_PS               #Peseudotime Score (PS) of selected cells
    SELECT_INDEX=OUT$SELECT_INDEX         #Index of selected cells in the expression matrix 
    SELECT_COL=OUT$COL[OUT$SELECT_INDEX]  #Colors
   
    #######################
    # Identify development related genes
    EXP=as.matrix(pbmc@assays$RNA@data)[which(rownames(pbmc) %in% VariableFeatures(pbmc)),SELECT_INDEX]
    COR=c()
    i=1
    while(i<=nrow(EXP)){
        this_cor=cor(SELECT_PS, EXP[i,],method='spearman')
        COR=c(COR,this_cor)
        if(i %%100==1){print(i)}
        i=i+1}
    names(COR)=rownames(EXP)
    head(sort(COR),n=10)     #Decreasing (top 10)
    tail(sort(COR),n=10)     #Increasing (top 10) 
    
    # Select one gene to draw figure
    show_gene=names(head(sort(COR),n=10))[1]
    show_gene.exp=EXP[which(rownames(EXP)==show_gene),]
    
    # Smooth expression value along pesudotime order (optional)
    show_gene.exp[order(SELECT_PS)]=smooth.spline(show_gene.exp[order(SELECT_PS)], df=5)$y    
    
    # Draw figure
    plot(jitter(SELECT_PS), show_gene.exp, pch=16,col=SELECT_COL, ylab=show_gene,xlab='PS')
    show_gene.fit=lm(show_gene.exp~SELECT_PS)
    abline(show_gene.fit,col='black',lwd=1)
    
    

## Other: Get UMAP and PCs from Monocle3. (cds: a Monocle object):
   
    # Get UMAP:
    VEC = cds@reducedDims$UMAP
    colnames(VEC) = c('UMAP_1','UMAP_2')
    
    # Get 150 PCs
    library(Seurat)
    DATA=as.matrix(cds@assays$data[[1]])
    pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)
    pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
    pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
    all.genes <- rownames(pbmc)
    pbmc <- ScaleData(pbmc, features = all.genes)
    pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs = 150)
    PCA = pbmc@reductions$pca@cell.embeddings


