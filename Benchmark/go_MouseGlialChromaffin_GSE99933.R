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




tiff(paste0("IMG/NEW_VECTOR_TRY.6.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL,AL=60,BD=FALSE,AW=1,AC='black')
dev.off()



tiff(paste0("IMG/NEW_VECTOR.6.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL,AL=40)
dev.off()


tiff(paste0("IMG/Marker.tiff"),width=5,height=5,units='in',res=600)
p <- FeaturePlot(pbmc, features=c('Sox10','Htr3a','Th','Cartpt'),order=TRUE, pt.size=0.3,combine = FALSE)

for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}

cowplot::plot_grid(plotlist = p)
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







