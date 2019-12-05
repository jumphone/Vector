source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

########################################
setwd('F:/Vector/data/MouseDentateGyrus_GSE104323/')

DATA=read.table(file='10X_expression_data.tab',sep='\t',header=TRUE,row.names=1)
saveRDS(DATA,file='DATA.RDS')


#################################
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseDentateGyrus_GSE104323/')


DATA=readRDS('DATA.RDS')
pbmc <- CreateSeuratObject(counts = DATA, project = "pbmc3k", min.cells = 0, min.features = 0)
#pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
BATCH=pbmc@meta.data$orig.ident
saveRDS(BATCH,file='BATCH.RDS')
###########################



mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=150, ROUND=1, GN=5000, SEED=1, COMBAT=TRUE, RMG=NULL)   
saveRDS(mybeer,file='mybeer.RDS')



########################
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseDentateGyrus_GSE104323/')


mybeer=readRDS(file='mybeer.RDS')
PCUSE=mybeer$select
COL=rep('black',length(mybeer$cor))
COL[PCUSE]='red'
plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,
    xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))


pbmc <- mybeer$seurat
PCUSE=mybeer$select   
pbmc=BEER.combat(pbmc) #Adjust PCs using ComBat
umap=BEER.bbknn(pbmc, PCUSE, NB=3, NT=10)
pbmc@reductions$umap@cell.embeddings=umap
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1,label=F)


############################

#colnames(pbmc)[1:3]
library(stringr)

META=read.table(file='GSE104323_metadata_barcodes_24185cells.txt',sep='\t',header=TRUE,row.names=NULL)
META=META[1:ncol(pbmc),]

META.NAME=str_replace(META[,1],'-','.')
META.NAME=str_replace(META.NAME,'10X','X10X')
META=META[match(colnames(pbmc), META.NAME),]
pbmc@meta.data$type=META$characteristics..cell.cluster

saveRDS(pbmc,file='pbmc.RDS')

#########################







##########################################
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseDentateGyrus_GSE104323/')
pbmc=readRDS(file='pbmc.RDS')




tiff(paste0("IMG/ALL.tiff"),width=5.5,height=5.5,units='in',res=600)
DimPlot(pbmc, group.by='type',label = TRUE)+NoLegend()+NoAxes()
dev.off()



VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= pbmc@reductions$pca@cell.embeddings




tiff(paste0("IMG/VECTOR.1.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.buildGrid(VEC, N=40,SHOW=TRUE)
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










####
#Change MS to NES
NES.EXP = pbmc@assays$RNA@data[which(rownames(pbmc) =='Nes'),]


OUT=vector.buildGrid(VEC, N=40,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)



OUT$VALUE=NES.EXP

tiff(paste0("IMG/NES.1.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.showValue(OUT)
dev.off()

tiff(paste0("IMG/NES.2.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.gridValue(OUT, SHOW=TRUE)
dev.off()

tiff(paste0("IMG/NES.3.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.autoCenter(OUT,UP=0.95,SHOW=TRUE)
dev.off()

tiff(paste0("IMG/NES.4.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)
dev.off()



###########
#Select starting point



OUT=vector.buildGrid(VEC, N=40,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
OUT=vector.gridValue(OUT,SHOW=TRUE)


OUT=vector.selectCenter(OUT)

tiff(paste0("IMG/MSP.2.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)
dev.off()





###########
#Select arrows


OUT=vector.buildGrid(VEC, N=40,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)

OUT=vector.reDrawArrow(OUT, COL=OUT$COL)
OUT=vector.selectRegion(OUT)



SELECT_SCORE=OUT$SELECT_SCORE
SELECT_INDEX=OUT$SELECT_INDEX

EXP=as.matrix(pbmc@assays$RNA@scale.data)[,SELECT_INDEX]


print(nrow(EXP))
#10646
COR=c()
i=1
while(i<=nrow(EXP)){
    this_cor=cor(SELECT_SCORE, EXP[i,],method='spearman')
    COR=c(COR,this_cor)
    if(i %%100==1){print(i)}
    i=i+1}


names(COR)=rownames(EXP)

head(sort(COR),n=10)
# Tubb4a      Aplp1      Stmn4       Plp1      Cryab     Slain1        App       Qdpr        Mbp 
#-0.8363392 -0.8004671 -0.7275945 -0.7040098 -0.6983233 -0.6754448 -0.6709485 -0.6688818 -0.6659579 
#     Cntn2 
#-0.6651465


tail(sort(COR),n=10)
#Rpl22l1      Sox9       Ran     Hmgb3     Hmgb2     Fabp7     H2afv    Slc1a3     Hmgn1     H2afz 
#0.5529178 0.5648042 0.5658522 0.5674023 0.6034539 0.6041718 0.6155595 0.6251123 0.6455395 0.7021575 

ORDER=order(-SELECT_SCORE)
SOX9=EXP[which(rownames(EXP)=='Sox9'),ORDER]
INDEX=1:length(SOX9)
USED=which(SOX9>0)

INDEX=INDEX[USED]
SOX9=SOX9[USED]


SOX9.fit=lm(SOX9~INDEX)


tiff(paste0("IMG/SOX9.tiff"),width=2,height=2,units='in',res=600)
par(mar=c(0,0,0,0))
plot(INDEX, SOX9, pch=16,col='grey70',cex=0.5)
abline(SOX9.fit,col='red',lwd=2)
dev.off()







