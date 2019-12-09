source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/Try/')

#####################################



pbmc=readRDS('pbmc.RDS')
PCA=pbmc@reductions$pca@cell.embeddings

VEC=pbmc@reductions$umap@cell.embeddings



MED=vector.medCurv(PCA)

PCA=PCA[,1:which(MED==max(MED))]


    OUT=vector.buildGrid(VEC, N=10,SHOW=TRUE)
    OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
    OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
    OUT=vector.gridValue(OUT,SHOW=TRUE)
    OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
    OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)






PCA.OUT.1=readRDS('F:/Vector/data/MouseIntestine_GSE92332/PCA.OUT_500.RDS')
PCA.OUT.2=readRDS('F:/Vector/data/MouseOligo_GSE75330/PCA.OUT_500.RDS')

PCA.1=PCA.OUT.1$x
PCA.2=PCA.OUT.2$x


MS.1=vector.calValue(PCA.1)$PCA.RC

COR.1=cor(MS.1)



vector.corCurv <- function(PCA,MAX=1000){
    PCA=PCA
    MAX=MAX
    MS=vector.calValue(PCA)$PCA.RC
    COR=cor(MS)
    MED=c()
    ALL=c()
    i=1
    while(i<=ncol(COR) & i<=MAX){
        this_col=COR[1:i,i]
        ALL=c(ALL,this_col)
        ALL[which(ALL==1)]=NA
        this_med=median(ALL,na.rm = TRUE)
        MED=c(MED,this_med)
        if(i %%100==1){print(i)}
        i=i+1}
    MED[1]=MED[2]
    return(MED)
    }







MS.2=vector.calValue(PCA.2)$PCA.RC

COR.2=cor(MS.2)



LLL.1=c(0)

i=2
while(i<=ncol(COR.1)){
    this_lll=median(unique(COR.1[1:i,1:i]))
    LLL.1=c(LLL.1,this_lll)

    i=i+1}

plot(LLL.1[20:length(LLL.1)])




LLL.2=c(0)

i=2
while(i<=ncol(COR.2)){
    this_lll=median(unique(COR.2[1:i,1:i]))
    LLL.2=c(LLL.2,this_lll)

    i=i+1}

plot(LLL.2[20:length(LLL.2)])






ROT.1=PCA.OUT.1$rotation
ROT.2=PCA.OUT.2$rotation





ROT.1.NAME.UP=toupper(rownames(ROT.1))
ROT.2.NAME.UP=toupper(rownames(ROT.2))

USED.1=which(ROT.1.NAME.UP %in% names(which(table(ROT.1.NAME.UP )==1)))
ROT.1=ROT.1[USED.1,]


USED.2=which(ROT.2.NAME.UP %in% names(which(table(ROT.2.NAME.UP )==1)))
ROT.2=ROT.2[USED.2,]


rownames(ROT.1)=toupper(rownames(ROT.1))
rownames(ROT.2)=toupper(rownames(ROT.2))

colnames(ROT.1)=paste0('D1_',colnames(ROT.1))
colnames(ROT.2)=paste0('D2_',colnames(ROT.2))


ROT.1=ROT.1[,1:500]
ROT.2=ROT.2[,1:500]



ROT.COM=.simple_combine(ROT.1,ROT.2)$combine

COR=cor(ROT.COM,method='spearman')
SSS=apply(abs(COR),2,mean)

USED=which(SSS<quantile(SSS,0.1))


ROT.COM=ROT.COM[,USED]
##########################################







source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/MouseIntestine_GSE92332/')
pbmc=readRDS(file='pbmc_smart.RDS')



ROT.COM=.simple_combine(ROT.1,ROT.2)$combine
COR=cor(ROT.COM,method='spearman')
SSS.1=apply(abs(COR[(ncol(ROT.1)+1):ncol(COR),]),2,mean)
SSS.2=apply(abs(COR[1:ncol(ROT.1),]),2,mean)
SSS=c(SSS.1[1:ncol(ROT.1)],SSS.2[(ncol(ROT.1)+1):ncol(COR)])
USED=which(SSS<quantile(SSS,0.6))


ROT.COM=ROT.COM[,USED]



EXP=as.matrix(pbmc@assays$RNA@data)
RN=toupper(rownames(EXP))
EXP=EXP[which(RN %in% names(which(table(RN)==1))),]

rownames(EXP)=toupper(rownames(EXP))

COM=.simple_combine(EXP,ROT.COM)$combine

EXP=COM[,1:ncol(EXP)]
ROT=COM[,(1+ncol(EXP)):ncol(COM)]

EXP=t(apply(t(EXP),2,scale))

TMP=t(EXP) %*% ROT




VEC=pbmc@reductions$umap@cell.embeddings
rownames(VEC)=colnames(pbmc)
PCA= TMP#pbmc@reductions$pca@cell.embeddings





OUT=vector.buildGrid(VEC, N=40,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)












