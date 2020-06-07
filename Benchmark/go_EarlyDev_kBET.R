source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

setwd('F:/Vector/data/kBET/')

#/home/zy/single_cell/mouse_EED
DATA=read.csv('Counts.csv',header=TRUE,row.names=1,sep=',')

DATA=t(DATA)


META=read.csv('pData.csv',header=TRUE,row.names=1,sep=',')
TYPE=META$sample
BATCH=META$batch

saveRDS(DATA,file='DATA.RDS')
saveRDS(META,file='META.RDS')

####################################################################
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/kBET/')

DATA=readRDS(file='DATA.RDS')
META=readRDS(file='META.RDS')

TYPE=META$sample
BATCH=META$batch

        library(sva)
        library(limma)
        pheno = data.frame(batch=as.matrix(BATCH))
        orig.data=DATA
        used.gene.index=1:nrow(DATA)#which(rownames(orig.data) %in% VARG)
        edata = as.matrix(orig.data)[used.gene.index,]
        batch = pheno$batch
        modcombat = model.matrix(~1, data=pheno)
        combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
        rownames(combat_edata)=rownames(edata)
        colnames(combat_edata)=colnames(edata)
        combat_edata=as.matrix(combat_edata)
        combat_edata[which(combat_edata<0)]=0



pbmc <- CreateSeuratObject(counts = combat_edata, project = "pbmc3k", min.cells = 0, min.features = 0)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs = 150)




pbmc@meta.data$type=TYPE


pbmc <- RunUMAP(pbmc, dims = 1:150)
DimPlot(pbmc, reduction = "umap",group.by='type',label=TRUE)+NoLegend()+NoAxes()



saveRDS(pbmc,file='pbmc.RDS')
######################################


source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
setwd('F:/Vector/data/kBET/')

pbmc=readRDS(file='pbmc.RDS')

PCA=pbmc@reductions$pca@cell.embeddings[,c(1:150)]
VEC=pbmc@reductions$umap@cell.embeddings#pbmc@reductions$pca@cell.embeddings[,c(1, 2)]


####
OUT=vector.buildGrid(VEC, N=20,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, OL=2,COL=OUT$COL)



################################################################
VVV=OUT$VALUE
TTT=pbmc@meta.data$type#rep(NA, length(VVV))
#DATA=pbmc@assays$RNA@counts
#TTT[which(DATA[which(rownames(DATA)=='Sox10'),] >0)]='Sox10 +'
#TTT[which(DATA[which(rownames(DATA)=='Th'),] >0)]='Th +'

#USED=which(!is.na(TTT))
#TTT=TTT[USED]
#VVV=VVV[USED]

TTT=as.factor(TTT)
TTT=factor(TTT , levels=c('2-cell','4-cell','8-cell','16-cell','32-cell'))

boxplot(VVV~TTT,outline=FALSE,xlab='',ylab='',las=2)
####################################################################



tiff(paste0("IMG/NEW_VECTOR_TRY.6.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL,AL=40,BD=F,OL=2,AW=1.5,AC='black',CEX=1)
dev.off()


tiff(paste0("IMG/NEW_VECTOR.6.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL,OL=2,AL=40,CEX=1)
dev.off()
###################




OUT=vector.buildGrid(VEC, N=20,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)







tiff(paste0("IMG/VECTOR.1.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.buildGrid(VEC, N=20,SHOW=TRUE)
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


tiff(paste0("IMG/TYPE.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
DimPlot(pbmc, reduction = "umap",group.by='type',label=TRUE,label.size=6)+NoLegend()+NoAxes()
dev.off()












vector.showValue<-function(OUT){
    VEC=OUT$VEC
    VALUE=OUT$VALUE
    #################################
    N.VALUE=.normX(VALUE)
    COL=vector.vcol(N.VALUE, c(0,0.5,1),c('#009FFF','#FFF200','#ec2F4B'))
    plot(VEC,col=COL,pch=16,cex=1)
    return(OUT)
    }


tiff(paste0("IMG/VECTOR.3.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.showValue(OUT)
dev.off()



vector.autoCenter <- function(OUT, UP=0.9, SHOW=TRUE){
    #####################
    OUT=OUT
    SHOW=SHOW
    UP=UP
    #############################
    DIST=OUT$DIST
    USED=OUT$USED
    USED_NAME=OUT$USED_NAME
    CENTER_VALUE=OUT$CENTER_VALUE
    CENTER_VEC=OUT$CENTER_VEC
    CENTER_INDEX=OUT$CENTER_INDEX
    INDEX_LIST=OUT$INDEX_LIST
    ############
    
    
    USED_CENTER_VALUE=CENTER_VALUE[USED]
    USED_CENTER_VEC=CENTER_VEC[USED,] 
    
       
    HIGH=USED[which(USED_CENTER_VALUE>=quantile(USED_CENTER_VALUE,UP))]
    HIGH_NAME=USED_NAME[which(USED_CENTER_VALUE>=quantile(USED_CENTER_VALUE,UP))]
    LOW_NAME=USED_NAME[which(USED_CENTER_VALUE<quantile(USED_CENTER_VALUE,UP))]
    #plot(OUT$CENTER_VEC[USED,])
    #points(CENTER_VEC[HIGH,], col='red',pch=16)
    #plot(CENTER_VEC[HIGH,], col='red')
    
    SUB=induced_subgraph(OUT$GRAPH, HIGH_NAME)
    SUB_CPT=components(SUB)
    CLUSTER=list()
    LENGTH=c()
    PCH=rep(1,length(HIGH))
    DIST_COR=c()
    DIST_MEAN=c()
    
    i=1
    while(i<=SUB_CPT$no){
        this_name=names(which(SUB_CPT$membership==i))
        this_index=as.numeric(str_replace(this_name,'P',''))  
        PCH[which(HIGH %in% this_index)]=as.character(i)
        LENGTH=c(LENGTH, length(this_index))
        if(length(this_index)==1){
            this_dist=DIST[USED, this_index]
        }else{
            this_dist=apply(DIST[USED, this_index],1,mean)
        }
        this_cor=cor(this_dist, CENTER_VALUE[USED],method='spearman')
        DIST_COR=c(DIST_COR,this_cor)
        
        
        DIST_MEAN=c(DIST_MEAN, mean(this_dist))
        CLUSTER=c(CLUSTER,list(this_index))       
        i=i+1}
    
    ####################
    
    #SELECT=which(DIST_COR==min(DIST_COR))[1]
    #SELECT=which(rank(DIST_MEAN)*rank(DIST_COR) == min(rank(DIST_MEAN)*rank(DIST_COR)) )
    
    TMP=10^LENGTH - DIST_COR
    #TMP=rank(-DIST_COR) * rank(LENGTH )
    SELECT=which(TMP==max(TMP))[1]
    #SELECT=which(LENGTH==max(LENGTH))[1]
    #SELECT=which( rank(-DIST_COR) * rank(LENGTH) == max( rank(-DIST_COR) * rank(LENGTH)  ) )[1]
    
    SUMMIT=CLUSTER[[SELECT]]
    #print(SELECT)
    #plot(OUT$CENTER_VEC[USED,])
    #points(OUT$CENTER_VEC[CLUSTER[[9]],],pch=16,col='red')
    
    
    SCORE=c()
    i=1
    while(i<=length(USED_NAME)){
        this_name=USED_NAME[i]
        this_dist=DIST[which(colnames(DIST)==this_name),
                       which(rownames(DIST) %in% paste0('P',SUMMIT))]
        this_dist=this_dist
        #this_value=CENTER_VALUE[USED]
        this_score=min(this_dist)#sum(rank(-this_dist) * rank(this_value))
        #this_cor=cor(-this_value, this_dist)#,method='spearman')
        SCORE=c(SCORE, this_score)
        i=i+1}
    SCORE=max(SCORE)-SCORE
    ##########################################
    
    VALUE=SCORE
    #plot(OUT$VEC, col='grey70',pch=16)
    N.VALUE=(VALUE-min(VALUE))/(max(VALUE)-min(VALUE))
    COL=vector.vcol(N.VALUE, c(0,0.5,1),c('#009FFF','#FFF200','#ec2F4B'))
    
    
    ###################################################
    OUT$COL=rep('grey70',nrow(OUT$VEC))
    OUT$ORIG.COL=rep('grey70',nrow(OUT$VEC))
    OUT$P.SCORE=rep(0,nrow(OUT$VEC))
    i=1
    while(i<=length(USED)){
        this_index=INDEX_LIST[[USED[i]]]
        OUT$COL[this_index]=COL[i]
        OUT$ORIG.COL[this_index]=OUT$ORIG.CENTER.COL[USED][i]
        OUT$P.SCORE[this_index]=SCORE[i]
        i=i+1
        }
    ###########
    
    
    if(SHOW==TRUE){
        #plot(OUT$VEC, col=OUT$COL, pch=16, cex=0.5 )
        plot(OUT$VEC, col=OUT$ORIG.COL, pch=16, cex=1)
        text(CENTER_VEC[HIGH,1], CENTER_VEC[HIGH,2],labels=PCH,cex=1,pos=2)
        points(CENTER_VEC[HIGH,1], CENTER_VEC[HIGH,2], col='black',pch=16,cex=1)
        points(CENTER_VEC[SUMMIT,1],CENTER_VEC[SUMMIT,2], col='black',pch=16,cex=2)
        points(CENTER_VEC[SUMMIT,1],CENTER_VEC[SUMMIT,2], col='red',pch=16,cex=1.5)  
        }
 
    ######################
    ############
    OUT$SCORE=SCORE
    OUT$SUMMIT=SUMMIT
    OUT$CLUSTER=CLUSTER
    OUT$LENGTH=LENGTH
    OUT$PCH=PCH
    OUT$DIST_COR=DIST_COR
    #OUT$DIST_MEAN=DIST_MEAN
    ################################
    
    #######################
    return(OUT)
    }
tiff(paste0("IMG/VECTOR.5.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
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
        plot(ALL_VEC,col=COL,pch=16,cex=1)
        }
    
    N.SCORE=.normX(SCORE)
    SCORE.COL=vector.vcol(N.SCORE, c(0,0.5,1),c('#009FFF','#FFF200','#ec2F4B'))
    
    A1_VEC=c()
    A2_VEC=c()
    A_LENGTH=c()
    one=1
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


tiff(paste0("IMG/VECTOR.6.tiff"),width=4,height=4,units='in',res=600)
par(mar=c(0,0,0,0))
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)
dev.off()









