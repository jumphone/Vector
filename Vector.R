###################
# Author: Feng Zhang
# Date: 11/25, 2019

##################
library('circlize')
library('gatepoints')
library('stringr')
library('igraph')
library('gmodels')
##################





vector.SeuratSelect <- function(pbmc){
    pbmc=pbmc
    ##################
    ppp=DimPlot(pbmc, reduction.use='umap', pt.size=0.5)
    used.cells <- CellSelector(plot = ppp)
    #################
    return(used.cells)
    }

vector.SeuratAddMetaByCell <- function(pbmc, used.cells){
    pbmc=pbmc
    used.cells=used.cells
    SELECT=rep('NO',ncol(pbmc))
    SELECT[which(colnames(pbmc) %in% used.cells)]='YES'
    pbmc@meta.data$select=SELECT
    return(pbmc)
    }




.normX <- function(x){
    y=(x-min(x))/(max(x)-min(x))
    return(y)
    }



vector.lcol <- function(TAG){
    TAG=as.factor(TAG)
    require(scales)
    my_color_palette <- hue_pal()(length(unique(TAG)))
    COL=my_color_palette[TAG]
    return(COL)
    }


vector.vcol<-function(VALUE, CV, CN){
    VALUE=VALUE 
    library('circlize')
    CRF=colorRamp2(CV, CN)
    COL=CRF(VALUE)
    return(COL)
    }




vector.showValue<-function(OUT){
    VEC=OUT$VEC
    VALUE=OUT$VALUE
    #################################
    N.VALUE=.normX(VALUE)
    COL=vector.vcol(N.VALUE, c(0,0.5,1),c('#009FFF','#FFF200','#ec2F4B'))
    plot(VEC,col=COL,pch=16,cex=0.5)
    return(OUT)
    }

vector.calValue <- function(PCA){
    OUT=list()
    PCA=PCA
    PCA.RC=apply(apply(PCA,2,rank), 2, .normX)
    PCA.RC=abs(PCA.RC-0.5)   
    VALUE=apply(PCA.RC,1,mean)
    #VALUE=rank(VALUE)
    ###########################
    OUT$VALUE=VALUE
    OUT$PCA.RC=PCA.RC
    return(OUT)
    }



vector.getValue <-function(OUT, PCA, SHOW=TRUE){
    OUT=OUT
    PCA=PCA
    SHOW=SHOW
    ############
    VALUE.OUT=vector.calValue(PCA)
    ###############
    OUT$VALUE=VALUE.OUT$VALUE
    OUT$PCA=PCA
    OUT$PCA.RC=VALUE.OUT$PCA.RC
    ###############
    if(SHOW==TRUE){
        vector.showValue(OUT)
        }
    #####
    return(OUT)
    }










vector.buildGrid <- function(VEC, N=30, SHOW=TRUE, COL='grey70'){
    #############
    VEC.E=VEC
    COL=COL
    SHOW=SHOW
    delta=0.000001
    N=N
    ##############
    if(SHOW==TRUE){
        plot(VEC.E,col=COL,pch=16,cex=0.2)
        }
    ###############
    X=VEC[,1]
    Y=VEC[,2]
    
    this_step_x=( max(X) - min(X) - delta)/N
    this_step_y=( max(Y) - min(Y) - delta)/N
    NUM_CUT=1
    ############
    
    INDEX_LIST=list()
    CENTER_LIST=list()

    this_x=min(X)
    while(this_x<max(X)){
        this_y=min(Y)
        while(this_y<max(Y)){  
            ################
            this_in_index = which(VEC.E[,1]>=this_x &  VEC.E[,1] <this_x+this_step_x &
                                  VEC.E[,2]>=this_y &  VEC.E[,2] <this_y+this_step_y)                                         
            this_center=c(this_x+this_step_x/2,this_y+this_step_y/2)                  
            if(length(this_in_index)>=NUM_CUT){
                #####################
                INDEX_LIST=c(INDEX_LIST, list(this_in_index))
                CENTER_LIST=c(CENTER_LIST, list(this_center))
                ######################
                if(SHOW==TRUE){
                    points(this_center[1],this_center[2],col= 'black',pch=16,cex=0.5)
                    }
                ######################
                }
                 
            this_y=this_y+this_step_y
            }
        #print(this_x)
        this_x=this_x+this_step_x
        }
    ############################################
    OUT=list()
    OUT$VEC=VEC.E
    OUT$INDEX_LIST=INDEX_LIST
    OUT$CENTER_LIST=CENTER_LIST
    OUT$this_step_x=this_step_x
    OUT$this_step_y=this_step_y
    ########################################
    return(OUT)
    }




vector.buildNet<-function(OUT,CUT=1,SHOW=TRUE,COL='grey70'){
    library(igraph)
    library(stringr)
    OUT=OUT
    SHOW=SHOW
    CUT=CUT
    COL=COL
    VEC=OUT$VEC
    INDEX_LIST=OUT$INDEX_LIST
    CENTER_LIST=OUT$CENTER_LIST
    this_step_x=OUT$this_step_x
    this_step_y=OUT$this_step_y
    delta=0.00001
    ###############################
    if(SHOW==TRUE){
        plot(VEC,col=COL,pch=16,cex=0.2)
        }
    #############################
    CENTER_VEC=c()
    i=1
    while(i<=length(CENTER_LIST)){
        CENTER_VEC=cbind(CENTER_VEC, CENTER_LIST[[i]])
    
        i=i+1}
    CENTER_VEC=t(CENTER_VEC)
    OUT$CENTER_VEC=CENTER_VEC
    ##############################
    p1=c()
    p2=c()
    
    #############################
    CNUM=length(CENTER_LIST)
    i=1
    while(i<=CNUM){
    
        this_p1_loc=CENTER_LIST[[i]] 
        this_p1=paste0('P',as.character(i))
    
        used_j=which( ( abs(CENTER_VEC[,1]-this_p1_loc[1]) <= this_step_x+delta) 
                     & ( abs(CENTER_VEC[,2]-this_p1_loc[2]) <= this_step_y+delta) )
        
        for(j in used_j){
            this_p2=paste0('P',as.character(j))
            
            this_p2_loc=CENTER_LIST[[j]]
 
            ######################
            
            if(length(INDEX_LIST[[i]])>=CUT & 
               length(INDEX_LIST[[j]])>=CUT &
               this_p1 != this_p2){
                
                p1=c(p1,this_p1)
                p2=c(p2,this_p2) 
                if(SHOW==TRUE){
                    segments(x0=this_p1_loc[1],x1=this_p2_loc[1],
                           y0=this_p1_loc[2],y1=this_p2_loc[2],
                       
                           col='black')
                    } 
                }
    
            ############
            }
       
            #if(i%%10==1){print(paste0(i,'/',CNUM))}
            i=i+1
       }
    ##########################  
    
    
    
    TAG=c()
    i=1
    while(i<=length(p1)){
        this_p1=p1[i]
        this_p2=p2[i]
        sorted_pair=sort(c(this_p1,this_p2))
        this_tag=paste0(sorted_pair[1],'|',sorted_pair[2])
        TAG=c(TAG, this_tag)
        i=i+1}
    TAG=unique(TAG)
    #NET=tapply(NET,1,sort)
    p1=c()
    p2=c()
    i=1
    while(i<=length(TAG)){
        this_p1=strsplit(TAG[i],'\\|')[[1]][1]
        this_p2=strsplit(TAG[i],'\\|')[[1]][2]
        p1=c(p1,this_p1)
        p2=c(p2,this_p2)
        i=i+1}
    
    ##########################
    
    library(igraph)
    OUT$p1=p1
    OUT$p2=p2
    NET = cbind(p1,p2) 
    g <- make_graph(t(NET),directed = FALSE)
    ALLNAME=paste0('P',1:CNUM)
    ADD=ALLNAME[which(! ALLNAME %in% as_ids(V(g)))]
    g <- g + ADD
    
    ##########################
    DIST=distances(g, v = V(g), to = V(g), mode = c("all"))
    library(stringr)
    DIST.NUM=as.numeric(str_replace(colnames(DIST),'P',''))
    DIST=DIST[order(DIST.NUM),order(DIST.NUM)]
    ###########################
    library(igraph)
    CPT=components(g)
    MAXC=which(CPT$csize==max(CPT$csize))[1]
    ############
    library(stringr)
    ###############
    USED_NAME=names(which(CPT$membership==MAXC))
    USED=as.numeric(str_replace(USED_NAME,'P',''))
    
    USED_NAME=USED_NAME[order(USED)]
    USED=USED[order(USED)]
    if(SHOW==TRUE){
        points(OUT$CENTER_VEC[USED,],col='red',pch=16, cex=0.5)
        }
    USED_INDEX=c()
    i=1
    while(i<=length(USED)){
        USED_INDEX=c(USED_INDEX,INDEX_LIST[[USED[i]]])
        i=i+1}
    
    
    ###################
    OUT$GRAPH=g
    OUT$DIST=DIST
    OUT$USED=USED
    OUT$USED_NAME=USED_NAME
    OUT$USED_INDEX=USED_INDEX
    return(OUT)
    ###########################
    }








vector.gridValue <- function(OUT, SHOW=TRUE){
    OUT=OUT
    INDEX_LIST=OUT$INDEX_LIST
    VALUE=OUT$VALUE
    SHOW=SHOW
    USED=OUT$USED
    ######################
    
    ################
    CENTER_VALUE=c()
    i=1
    while(i<=length(INDEX_LIST)){
        this_value = mean(VALUE[INDEX_LIST[[i]]])
        CENTER_VALUE=c(CENTER_VALUE,this_value)
        i=i+1}
    ############
    CENTER_VEC=OUT$CENTER_VEC

    ##################
    VALUE=VALUE
    N.VALUE=(VALUE-min(VALUE))/(max(VALUE)-min(VALUE))
    VALUE.COL=vector.vcol(N.VALUE, c(0,0.5,1),c('#009FFF','#FFF200','#ec2F4B'))
    ####################
    
    #################################
    VALUE=CENTER_VALUE
    N.VALUE=(VALUE-min(VALUE))/(max(VALUE)-min(VALUE))
    COL=vector.vcol(N.VALUE, c(0,0.5,1),c('#009FFF','#FFF200','#ec2F4B'))
    
    ################
    if(SHOW==TRUE){
        
        plot(OUT$VEC,col='grey80',pch=16,cex=0.5)
        #points(OUT$CENTER_VEC,col=COL,pch=16)
        points(OUT$CENTER_VEC[USED,],col=COL[USED],pch=15, cex=1.5)
        }
    ################    
    OUT$CENTER_VALUE=CENTER_VALUE     
    OUT$ORIG.CENTER.COL=COL
    #OUT$VALUE=VALUE
    OUT$VALUE.COL=VALUE.COL
    return(OUT)
    }



vector.nonCenter<-function(OUT){
    OUT=OUT
    OUT$SCORE=OUT$CENTER_VALUE[OUT$USED]
    OUT$COL=OUT$VALUE.COL
    return(OUT)
    }




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
    PS=c()
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
        PS=c(PS, this_score)
        i=i+1}
    SCORE=max(SCORE)-SCORE
    ##########################################
    
    VALUE=SCORE
    #plot(OUT$VEC, col='grey70',pch=16)
    N.VALUE=(VALUE-min(VALUE))/(max(VALUE)-min(VALUE))
    #COL=vector.vcol(N.VALUE, c(0,0.5,1),c('#009FFF','#FFF200','#ec2F4B')) # too strong
    #COL=vector.vcol(N.VALUE,c(0,0.5,1),c('#009FFF','#FFF200','#ffdde1')) # too weak
    COL=vector.vcol(N.VALUE,c(0,0.5,1),c('#009FFF','#FFF200','#ee9ca7'))
    
    ###################################################
    OUT$COL=rep('grey70',nrow(OUT$VEC))
    OUT$ORIG.COL=rep('grey70',nrow(OUT$VEC))
    OUT$P.SCORE=rep(0,nrow(OUT$VEC))
    OUT$P.PS=rep(NA,nrow(OUT$VEC))
    i=1
    while(i<=length(USED)){
        this_index=INDEX_LIST[[USED[i]]]
        OUT$COL[this_index]=COL[i]
        OUT$ORIG.COL[this_index]=OUT$ORIG.CENTER.COL[USED][i]
        OUT$P.SCORE[this_index]=SCORE[i]
        OUT$P.PS[this_index]=PS[i]
        i=i+1
        }
    ###########
    
    
    if(SHOW==TRUE){
        #plot(OUT$VEC, col=OUT$COL, pch=16, cex=0.5 )
        plot(OUT$VEC, col=OUT$ORIG.COL, pch=16, cex=0.5)
        text(CENTER_VEC[HIGH,1], CENTER_VEC[HIGH,2],labels=PCH,cex=1,pos=2)
        points(CENTER_VEC[HIGH,1], CENTER_VEC[HIGH,2], col='black',pch=16,cex=1)
        points(CENTER_VEC[SUMMIT,1],CENTER_VEC[SUMMIT,2], col='black',pch=16,cex=1.5)
        points(CENTER_VEC[SUMMIT,1],CENTER_VEC[SUMMIT,2], col='red',pch=16,cex=1)  
        }
 
    ######################
    ############
    OUT$SCORE=SCORE
    OUT$SUMMIT=SUMMIT
    OUT$CLUSTER=CLUSTER
    OUT$LENGTH=LENGTH
    OUT$PCH=PCH
    OUT$DIST_COR=DIST_COR
    OUT$PS=PS
    #OUT$DIST_MEAN=DIST_MEAN
    ################################
    
    #######################
    return(OUT)
    }




vector.drawArrow <- function(OUT, P=0.9, SHOW=TRUE, COL='grey70',OL=1.5,AL=60,CEX=0.5,AW=2, BD=TRUE, AC='grey20', SHOW.SUMMIT=TRUE){
    ################
    AW=AW
    BD=BD
    AC=AC
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
    OL=OL
    CEX=CEX
    SHOW.SUMMIT=SHOW.SUMMIT
    #####################
    one=min(dist(USED_CENTER_VEC)) * OL
    #####################
    ###################
    .norm_one <-function(x,one=1){
        one=one
        if(var(x)!=0){
            x=x/sqrt(sum(x^2)) * one }
        return(x)
        }
    
    DIV=1/P
    ####################
    if(SHOW==TRUE){
        #plot(ALL_VEC,col=COL,pch=16,cex=CEX)
        if(BD==TRUE){
            plot(ALL_VEC,col=COL,pch=16,cex=CEX, 
                  xlim=c(min(ALL_VEC[,1])-one, max(ALL_VEC[,1])+one ),
                 ylim=c(min(ALL_VEC[,2])-one, max(ALL_VEC[,2])+one ) 
                )
            }else{
            plot(ALL_VEC,col=COL,pch=16,cex=CEX, 
                 xlim=c(min(ALL_VEC[,1])-one, max(ALL_VEC[,1])+one ),
                 ylim=c(min(ALL_VEC[,2])-one, max(ALL_VEC[,2])+one ) , yaxt="n", axes=F
                )
            }
        }
    ##########################
    N.SCORE=.normX(SCORE)
    SCORE.COL=vector.vcol(N.SCORE, c(0,0.5,1),c('#009FFF','#FFF200','#ec2F4B'))
    
    A1_VEC=c()
    A2_VEC=c()
    A_LENGTH=c()
    #####################
    #one=min(dist(USED_CENTER_VEC)) * OL
    #####################
    i=1
    while(i<=length(USED)){
        this_p1_loc=USED_CENTER_VEC[i,]
        
        vector_list=cbind(USED_CENTER_VEC[,1]-this_p1_loc[1],USED_CENTER_VEC[,2]-this_p1_loc[2])
        vector_list_norm=t(apply(vector_list,1,.norm_one, one))
        
        vector_weight_1= DIV^-(rank(USED_DIST[i,])-1)    # equals to P ^ (rank(USED_DIST[i,])-1)
        vector_weight_2= SCORE[i]-SCORE                  # PS - PS[i]; "PS_i - PS_j" in  the paper; In paper, PS[i] is "PS_j".
        
        vector_weight = vector_weight_1 * vector_weight_2        
        vector_weight = vector_weight/sum(abs(vector_weight))
        
        final_vec=t(vector_list_norm) %*% vector_weight
        
        this_p2_loc=c(this_p1_loc[1]+final_vec[1],this_p1_loc[2]+final_vec[2])
                
        #plot(ALL_VEC,col=COL,pch=16,cex=0.2)
        #this_arrow_length=0.1*sqrt(sum(final_vec^2))
        #this_arrow_length=sqrt(sum(final_vec^2))#0.1 #* (1+sqrt(sum(final_vec^2)))
        this_arrow_length=dev.size()[1]/AL *  sqrt(sum(final_vec^2))/one   # * sqrt(sum(final_vec^2)) #0.25
        if(SHOW==TRUE){
            this_arrow_col=AC#'grey20'
            #if(AC==TRUE){this_arrow_col=SCORE.COL[i]}
            arrows(x0=this_p1_loc[1],y0=this_p1_loc[2],
                   x1=this_p2_loc[1],y1=this_p2_loc[2],
                   lwd=AW, length=this_arrow_length,
                   col=this_arrow_col
                   #col=this_arrow_col#SCORE.COL[i]
                   )
            
            }
        A1_VEC=cbind(A1_VEC,this_p1_loc)
        A2_VEC=cbind(A2_VEC,this_p2_loc)
        A_LENGTH=c(A_LENGTH,this_arrow_length)
        i=i+1}
    
    ###############################
    if(SHOW==TRUE & SHOW.SUMMIT==TRUE){   
        #points(OUT$CENTER_VEC[OUT$SUMMIT,],col='black',pch=16,cex=1)
        #points(OUT$CENTER_VEC[OUT$SUMMIT,],col='red',pch=16,cex=0.5)
        X1=min(OUT$CENTER_VEC[OUT$SUMMIT,1])-one/10
        X2=max(OUT$CENTER_VEC[OUT$SUMMIT,1])+one/10
        Y1=min(OUT$CENTER_VEC[OUT$SUMMIT,2])-one/10
        Y2=max(OUT$CENTER_VEC[OUT$SUMMIT,2])+one/10
        
        rect(xleft=X1, ybottom=Y1, xright=X2, ytop=Y2, angle = 45,
        col = NA, border = '#009FFF', lty = 1, lwd = 2)
        #rect(xleft=X1, ybottom=Y1, xright=X2, ytop=Y2, angle = 45,
        #col = NA, border = 'black', lty = 2, lwd = 2)
        }
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




#################################################################
# Manually do something

vector.selectPoint<-function(VEC,CEX=0.5){
    VEC=VEC
    CEX=CEX
    ############
    library(gatepoints)
    selectedPoints <- fhs(VEC, pch=16,col='red3',cex=CEX,mark = TRUE)
    return(selectedPoints)
    }


vector.selectCenter <- function(OUT){
    #####################
    OUT=OUT
    SHOW=TRUE
    #############################
    DIST=OUT$DIST
    USED=OUT$USED
    USED_NAME=OUT$USED_NAME
    CENTER_VALUE=OUT$CENTER_VALUE
    CENTER_VEC=OUT$CENTER_VEC
    CENTER_INDEX=OUT$CENTER_INDEX
    INDEX_LIST=OUT$INDEX_LIST
    ############
    
    plot(OUT$VEC, col='grey80',pch=16,cex=0.5)
    points(OUT$CENTER_VEC[USED,], col=OUT$ORIG.CENTER.COL[USED], pch=16, cex=1)
    SELECT_NAME=vector.selectPoint(OUT$CENTER_VEC[USED,],CEX=1)
    SUMMIT=USED[as.numeric(SELECT_NAME)]
    #print(SELECT)
    #plot(OUT$CENTER_VEC[USED,])
    #points(OUT$CENTER_VEC[CLUSTER[[9]],],pch=16,col='red')
    
    
    SCORE=c()
    
    i=1
    while(i<=length(USED_NAME)){
        this_name=USED_NAME[i]
        this_dist=DIST[which(colnames(DIST)==this_name),which(rownames(DIST) %in% paste0('P',SUMMIT))]
        this_dist=this_dist
        #this_value=CENTER_VALUE[USED]
        this_score=min(this_dist)#sum(rank(-this_dist) * rank(this_value))
        #this_cor=cor(-this_value, this_dist)#,method='spearman')
        SCORE=c(SCORE, this_score)
        i=i+1}
    
    PS=SCORE
    
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
    OUT$P.PS=rep(NA,nrow(OUT$VEC))
    i=1
    while(i<=length(USED)){
        this_index=INDEX_LIST[[USED[i]]]
        OUT$COL[this_index]=COL[i]
        OUT$ORIG.COL[this_index]=OUT$ORIG.CENTER.COL[USED][i]
        OUT$P.SCORE[this_index]=SCORE[i]
        OUT$P.PS[this_index]=PS[i]
        i=i+1
        }
    ###########
    
    
    if(SHOW==TRUE){
        plot(OUT$VEC, col=OUT$COL, pch=16, cex=0.5 )
        points(CENTER_VEC[SUMMIT,], col='black',pch=16,cex=1.5)
        points(CENTER_VEC[SUMMIT,], col='red',pch=16,cex=1)  
        }
    
    ######################
    ############
    OUT$SCORE=SCORE
    OUT$PS=PS
    OUT$SUMMIT=SUMMIT
    ################################
    
    #######################
    return(OUT)
    }





vector.reDrawArrow <- function(OUT, COL='grey70'){
    ################
    A1_VEC=OUT$A1_VEC
    A2_VEC=OUT$A2_VEC
    A_LENGTH=OUT$A_LENGTH
    #A_COL=OUT$A_COL
    VEC=OUT$VEC
    COL=COL
    #####################
    plot(x=VEC[,1],y=VEC[,2], col=COL,cex=0.5, pch=16)
    i=1
    while(i<=length(A_LENGTH)){
        arrows(x0=A1_VEC[i,1],y0=A1_VEC[i,2],
                   x1=A2_VEC[i,1],y1=A2_VEC[i,2],
                   lwd=2, length=A_LENGTH[i],
                   col='black'
                   #col=A_COL[i]
                   )
        i=i+1
        }
    ############################
    return(OUT)
    ############################
    }



vector.selectRegion <- function(OUT){
    #######################
    P.SCORE=OUT$P.SCORE
    ################
    A1_VEC=OUT$A1_VEC
    A2_VEC=OUT$A2_VEC
    A_LENGTH=OUT$A_LENGTH
    VEC=OUT$VEC
    P.PS=OUT$P.PS
    #####################
    CENTER_LIST=OUT$CENTER_LIST
    INDEX_LIST=OUT$INDEX_LIST
    USED=OUT$USED  
    ################################
    SELECT_NAME=vector.selectPoint(VEC,CEX=0.1)
    #########################
    SELECT_INDEX=which(rownames(VEC) %in% SELECT_NAME)
    SELECT_INDEX=SELECT_INDEX[which(SELECT_INDEX %in% OUT$USED_INDEX)]
    #########################
    #########################
    A_USED=c()
    i=1
    while(i<=length(USED)){
        if( length(which(INDEX_LIST[[USED[i]]] %in% SELECT_INDEX )) > 0 ){
            A_USED=c(A_USED, i)
            }
                         
        i=i+1}
    
    ##########################
    #Draw new
    COL=OUT$COL
    COL[which(!rownames(VEC) %in% SELECT_NAME)]='grey70'
    plot(x=VEC[,1],y=VEC[,2], col=COL,cex=0.5, pch=16)
    #points(x=VEC[,1],y=VEC[,2], col=COL,cex=0.5, pch=16)
    #########################
    i=1
    while(i<=length(A_LENGTH)){
        arrows(x0=A1_VEC[i,1],y0=A1_VEC[i,2],
                   x1=A2_VEC[i,1],y1=A2_VEC[i,2],
                   lwd=2, length=A_LENGTH[i],
                   col='black'
                   )
        i=i+1
        }
    ##########################
    #points(x=VEC[,1],y=VEC[,2], col='grey70',cex=0.5, pch=16)
    for(i in A_USED){
        arrows(x0=A1_VEC[i,1],y0=A1_VEC[i,2],
                   x1=A2_VEC[i,1],y1=A2_VEC[i,2],
                   lwd=2, length=A_LENGTH[i],
                   col='red'
                   )
        
        }
    #########################
    OUT$A_USED=A_USED
    OUT$SELECT_NAME=SELECT_NAME
    OUT$SELECT_INDEX=SELECT_INDEX
    OUT$SELECT_SCORE=P.SCORE[SELECT_INDEX]
    OUT$SELECT_PS=P.PS[SELECT_INDEX]
    ########################
    return(OUT)
    }



vector.reDrawRegion <- function(OUT){
    #######################
    P.SCORE=OUT$P.SCORE
    ################
    A1_VEC=OUT$A1_VEC
    A2_VEC=OUT$A2_VEC
    A_LENGTH=OUT$A_LENGTH
    VEC=OUT$VEC
    #####################
    CENTER_LIST=OUT$CENTER_LIST
    INDEX_LIST=OUT$INDEX_LIST
    USED=OUT$USED  
    ################################
    SELECT_NAME=OUT$SELECT_NAME#vector.selectPoint(VEC,CEX=0.1)
    #########################
    SELECT_INDEX=which(rownames(VEC) %in% SELECT_NAME)
    #########################
    #########################
    A_USED=c()
    i=1
    while(i<=length(USED)){
        if( length(which(INDEX_LIST[[USED[i]]] %in% SELECT_INDEX )) > 0 ){
            A_USED=c(A_USED, i)
            }
                         
        i=i+1}
    
    ##########################
    #Draw new
    COL=OUT$COL
    COL[which(!rownames(VEC) %in% SELECT_NAME)]='grey70'
    plot(x=VEC[,1],y=VEC[,2], col=COL,cex=0.5, pch=16)
    #points(x=VEC[,1],y=VEC[,2], col=COL,cex=0.5, pch=16)
    #########################
    i=1
    while(i<=length(A_LENGTH)){
        arrows(x0=A1_VEC[i,1],y0=A1_VEC[i,2],
                   x1=A2_VEC[i,1],y1=A2_VEC[i,2],
                   lwd=2, length=A_LENGTH[i],
                   col='black'
                   )
        i=i+1
        }
    ##########################
    #points(x=VEC[,1],y=VEC[,2], col='grey70',cex=0.5, pch=16)
    for(i in A_USED){
        arrows(x0=A1_VEC[i,1],y0=A1_VEC[i,2],
                   x1=A2_VEC[i,1],y1=A2_VEC[i,2],
                   lwd=2, length=A_LENGTH[i],
                   col='red'
                   )
        
        }
    #########################
    OUT$A_USED=A_USED
    OUT$SELECT_NAME=SELECT_NAME
    OUT$SELECT_INDEX=SELECT_INDEX
    OUT$SELECT_SCORE=P.SCORE[SELECT_INDEX]
    OUT$SELECT_PS=P.PS[SELECT_INDEX]
    ########################
    return(OUT)
    }





######################################

vector.RPPCA <- function(PCA){
    PCA=PCA
    ##################
    R.PCA=apply(PCA,2,rank)
    library(gmodels)
    D=t(R.PCA)
    PCA.OUT=fast.prcomp(t(D), retx = TRUE, center = FALSE, scale. = FALSE, tol = NULL)
    PPCA=PCA.OUT$x
    return(PPCA)
    }




vector.SeuratPCA <-function(pbmc, CUT=0.5){
    pbmc=pbmc
    CUT=CUT
    #####################
    library(gmodels)
    D=as.matrix(pbmc@assays$RNA@scale.data)
    PCA.OUT=fast.prcomp(t(D), retx = TRUE, center = FALSE, scale. = FALSE, tol = NULL)
    EXP=(cumsum(PCA.OUT$sdev^2)/sum(PCA.OUT$sdev^2))
    N=min(which( cumsum(PCA.OUT$sdev^2)/sum(PCA.OUT$sdev^2) > CUT))
    #####################
    PCA.OUT$EXP=EXP
    PCA.OUT$CUT=CUT
    PCA.OUT$N=N
    #N=min(which( cumsum(PCA.OUT$sdev^2)/sum(PCA.OUT$sdev^2) > 0.7))
    return(PCA.OUT)
    }


vector.SeuratRandomPCA <-function(pbmc, RN=1000, CUT=0.5){
    pbmc=pbmc
    CUT=CUT
    RN=RN
    #####################
    library(gmodels)
    D=as.matrix(pbmc@assays$RNA@scale.data)
    R_INDEX=sample(1:ncol(D), RN, replace =TRUE)
    ###########################
    
    RD=D[,R_INDEX]
    PCA.OUT=fast.prcomp(t(RD), retx = TRUE, center = FALSE, scale. = FALSE, tol = NULL)
    EXP=(cumsum(PCA.OUT$sdev^2)/sum(PCA.OUT$sdev^2))
    N=min(which( cumsum(PCA.OUT$sdev^2)/sum(PCA.OUT$sdev^2) > CUT))
    
    PRED.PCA=t(D) %*% PCA.OUT$rotation
    
    #####################
    PCA.OUT$EXP=EXP
    PCA.OUT$CUT=CUT
    PCA.OUT$N=N
    PCA.OUT$RN=RN
    PCA.OUT$R_INDEX=R_INDEX
    PCA.OUT$PRED.PCA=PRED.PCA
    #N=min(which( cumsum(PCA.OUT$sdev^2)/sum(PCA.OUT$sdev^2) > 0.7))
    return(PCA.OUT)
    }




vector.medCurv <- function(PCA,MAX=1000){
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



# 2020.01.03
vector.smoothOut <- function(X, Z){
    Y=X
    Y[order(Z)]=X[order(Z)]-smooth(X[order(Z)])
    return(Y)
    }

vector.regressOut <- function(X, Z){
    FIT=lm(X~Z)
    Y=X-predict(FIT)    
    return(Y)
    }



vector.removeOut <- function(X){   
    X=X
    ################
    Q1=quantile(X,0.25)
    Q3=quantile(X,0.75)
    IQR=Q3-Q1
    LW=Q1-1.5*IQR
    UP=Q3+1.5*IQR
    ###############################
    X[which(X>UP)]=UP
    X[which(X<LW)]=LW
    ################    
    return(X)
    }


vector.rankPCA <- function(PCA){
    PCA=PCA #ncol=num(PC)
    R.PCA=apply(PCA,2,rank)
    library(gmodels)
    PCA.OUT=fast.prcomp(R.PCA, retx = TRUE, center = TRUE, scale. = TRUE, tol = NULL)
    N.PCA=PCA.OUT$x
    return(N.PCA)
    }



###################
#2020.1.12

.normNew <- function(x){
    pos_index=which(x>0)
    neg_index=which(x<0)
    x[pos_index]=rank(x[pos_index])
    x[neg_index]=rank(-x[neg_index])
    return(x)
    }



vector.calValueNew <- function(PCA){
    OUT=list()
    PCA=PCA
    PCA.RC=apply(PCA, 2, .normNew)   
    VALUE=apply(PCA.RC,1,mean)
    ###########################
    OUT$VALUE=VALUE
    OUT$PCA.RC=PCA.RC
    return(OUT)
    }



vector.getValueNew <-function(OUT, PCA, SHOW=TRUE){
    OUT=OUT
    PCA=PCA
    SHOW=SHOW
    ############
    VALUE.OUT=vector.calValueNew(PCA)
    ###############
    OUT$VALUE=VALUE.OUT$VALUE
    OUT$PCA=PCA
    OUT$PCA.RC=VALUE.OUT$PCA.RC
    ###############
    if(SHOW==TRUE){
        vector.showValue(OUT)
        }
    #####
    return(OUT)
    }

###############################################

#20201030





vector.autoCenterCor <- function(OUT, UP=0.9, SHOW=TRUE){
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
    
    #TMP=10^LENGTH - DIST_COR
    TMP=-DIST_COR
    #TMP=rank(-DIST_COR) * rank(LENGTH )
    SELECT=which(TMP==max(TMP))[1]
    #SELECT=which(LENGTH==max(LENGTH))[1]
    #SELECT=which( rank(-DIST_COR) * rank(LENGTH) == max( rank(-DIST_COR) * rank(LENGTH)  ) )[1]
    
    SUMMIT=CLUSTER[[SELECT]]
    #print(SELECT)
    #plot(OUT$CENTER_VEC[USED,])
    #points(OUT$CENTER_VEC[CLUSTER[[9]],],pch=16,col='red')
    
    
    SCORE=c()
    PS=c()
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
        PS=c(PS, this_score)
        i=i+1}
    SCORE=max(SCORE)-SCORE
    ##########################################
    
    VALUE=SCORE
    #plot(OUT$VEC, col='grey70',pch=16)
    N.VALUE=(VALUE-min(VALUE))/(max(VALUE)-min(VALUE))
    #COL=vector.vcol(N.VALUE, c(0,0.5,1),c('#009FFF','#FFF200','#ec2F4B')) # too strong
    #COL=vector.vcol(N.VALUE,c(0,0.5,1),c('#009FFF','#FFF200','#ffdde1')) # too weak
    COL=vector.vcol(N.VALUE,c(0,0.5,1),c('#009FFF','#FFF200','#ee9ca7'))
    
    ###################################################
    OUT$COL=rep('grey70',nrow(OUT$VEC))
    OUT$ORIG.COL=rep('grey70',nrow(OUT$VEC))
    OUT$P.SCORE=rep(0,nrow(OUT$VEC))
    OUT$P.PS=rep(NA,nrow(OUT$VEC))
    i=1
    while(i<=length(USED)){
        this_index=INDEX_LIST[[USED[i]]]
        OUT$COL[this_index]=COL[i]
        OUT$ORIG.COL[this_index]=OUT$ORIG.CENTER.COL[USED][i]
        OUT$P.SCORE[this_index]=SCORE[i]
        OUT$P.PS[this_index]=PS[i]
        i=i+1
        }
    ###########
    
    
    if(SHOW==TRUE){
        #plot(OUT$VEC, col=OUT$COL, pch=16, cex=0.5 )
        plot(OUT$VEC, col=OUT$ORIG.COL, pch=16, cex=0.5)
        text(CENTER_VEC[HIGH,1], CENTER_VEC[HIGH,2],labels=PCH,cex=1,pos=2)
        points(CENTER_VEC[HIGH,1], CENTER_VEC[HIGH,2], col='black',pch=16,cex=1)
        points(CENTER_VEC[SUMMIT,1],CENTER_VEC[SUMMIT,2], col='black',pch=16,cex=1.5)
        points(CENTER_VEC[SUMMIT,1],CENTER_VEC[SUMMIT,2], col='red',pch=16,cex=1)  
        }
 
    ######################
    ############
    OUT$SCORE=SCORE
    OUT$SUMMIT=SUMMIT
    OUT$CLUSTER=CLUSTER
    OUT$LENGTH=LENGTH
    OUT$PCH=PCH
    OUT$DIST_COR=DIST_COR
    OUT$PS=PS
    #OUT$DIST_MEAN=DIST_MEAN
    ################################
    
    #######################
    return(OUT)
    }




