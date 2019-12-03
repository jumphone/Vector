###################
# Author: Feng Zhang
# Date: 11/25, 2019

##################
library('circlize')
library('gatepoints')
library('stringr')
library('igraph')
##################


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



vector.getValue <-function(OUT, PCA){
    OUT=OUT
    PCA=PCA
    PCA.RC=apply(apply(PCA,2,rank), 2, .normX)
    PCA.RC=abs(PCA.RC-0.5)   
    VALUE=apply(PCA.RC,1,mean)
    #####
    OUT$VALUE=VALUE
    OUT$PCA=PCA
    OUT$PCA.RC=PCA.RC
    return(OUT)
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




vector.buildGrid <- function(VEC, N=20, SHOW=TRUE, COL='grey70'){
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
    while(i<CNUM){
    
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
    ##########################
    
    library(igraph)
    OUT$p1=p1
    OUT$p2=p2
    NET = cbind(p1,p2) 
    g <- make_graph(t(NET),directed = FALSE)
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
    ###################
    OUT$GRAPH=g
    OUT$DIST=DIST
    OUT$USED=USED
    OUT$USED_NAME=USED_NAME
    return(OUT)
    ###########################
    }






vector.gridValue <- function(OUT, VALUE， SHOW=TRUE){
    OUT=OUT
    INDEX_LIST=OUT$INDEX_LIST
    VALUE=VALUE
    SHOW=SHOW
    USED=OUT$USED
    ################
    CENTER_VALUE=c()
    i=1
    while(i<=length(INDEX_LIST)){
        CENTER_VALUE=c(CENTER_VALUE,mean(VALUE[INDEX_LIST[[i]]]))
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
        
        plot(OUT$VEC,col='grey70',pch=16,cex=0.5)
        points(OUT$CENTER_VEC,col=COL,pch=16)
        points(OUT$CENTER_VEC[USED,],col='black',pch=0, cex=1.5)
        }
    ################    
    OUT$CENTER_VALUE=CENTER_VALUE     
    OUT$ORIG.CENTER.COL=COL
    OUT$VALUE=VALUE
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
    i=1
    while(i<=SUB_CPT$no){
        this_name=names(which(SUB_CPT$membership==i))
        this_index=as.numeric(str_replace(this_name,'P',''))  
        PCH[which(HIGH %in% this_index)]=as.character(i)
        LENGTH=c(LENGTH, length(this_index))
        CLUSTER=c(CLUSTER,list(this_index))       
        i=i+1}
   
    ####################
    SELECT=which(LENGTH==max(LENGTH))[1]
    SUMMIT=CLUSTER[[SELECT]]
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
        plot(OUT$VEC, col=OUT$ORIG.COL, pch=16, cex=0.5)
        text(CENTER_VEC[HIGH,],labels=PCH,cex=1,pos=2)
        points(CENTER_VEC[HIGH,], col='black',pch=16,cex=1)
        points(CENTER_VEC[SUMMIT,], col='black',pch=16,cex=1.5)
        points(CENTER_VEC[SUMMIT,], col='red',pch=16,cex=1)  
        }
    
    ######################
    ############
    OUT$SCORE=SCORE
    OUT$SUMMIT=SUMMIT
    OUT$CLUSTER=CLUSTER
    OUT$LENGTH=LENGTH
    OUT$PCH=PCH
    ################################
    
    #######################
    return(OUT)
    }




vector.drawArrow <- function(OUT, P=0.9, SHOW=TRUE, COL='grey70'){
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
    

    ###################
    .norm_one <-function(x){
        if(var(x)!=0){
            x=x/sqrt(sum(x^2))}
        return(x)
        }
    
    DIV=1/P
    if(SHOW==TRUE){
        plot(ALL_VEC,col=COL,pch=16,cex=0.2)
        }
    
    A1_VEC=c()
    A2_VEC=c()
    A_LENGTH=c()
    i=1
    while(i<=length(USED)){
        this_p1_loc=USED_CENTER_VEC[i,]
        
        vector_list=cbind(USED_CENTER_VEC[,1]-this_p1_loc[1],USED_CENTER_VEC[,2]-this_p1_loc[2])
        vector_list_norm=t(apply(vector_list,1,.norm_one))
        
        vector_weight_1= DIV^-(rank(USED_DIST[i,])-1)  
        vector_weight_2= SCORE[i]-SCORE
        
        vector_weight = vector_weight_1 * vector_weight_2        
        vector_weight = vector_weight/sum(abs(vector_weight))
        
        final_vec=t(vector_list_norm) %*% vector_weight
        
        this_p2_loc=c(this_p1_loc[1]+final_vec[1],this_p1_loc[2]+final_vec[2])
        this_arrow_length=0.1*sqrt(sum(final_vec^2))
        
        if(SHOW==TRUE){
            arrows(x0=this_p1_loc[1],y0=this_p1_loc[2],
                   x1=this_p2_loc[1],y1=this_p2_loc[2],
                   lwd=2, length=this_arrow_length,
                   col='black'
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
    ###########
    return(OUT)
    }


VECTOR <- function(VEC, PCA, N=20){
    N=N
    par(mfrow=c(3,2))
    ########################
    
    OUT=vector.buildGrid(VEC, N=N,SHOW=TRUE)
    OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
    VALUE=vector.getValue(PCA)
    OUT=vector.gridValue(OUT,VALUE, SHOW=TRUE)
    OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
    OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)

    par(mfrow=c(1,1))
    return(OUT)
    #####
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
        plot(OUT$VEC, col=OUT$COL, pch=16, cex=0.5 )
        points(CENTER_VEC[SUMMIT,], col='black',pch=16,cex=1.5)
        points(CENTER_VEC[SUMMIT,], col='red',pch=16,cex=1)  
        }
    
    ######################
    ############
    OUT$SCORE=SCORE
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
    #####################
    CENTER_LIST=OUT$CENTER_LIST
    INDEX_LIST=OUT$INDEX_LIST
    USED=OUT$USED  
    ################################
    SELECT_NAME=vector.selectPoint(VEC,CEX=0.1)
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
    ########################
    return(OUT)
    }







