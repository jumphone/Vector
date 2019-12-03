
vector.getValue_old <-function(PCA){
    PCA=PCA
    r.pca=apply(abs(PCA), 2, rank)
    mean.r.pca=apply(r.pca,1,mean)
    return(mean.r.pca)
    }


vector.autoCenterNew <- function(OUT,  SHOW=TRUE){
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
    .tmp<-function(x){
        return(cor(x,CENTER_VALUE[USED],method='spearman'))
        }
    DIST_USED.COR=apply(DIST[USED,USED],2,.tmp)
    SUMMIT=USED[which(DIST_USED.COR==min(DIST_USED.COR))]
    SUMMIT_NAME=paste0('P',SUMMIT)
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
    ######################
    ############
    OUT$SCORE=SCORE
    OUT$SUMMIT=SUMMIT
    ################################
    if(SHOW==TRUE){
        #plot(OUT$VEC, col=OUT$COL, pch=16, cex=0.5 )
        plot(OUT$VEC, col=OUT$ORIG.COL, pch=16, cex=0.5)
        #text(CENTER_VEC[HIGH,],labels=PCH,cex=1,pos=2)
        #points(CENTER_VEC[HIGH,], col='black',pch=16,cex=1)
        points(CENTER_VEC[SUMMIT,1], CENTER_VEC[SUMMIT,2], col='black',pch=16,cex=1.5)
        points(CENTER_VEC[SUMMIT,1], CENTER_VEC[SUMMIT,2], col='red',pch=16,cex=1)  
        }
    #######################
    return(OUT)


     }





vector.calScore <-function(OUT,PCA,SHOW=TRUE){
    OUT=OUT
    INDEX_LIST=OUT$INDEX_LIST
    #PCA=t(DATA[which(rownames(DATA) %in% VariableFeatures(pbmc)),])
    PCA=PCA
    SHOW=SHOW
    DIST=OUT$DIST
    CENTER_VEC=OUT$CENTER_VEC
    #CENTER_PCA=c()
    USED=OUT$USED
    USED_NAME=OUT$USED_NAME
    ################
    CENTER_PCA=c()
    i=1
    while(i<=length(INDEX_LIST)){
        this_index=INDEX_LIST[[i]]
        if(length(this_index)==1){
            CENTER_PCA=cbind(CENTER_PCA, PCA[this_index,])
            }else{
            CENTER_PCA=cbind(CENTER_PCA, apply(PCA[this_index,],2,mean))
            }
        #print(i)
        i=i+1}
    ######################################
    CENTER_PCA=t(CENTER_PCA)
    rownames(CENTER_PCA)=colnames(DIST)
    colnames(CENTER_PCA)=colnames(PCA)
    ####################################

    used_p_index=which(OUT$p1 %in% USED_NAME & OUT$p2 %in% USED_NAME)
    used_p1=OUT$p1[used_p_index]
    used_p2=OUT$p2[used_p_index]
    
    CENTER_PCA.N=apply(CENTER_PCA,2,.normX)
    PCA.N=apply(PCA,2,rank)
    
    LENGTH=c()
    SSS=c()
    i=1
    while(i<=length(USED)){
        #this_pca_dist=#CENTER_PCA_SCORE[USED]-CENTER_PCA_SCORE[USED[i]]#CENTER_PCA_DIST[USED[i],USED]
        #CENTER_PCA[USED[i],]
        
        this_net_dist=DIST[USED[i],USED]
        
        this_net_dist_cluster=round(.normX(rank(this_net_dist)),1)
        step_list=sort(unique(this_net_dist_cluster))
        
        step_score=c()
        used_list=c()
        j=1
        while(j<=length(step_list)){
            this_step=step_list[j]
            this_index=which(this_net_dist_cluster==this_step)
            this_orig_index=c()
            t=1
            while(t<=length(this_index)){
                this_orig_index=c(this_orig_index, INDEX_LIST[[USED[this_index[t]]]])
                t=t+1}
               
            if(length(this_orig_index)>1){
                this_tmp=abs(cor(PCA[this_orig_index,],method='spearman'))
                this_step_score = mean(unique(this_tmp))
                step_score=c(step_score,this_step_score)
                used_list=c(used_list, this_step)
                }
            j=j+1}
        
        this_score=cor(used_list, step_score,method='spearman')#mean(this_mut_list)
        SSS=c(SSS,this_score)
        
        if(i %%10==1){print(paste0(i,'/',length(USED)))}
        i=i+1}
        
        
    #plot(SSS)
  
    #################################
    VALUE=SSS#TMP
    N.VALUE=(VALUE-min(VALUE))/(max(VALUE)-min(VALUE))
    COL=vector.vcol(N.VALUE, c(0,0.5,1),c('#009FFF','#FFF200','#ec2F4B'))
    
    
    SSS.N=.normX(SSS)
    
    SCORE=rep(0,nrow(OUT$VEC))
    i=1
    while(i<=length(USED)){
        this_index=INDEX_LIST[[USED[i]]]
        SCORE[this_index]=SSS.N[i]
        if(is.na(SSS.N[USED[i]])){print(i)}
        #SUMMIT_NEG_DIST[USED[i]]
        i=i+1
        }
    
    names(SCORE)=rownames(OUT$VEC)
    
    VALUE=SCORE
    N.VALUE=(VALUE-min(VALUE))/(max(VALUE)-min(VALUE))
    COL=vector.vcol(N.VALUE, c(0,0.5,1),c('#009FFF','#FFF200','#ec2F4B'))
    ################
    if(SHOW==TRUE){
        
        plot(OUT$VEC,col=COL,pch=16,cex=0.5)

        }
   
    return(SCORE)

    }




