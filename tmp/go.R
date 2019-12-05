

vector.gridValueSmooth <- function(OUT, SHOW=TRUE){
    OUT=OUT
    SHOW=SHOW
    ####################################
    #INDEX_LIST=OUT$INDEX_LIST
    #VALUE=OUT$VALUE
    USED=OUT$USED
    CENTER_VALUE=OUT$CENTER_VALUE
    CENTER_VEC=OUT$CENTER_VEC
    p1=OUT$p1
    p2=OUT$p2
    p1_index=as.numeric(str_replace(p1,'P',''))
    p2_index=as.numeric(str_replace(p2,'P',''))
    library(igraph)
    #################
    DEG=degree(OUT$GRAPH,v = V(OUT$GRAPH))
    names(DEG)=as_ids(V(OUT$GRAPH))
    #W=DEG/max(DEG)
    #W=rank(DEG)/length(DEG)
    DEG=DEG[order(as.numeric(str_replace(names(DEG),'P','')))]
    
    #################
    NEW_CENTER_VALUE=CENTER_VALUE
    NB_VALUE_MEAN=c()
    NB_VALUE_MIN=c()
    NB_VALUE_MAX=c()
    NB_VALUE_SD=c()
    NB_DEG_MEAN=c()
    ###########################################
    
    i=1
    while(i<=length(CENTER_VALUE)){
        this_value=CENTER_VALUE[i]
        neighbor_p1_index=p1_index[which(p2_index==i)]
        neighbor_p2_index=p2_index[which(p1_index==i)]
        neighbor_index=unique(sort(c(neighbor_p1_index,neighbor_p2_index)))
        nb_value=CENTER_VALUE[neighbor_index]
        nb_value_sd=sd(nb_value)
        nb_value_mean=mean(c(nb_value,this_value))
        nb_value_min=min(c(nb_value,this_value))
        nb_value_max=max(c(nb_value,this_value))
        #this_value=CENTER_VALUE[i]
        #this_value= nb_value_mean + (this_value-nb_value_mean)
        NB_VALUE_MEAN=c(NB_VALUE_MEAN,nb_value_mean)
        NB_VALUE_MIN=c(NB_VALUE_MIN,nb_value_min)
        NB_VALUE_MAX=c(NB_VALUE_MAX,nb_value_max)
        NB_VALUE_SD=c(NB_VALUE_SD,nb_value_sd) 
        NB_DEG_MEAN=c(NB_DEG_MEAN, mean(DEG[c(neighbor_index,i)]) )
        i=i+1}
    
    #RATIO = NB_VALUE_SD / max(NB_VALUE_SD)
    #W = NB_DEG_MEAN/max(NB_DEG_MEAN)
    #V =  NB_VALUE_MEAN / max(NB_VALUE_MEAN)
    
    RATIO=1#rank(NB_VALUE_SD)
    W=rank(NB_DEG_MEAN)
    V=rank(CENTER_VALUE)
    
    
    NEW_CENTER_VALUE =  V *  W  * RATIO# + NB_VALUE_MIN  #* W
    
    
    #####################################
    N.VALUE=.normX(NEW_CENTER_VALUE)#(VALUE-min(VALUE))/(max(VALUE)-min(VALUE))
    ORIG.CENTER.COL=vector.vcol(N.VALUE, c(0,0.5,1),c('#009FFF','#FFF200','#ec2F4B')) 
    #######################
    if(SHOW==TRUE){
            plot(OUT$VEC[,1],OUT$VEC[,2],col='grey80',cex=0.5, pch=16)
            points(CENTER_VEC[USED,1],CENTER_VEC[USED,2], col=ORIG.CENTER.COL[USED], pch=15, cex=1.5)
            }
    #################
    OUT$CENTER_VALUE=NEW_CENTER_VALUE       
    OUT$ORIG.CENTER.COL=ORIG.CENTER.COL
    return(OUT)
    }



vector.gridValueSmooth <- function(OUT,CUT=0.95, SHOW=TRUE, MAX=2000){
    OUT=OUT
    SHOW=SHOW
    CUT=CUT
    MAX=MAX
    INDEX_LIST=OUT$INDEX_LIST
    VALUE=OUT$VALUE
    USED=OUT$USED
    ################
    CENTER_VALUE=OUT$CENTER_VALUE
    CENTER_VEC=OUT$CENTER_VEC
    p1=OUT$p1
    p2=OUT$p2
    
    
    vector.smooth <- function(CUT, MAX, CENTER_VALUE, p1, p2){
    
        NEW_CENTER_VALUE=CENTER_VALUE
        p1_value=c()
        p2_value=c()
        DIFF=c()
        i=1
        while(i<=length(p1)){
            this_p1=p1[i]
            this_p2=p2[i]
            this_p1_index=as.numeric(str_replace(this_p1,'P',''))
            this_p2_index=as.numeric(str_replace(this_p2,'P',''))
            this_p1_value=NEW_CENTER_VALUE[this_p1_index]
            this_p2_value=NEW_CENTER_VALUE[this_p2_index]
            p1_value=c(p1_value, this_p1_value)
            p2_value=c(p2_value, this_p2_value)
            this_diff=this_p1_value-this_p2_value
            DIFF=c(DIFF, this_diff)
            i=i+1}    
        ABS_DIFF=abs(DIFF)
        POS_NUM=length(which(ABS_DIFF>0))
        ##############################################
        ABS_DIFF_COR=cor(1:length(ABS_DIFF),sort(ABS_DIFF))
        ############################
        COR_HIST=c()
        TIME=1
        while(ABS_DIFF_COR < CUT & TIME <= MAX){
            ################
            target_abs_diff= median(ABS_DIFF)
        
            ###############################
        
            this_max_index=which(ABS_DIFF==max(ABS_DIFF))[1]

            max_p1_index=as.numeric(str_replace( p1[this_max_index] ,'P',''))
            max_p2_index=as.numeric(str_replace( p2[this_max_index] ,'P',''))
            
            ###########################
            #NEW_CENTER_VALUE=NEW_CENTER_VALUE-target_abs_diff/2
            #NEW_CENTER_VALUE[which(NEW_CENTER_VALUE<0)]=0
            #################################
            if(NEW_CENTER_VALUE[max_p1_index] < NEW_CENTER_VALUE[max_p2_index]){
            
                NEW_CENTER_VALUE[max_p2_index]=NEW_CENTER_VALUE[max_p1_index] + target_abs_diff #/2          
                #NEW_CENTER_VALUE[max_p1_index] = max(0, NEW_CENTER_VALUE[max_p1_index]- target_abs_diff/2)
                ##########################################
                p1_changed_index=which( p1==p2[this_max_index] )
                ABS_DIFF[p1_changed_index]=abs( p1_value[p1_changed_index]- target_abs_diff - p2_value[p1_changed_index] )
                p2_changed_index=which( p2==p2[this_max_index] )
                ABS_DIFF[p2_changed_index]=abs( p1_value[p2_changed_index]+ target_abs_diff - p2_value[p2_changed_index] )
                ######################################
                }else{
            
                NEW_CENTER_VALUE[max_p1_index]=NEW_CENTER_VALUE[max_p2_index] + target_abs_diff #/2
                #NEW_CENTER_VALUE[max_p2_index]= max(0, NEW_CENTER_VALUE[max_p2_index]- target_abs_diff/2)
                ##########################################
                p1_changed_index=which( p1==p1[this_max_index] )
                ABS_DIFF[p1_changed_index]=abs( p1_value[p1_changed_index]+ target_abs_diff - p2_value[p1_changed_index] )
                
                p2_changed_index=which( p2==p1[this_max_index] )
                ABS_DIFF[p2_changed_index]=abs( p1_value[p2_changed_index]- target_abs_diff - p2_value[p2_changed_index] )
                ######################################
                }
        
            #ABS_DIFF[this_max_index]= target_abs_diff
            
            ##############################################
         
            ###############################
            ABS_DIFF_COR=cor(1:length(ABS_DIFF),sort(ABS_DIFF))
        
            ############################
            COR_HIST=c(COR_HIST, ABS_DIFF_COR)
            POS_NUM=length(which(ABS_DIFF>0))
            TIME=TIME+1
        }
        ###################################################
    
        NEW_CENTER_VALUE[which(!c(1:nrow(OUT$CENTER_VEC)) %in% USED)]=0
        NEW_CENTER_VALUE[USED]=.normX(NEW_CENTER_VALUE[USED])
        
        N.VALUE=.normX(NEW_CENTER_VALUE)#(VALUE-min(VALUE))/(max(VALUE)-min(VALUE))
        ORIG.CENTER.COL=vector.vcol(N.VALUE, c(0,0.5,1),c('#009FFF','#FFF200','#ec2F4B')) 
        #####################
        OUT=list()
        OUT$NEW_CENTER_VALUE=NEW_CENTER_VALUE
        OUT$TIME=TIME
        OUT$COR_HIST=COR_HIST
        OUT$ABS_DIFF_COR=ABS_DIFF_COR
        OUT$ORIG.CENTER.COL=ORIG.CENTER.COL
        return(OUT)
        }        
        
    SMOOTH.OUT= vector.smooth(CUT, MAX, CENTER_VALUE, p1, p2)
    ABS_DIFF_COR=SMOOTH.OUT$ABS_DIFF_COR
    NEW_CENTER_VALUE=SMOOTH.OUT$NEW_CENTER_VALUE
    ORIG.CENTER.COL=SMOOTH.OUT$ORIG.CENTER.COL
    COR_HIST=SMOOTH.OUT$COR_HIST
    TIME=SMOOTH.OUT$TIME
    ###############
        
    ####################################
    if(ABS_DIFF_COR<CUT){
        
        print('CUT is too high!!!')
        print(paste0('Max CUT should be less than: ', max(COR_HIST) ) )
        print(paste0('CUT is changed to: ',  max(COR_HIST) ) )
        ##########################
        SMOOTH.OUT= vector.smooth(max(COR_HIST), MAX, CENTER_VALUE, p1, p2)
        ABS_DIFF_COR=SMOOTH.OUT$ABS_DIFF_COR
        NEW_CENTER_VALUE=SMOOTH.OUT$NEW_CENTER_VALUE
        ORIG.CENTER.COL=SMOOTH.OUT$ORIG.CENTER.COL
        COR_HIST=SMOOTH.OUT$COR_HIST
        TIME=SMOOTH.OUT$TIME
        
        }
    
    if(SHOW==TRUE){
            plot(OUT$VEC[,1],OUT$VEC[,2],col='grey80',cex=0.5, pch=16)
            points(CENTER_VEC[USED,1],CENTER_VEC[USED,2], col=ORIG.CENTER.COL[USED], pch=15, cex=1.5)
            }   
    OUT$CENTER_VALUE=NEW_CENTER_VALUE
    OUT$ORIG.CENTER.COL=ORIG.CENTER.COL
    OUT$COR_HIST=COR_HIST
    return(OUT)
    }



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




