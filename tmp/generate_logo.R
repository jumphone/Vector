

vector.selectPoint<-function(VEC,CEX=0.5){
    VEC=VEC
    CEX=CEX
    ############
    library(gatepoints)
    selectedPoints <- fhs(VEC, pch=16,col='red3',cex=CEX,mark = TRUE)
    return(selectedPoints)
    }



set.seed(123)
X=runif(10000)
Y=runif(10000)

VEC=cbind(X,Y)
rownames(VEC)=paste0(1:nrow(VEC))
plot(X,Y,pch=16,cex=0.5,col='grey60')
OUT=vector.selectPoint(VEC)

plot(VEC[as.numeric(OUT),])

VEC=VEC[as.numeric(OUT),]
saveRDS(VEC, file='F:/Vector/data/logoVEC.RDS')
##################################


VEC=readRDS( file='F:/Vector/data/logoVEC.RDS')

#VEC[,1]=scale(VEC[,1])
#VEC[,2]=scale(VEC[,2])




    OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
    OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
    #OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
    OUT$VALUE=abs(rnorm(nrow(VEC)))
    
    OUT=vector.gridValue(OUT,SHOW=TRUE)
    OUT=vector.selectCenter(OUT)
    #OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
    OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL)   
    OUT=vector.reDrawArrow(OUT,SHOW=TRUE) 

    #tiff(paste0("F:/Vector/data/TMP.tiff"),width=4,height=4,units='in',res=600)
    #par(mar=c(0,0,0,0))
    
    #OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL,AL=50)
    #dev.off()


    
    



