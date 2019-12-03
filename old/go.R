
vector.getValue_old <-function(PCA){
    PCA=PCA
    r.pca=apply(abs(PCA), 2, rank)
    mean.r.pca=apply(r.pca,1,mean)
    return(mean.r.pca)
    }
