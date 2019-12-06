
setwd('F:/Vector/data/HumanHepDiff_GSE81252')

D=read.csv('GSE81252_data.cast.log2.lineage.csv',header=T,row.names=1)

DATA=t(D[,2:ncol(D)])
TYPE=D[,1]

