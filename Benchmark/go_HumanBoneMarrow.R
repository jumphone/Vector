#Citation: https://data.humancellatlas.org/explore/projects/cc95ff89-2e68-4a08-a234-480eca21ce79



#Python:

fi=open('F:/Vector/data/HumanBoneMarrow/5635f5a8-532e-483b-81a2-9c53495dbd90.csv/expression.csv')
fo=open('F:/Vector/data/HumanBoneMarrow/5635f5a8-532e-483b-81a2-9c53495dbd90.csv/expression_30000.csv','w')
l1=fi.readline()
fo.write(l1)
i=1
while i<=30000:
    fo.write(fi.readline())
    print(i)
    i=i+1
fi.close()
fo.close()
#################################

#R:
setwd('F:/Vector/data/HumanBoneMarrow/5635f5a8-532e-483b-81a2-9c53495dbd90.csv/')
DATA=read.csv(file='expression_30000.csv',sep=',',header=TRUE,row.names=1)
DATA=t(DATA)

