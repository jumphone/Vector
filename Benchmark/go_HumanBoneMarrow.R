#Citation: https://data.humancellatlas.org/explore/projects/cc95ff89-2e68-4a08-a234-480eca21ce79


#Python:

fi=open('F:/Vector/data/HumanBoneMarrow/5635f5a8-532e-483b-81a2-9c53495dbd90.csv/expression.csv')
i=1
this_line=fi.readline()
while this_line !='':
    this_line=fi.readline()
    i=i+1
print i
#
fi.close()

