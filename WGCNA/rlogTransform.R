library(DESeq2)

counts=read.table("AllCountsNH_Combined.tab",header=TRUE,row.names=1) 

# Remove AS sample
counts<-counts[,-c(1:32)]

traits=read.csv('colData_Combined.csv',row.names = 1)
traits$Type = as.factor(traits$Type)

type = traits$Type
type

conditions=data.frame(type)
conditions

########################### creating DESeq2 dataset - for WGCNA based pipeline

#Remove isogroups with low counts from dataset - count less than 2 in more than 90% of samples
counts$low = apply(counts[,1:76],1,function(x){sum(x<=2)})  #making new column counting number of samples with counts <=2 within each isogroup (host) 

counts2<-counts[-which(counts$low>68),] #68 is ~90% of 76 samples - get rid of count less than 2 in more than 90% of samples

# #now transform counts data using DESeq2
ddsCOUNTS<-DESeqDataSetFromMatrix(countData=counts2[,1:76],colData=conditions,design=~type)

# #rlog transform
rlogCOUNTS<-rlog(ddsCOUNTS,blind=TRUE) #use blind=TRUE to not account for experimental design


save(rlogCOUNTS, file="rlogCountsGreedy.RData") 