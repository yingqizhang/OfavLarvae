library(WGCNA)
library(flashClust)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(diagram)
library(ggplot2)
library(cowplot)

# For PCA analysis
#install.packages(c("FactoMineR", "factoextra"))
library("FactoMineR")
library("factoextra")

CombCts=read.csv('AllCountsNH_Combined.csv',row.names=1,check.names=FALSE)
head(CombCts)
ncol(CombCts) #108 samples in total
nrow(CombCts) #69109 isogroups in total

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

##############Normalize reads#################
# Ran rlogTransform.R on CARC

load(file = "rlogCountsGreedy.RData")
head(assay(rlogCOUNTS))

dat=as.data.frame(assay(rlogCOUNTS))

nrow(dat) # 28675 genes after removing isogroups with count less than 2 in more than 90% of samples; started with 69109 genes

# Reorganize data so that row corresponds to sample and column corresponds to gene
datExpr0 = as.data.frame(t(dat[,1:76]))

# Check for genes outliers and samples with too many missing values
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK # TRUE

#Read in trait data
datTraits= read.csv('Trait.csv',row.names = 1)
dim(datTraits)
names(datTraits)

table(rownames(datTraits[1:46,])==rownames(datExpr0[1:46,])) # first 46 samples have matching names

rownames(datExpr0) = rownames(datTraits)

table(rownames(datTraits)==rownames(datExpr0))

# Cluster samples by expression
A=adjacency(t(datExpr0),type="signed") # SELECT SIGNED OR UNSIGNED HERE # this calculates the whole network connectivity
k = as.numeric(apply(A,2,sum))-1 # standardized connectivity
Z.k = scale(k)
thresholdZ.k = -2.5 # often -2.5
outlierColor = ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average")

quartz()

par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

# Convert traits to a color representation: red means high
traitColors = data.frame(numbers2colors(datTraits,signed=FALSE))
dimnames(traitColors)[[2]] = paste(names(datTraits))
datColors = data.frame(outlier = outlierColor,traitColors) # Combine outlier color info with trait info
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree,groupLabels=names(datColors),colors=datColors,main="Sample Dendrogram and Trait Heatmap")

# Remove outlying samples 
remove.samples = Z.k<thresholdZ.k | is.na(Z.k)
datExprOut = datExpr0[!remove.samples,]
datTraitsOut = datTraits[!remove.samples,]

datExprOut2019= datExprOut[1:45,]
datExprOut2021= datExprOut[46:74,]
datTraitsOut2019 = datTraitsOut[1:45,-c(4:5)]
datTraitsOut2021 = datTraitsOut[46:74,-c(1,4:5)]

save(datExprOut, datTraitsOut, datExprOut2019, datExprOut2021, datTraitsOut2019, datTraitsOut2021, file="SamplesAndTraits_OutliersRemoved.RData")

#######PCA analysis########

write.table(datTraitsOut2019,"datTraitsOut2019_PCA.csv", row.names=F, sep=",")
write.table(datTraitsOut2021,"datTraitsOut2021_PCA.csv", row.names=F, sep=",")

datTraitsOut2019_PCA = read.csv("datTraitsOut2019_PCA.csv")
rownames(datTraitsOut2019_PCA) = rownames(datTraitsOut2019)
datTraitsOut2021_PCA = read.csv("datTraitsOut2021_PCA.csv")
rownames(datTraitsOut2021_PCA) = rownames(datTraitsOut2021)

lar2019.pca = PCA(datExprOut2019, graph = FALSE)
lar2021.pca = PCA(datExprOut2021, graph = FALSE)

ind <- get_pca_ind(lar2019.pca)
#PCcoord2019 = data.frame(ind$coord[,1:2])
PCcoord2019 = data.frame(ind$coord[,2:3])
#names(PCcoord2019) = c("PC1", "PC2")
names(PCcoord2019) = c("PC2", "PC3")
table(row.names(PCcoord2019) == row.names(datTraitsOut2019_PCA)) # TURE
PCcoord2019 = cbind(PCcoord2019, datTraitsOut2019_PCA)
PCcoord2019$Type = factor(PCcoord2019$Type, levels = c("Inshore", "Offshore", "Hybrid"))
PCcoord2019$Trmt = as.factor(PCcoord2019$Trmt)
#names(PCcoord2019) = c("PC1", "PC2", "Origin", "Trmt")
names(PCcoord2019) = c("PC2", "PC3", "Origin", "Trmt")

ind <- get_pca_ind(lar2021.pca)
PCcoord2021 = data.frame(ind$coord[,1:2])
names(PCcoord2021) = c("PC1", "PC2")
table(row.names(PCcoord2021) == row.names(datTraitsOut2021_PCA)) # TURE
PCcoord2021 = cbind(PCcoord2021, datTraitsOut2021_PCA[, 2:4])
PCcoord2021$Type = factor(PCcoord2021$Type, levels = c("Offshore", "Cross1", "Cross2"))
PCcoord2021$Trmt = as.factor(PCcoord2021$Trmt)
names(PCcoord2021) = c("PC1", "PC2", "Origin0", "Origin", "Trmt")

quartz()
g1=ggplot(PCcoord2019, aes(x=PC1, y=PC2,colour=Origin,shape=Trmt)) +
  geom_point(aes(colour=Origin,shape=Trmt),size=3)+
  scale_color_manual(values=c("#FC4E07","#00AFBB","#E7B800"),
                     labels=c("CR x CR", "HR x HR", "Hybrid"))+
  scale_shape_manual(values = c(1,16),
                     labels=c("Ctrl", "Heat"))+
  xlab(paste0("PC1: ",round(lar2019.pca$eig[1,2], 1),"% variance"))+
  ylab(paste0("PC2: ",round(lar2019.pca$eig[2,2], 1),"% variance"))+
  #theme_bw()+
  guides(shape = guide_legend(order = 2),col = guide_legend(order = 1))+
  ggtitle("2019")

# plotting PC2 vs. PC3
ggplot(PCcoord2019, aes(x=PC2, y=PC3,colour=Origin,shape=Trmt)) +
  geom_point(aes(colour=Origin,shape=Trmt),size=3)+
  scale_color_manual(values=c("#FC4E07","#00AFBB","#E7B800"),
                     labels=c("CR x CR", "HR x HR", "Hybrid"))+
  scale_shape_manual(values = c(1,16),
                     labels=c("Ctrl", "Heat"))+
  xlab(paste0("PC2: ",round(lar2019.pca$eig[2,2], 1),"% variance"))+
  ylab(paste0("PC3: ",round(lar2019.pca$eig[3,2], 1),"% variance"))+
  #theme_bw()+
  guides(shape = guide_legend(order = 2),col = guide_legend(order = 1))

g2=ggplot(PCcoord2021, aes(x=PC1, y=PC2,colour=Origin,shape=Trmt)) +
  geom_point(aes(colour=Origin,shape=Trmt),size=3)+
  scale_color_manual(values=c("#00AFBB","#FCD2B2","#FF8D33"),
                     labels=c("HR x HR", "CR x HR", "HR x CR"))+
  scale_shape_manual(values = c(1,16),
                     labels=c("Ctrl", "Heat"))+
  xlab(paste0("PC1: ",round(lar2021.pca$eig[1,2], 1),"% variance"))+
  ylab(paste0("PC2: ",round(lar2021.pca$eig[2,2], 1),"% variance"))+
  #theme_bw()+
  guides(shape = guide_legend(order = 2),col = guide_legend(order = 1))+
  ggtitle("2021")

quartz()
plot_grid(g1, g2, labels = c("a", "b"))

# Choose a soft threshold power
powers = c(seq(1,20,by=1)) #may need to adjust these power values to hone in on proper sft value
sft = pickSoftThreshold(datExprOut, powerVector=powers, verbose =5, networkType="signed") #call network topology analysis function

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, col="red")
abline(h=0.90, col="red")
# Choose power of 4

############Correlating general network properties#############
softPower = 4
rankExp2019= rank(colMeans(datExprOut[1:45,]))
rankExp2021= rank(colMeans(datExprOut[46:74,]))
random5000= sample(colnames(datExprOut),5000)
rankConn2019= rank(softConnectivity(datExprOut[1:45, random5000],type="signed",power=softPower)) 
rankConn2021= rank(softConnectivity(datExprOut[46:74, random5000],type="signed",power=softPower))

pdf("generalNetworkProperties.pdf", height=10, width=9)
par(mfrow=c(1,2))
verboseScatterplot(rankExp2019,rankExp2021, xlab="Ranked Expression (2019)",
                   ylab="Ranked Expression (2021)") 
verboseScatterplot(rankConn2019,rankConn2021, xlab="Ranked Connectivity (2019)",
                   ylab="Ranked Connectivity (2021)")

############Run WGCNA on the data sets#############
# Ran TOMstep.R on /scratch1/yingqizh/WGCNA2022
# Had to split up the job due to high computational demand

# Takes a min to load...
load(file = "TOM_Output2019.RData")
load(file = "TOM_Output2021.RData")

pdf("dendrogram.pdf",height=6,width=16)
par(mfrow=c(1,2))
plot(geneTree2019,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (2019)", labels=FALSE,hang=0.04)
plot(geneTree2021,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (2021)", labels=FALSE,hang=0.04)
dev.off()

mColorh=NULL
for (ds in 0:3){
  tree = cutreeHybrid(dendro = geneTree2019, pamStage=FALSE, minClusterSize = (30-3*ds), cutHeight = 0.99, deepSplit = ds, distM = dissTOM2019)
  mColorh=cbind(mColorh,labels2colors(tree$labels)); }
pdf("Module_choices.pdf", height=10,width=25);
plotDendroAndColors(geneTree2019, mColorh, paste("dpSplt =",0:3), main = "",dendroLabels=FALSE); dev.off()

save(mColorh, file="mColorh.RData")
modules2019 = mColorh[,1] # (Chosen deepsplit = 0, based on pdf plot generated)

PCs1_2019= moduleEigengenes(datExprOut2019, colors=modules2019)
ME1_2019= PCs1_2019$eigengenes
distPC1_2019= 1-abs(cor(ME1_2019,use="p"))
distPC1_2019= ifelse(is.na(distPC1_2019), 0, distPC1_2019)
pcTree1_2019= hclust(as.dist(distPC1_2019),method="a")
MDS1_2019 = cmdscale(as.dist(distPC1_2019),2)
colors2019 = names(table(modules2019))

quartz()
plot(pcTree1_2019, xlab="",ylab="",main="",sub="")
plot(MDS1_2019, col= colors2019, main="MDS plot", cex=2, pch=19)

# Qualitatively and quantitatively measure network preservation at the module level.

pdf("Final_modules.pdf",height=8,width=12)
plotDendroAndColors(geneTree2019, modules2019, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE,
                    guideHang=0.05, main="Gene dendrogram and module colors (2019)") 
plotDendroAndColors(geneTree2021, modules2019, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE,
                    guideHang=0.05, main="Gene dendrogram and module colors (2021)") 
dev.off()

multiExpr = list(A1=list(data=datExprOut2019),A2=list(data=datExprOut2021))
multiColor = list(A1 = modules2019) 
mp=modulePreservation(multiExpr,multiColor,referenceNetworks=1,verbose=3,networkType="signed",
                      nPermutations=30,maxGoldModuleSize=100,maxModuleSize=400) 
stats = mp$preservation$Z$ref.A1$inColumnsAlsoPresentIn.A2 
stats[order(-stats[,2]),c(1:2)]

# 5<Z<10 indicates moderate preservation, while Z>10 indicates high preservation.

# moduleSize Zsummary.pres
# black              400     39.160347
# tan                196     22.980468
# cyan               152     19.996518
# brown              400     19.568553
# green              400     15.428211
# turquoise          400     13.408618
# blue               400     10.448824
# magenta            283      8.037688
# greenyellow        268      8.016037
# pink               400      6.768167
# yellow             400      6.252689
# red                400      6.002128
# salmon             160      5.780293
# purple             268      5.109757
# gold               100      2.809810
# grey               400      2.652799

# Module membership (kME) and its use in comparing networks

geneModuleMembership1 = signedKME(datExprOut2019, ME1_2019)
colnames(geneModuleMembership1)=paste("PC",colors2019,".cor",sep="")
MMPvalue1=corPvalueStudent(as.matrix(geneModuleMembership1),dim(t(datExprOut2019))[[2]])
colnames(MMPvalue1)=paste("PC",colors2019,".pval",sep="")
Gene = rownames(t(datExprOut2019)) 
kMEtable1 = cbind(Gene,Gene,modules2019)

for (i in 1:length(colors2019)){
  kMEtable1 = cbind(kMEtable1, geneModuleMembership1[,i], MMPvalue1[,i]) 
}
colnames(kMEtable1)=c("PSID","Gene","Module",sort(c(colnames(geneModuleMembership1), colnames(MMPvalue1))))
write.csv(kMEtable1,"kMEtable1.csv",row.names=FALSE)

# Now repeat for 2021, using the module assignments from 2019 to determine kME values.

# First calculate MEs for 2021, since we haven't done that yet
PCs2_2021 = moduleEigengenes(datExprOut2021, colors=modules2019)
ME2_2021 = PCs2_2021$eigengenes

geneModuleMembership2 = signedKME(datExprOut2021, ME2_2021)
colnames(geneModuleMembership2)=paste("PC",colors2019,".cor",sep="")
MMPvalue2=corPvalueStudent(as.matrix(geneModuleMembership2),dim(t(datExprOut2021))[[2]])
colnames(MMPvalue2)=paste("PC",colors2019,".pval",sep="")
kMEtable2 = cbind(Gene,Gene,modules2019)

for (i in 1:length(colors2019)){
  kMEtable2 = cbind(kMEtable2, geneModuleMembership2[,i], MMPvalue2[,i]) 
}
colnames(kMEtable2)=colnames(kMEtable1)

write.csv(kMEtable2,"kMEtable2.csv",row.names=FALSE)

# plot the kME values of each gene in 2019 against the corresponding kME values of each gene in 2021
pdf("all_kMEtable2_vs_kMEtable1.pdf",height=8,width=8) 
for (c in 1:length(colors2019)){
  verboseScatterplot(geneModuleMembership2[,c],geneModuleMembership1[,c],main=colors2019[c], xlab="kME in 2021",ylab="kME in 2019")
}; dev.off()

pdf("inModule_kMEtable2_vs_kMEtable1.pdf",height=8,width=8) 
for (c in 1:length(colors2019)){
  inMod = modules2019== colors2019[c] 
  verboseScatterplot(geneModuleMembership2[inMod,c],geneModuleMembership1[inMod,c],main=colors2019[c], xlab="kME in 2021",ylab="kME in 2019")
}; dev.off()  

# determine which genes have extremely high kME values in both networks
topGenesKME = NULL
for (c in 1:length(colors2019)){
  kMErank1 = rank(-geneModuleMembership1[,c])
  kMErank2 = rank(-geneModuleMembership2[,c])
  maxKMErank = rank(apply(cbind(kMErank1,kMErank2+.00001),1,max)) 
  topGenesKME = cbind(topGenesKME,Gene[maxKMErank<=10])
}; colnames(topGenesKME) = colors2019 
topGenesKME

# black                           blue                           
# [1,] "isogroupOFAVBQ_DN187435_c0_g1" "isogroupOFAVBQ_DN108629_c0_g1"
# [2,] "isogroupOFAVBQ_DN194548_c0_g1" "isogroupOFAVBQ_DN118079_c0_g1"
# [3,] "isogroupOFAVBQ_DN205548_c0_g1" "isogroupOFAVBQ_DN196808_c1_g2"
# [4,] "isogroupOFAVBQ_DN206940_c2_g2" "isogroupOFAVBQ_DN203360_c3_g2"
# [5,] "isogroupOFAVBQ_DN207439_c2_g1" "isogroupOFAVBQ_DN203405_c2_g1"
# [6,] "isogroupOFAVBQ_DN212504_c1_g1" "isogroupOFAVBQ_DN210180_c2_g2"
# [7,] "isogroupOFAVBQ_DN214979_c0_g1" "isogroupOFAVBQ_DN214046_c1_g3"
# [8,] "isogroupOFAVBQ_DN215260_c2_g3" "isogroupOFAVBQ_DN215917_c6_g1"
# [9,] "isogroupOFAVBQ_DN220734_c1_g2" "isogroupOFAVBQ_DN223273_c0_g1"
# [10,] "isogroupOFAVBQ_DN223324_c1_g1" "isogroupOFAVBQ_DN65173_c0_g1" 
# brown                           cyan                           
# [1,] "isogroupOFAVBQ_DN189162_c0_g1" "isogroupOFAVBQ_DN197987_c8_g4"
# [2,] "isogroupOFAVBQ_DN196211_c0_g1" "isogroupOFAVBQ_DN203690_c2_g7"
# [3,] "isogroupOFAVBQ_DN198065_c3_g1" "isogroupOFAVBQ_DN203774_c7_g1"
# [4,] "isogroupOFAVBQ_DN203850_c4_g1" "isogroupOFAVBQ_DN204704_c2_g1"
# [5,] "isogroupOFAVBQ_DN205455_c0_g1" "isogroupOFAVBQ_DN207158_c4_g5"
# [6,] "isogroupOFAVBQ_DN209721_c0_g1" "isogroupOFAVBQ_DN208767_c3_g4"
# [7,] "isogroupOFAVBQ_DN216108_c0_g1" "isogroupOFAVBQ_DN210060_c0_g1"
# [8,] "isogroupOFAVBQ_DN223830_c3_g4" "isogroupOFAVBQ_DN211318_c2_g4"
# [9,] "isogroupOFAVBQ_DN36513_c0_g1"  "isogroupOFAVBQ_DN214858_c9_g2"
# [10,] "isogroupOFAVBQ_DN41324_c0_g2"  "isogroupOFAVBQ_DN216465_c0_g1"
# green                           greenyellow                    
# [1,] "isogroupOFAVBQ_DN180621_c0_g1" "isogroupOFAVBQ_DN198879_c2_g4"
# [2,] "isogroupOFAVBQ_DN194753_c0_g1" "isogroupOFAVBQ_DN200134_c2_g1"
# [3,] "isogroupOFAVBQ_DN197912_c0_g1" "isogroupOFAVBQ_DN204228_c0_g1"
# [4,] "isogroupOFAVBQ_DN198782_c0_g1" "isogroupOFAVBQ_DN204444_c5_g5"
# [5,] "isogroupOFAVBQ_DN200267_c0_g1" "isogroupOFAVBQ_DN205975_c5_g2"
# [6,] "isogroupOFAVBQ_DN201059_c2_g2" "isogroupOFAVBQ_DN208743_c4_g2"
# [7,] "isogroupOFAVBQ_DN208089_c0_g2" "isogroupOFAVBQ_DN213543_c7_g2"
# [8,] "isogroupOFAVBQ_DN208695_c4_g3" "isogroupOFAVBQ_DN218245_c4_g3"
# [9,] "isogroupOFAVBQ_DN211108_c0_g2" "isogroupOFAVBQ_DN222736_c0_g1"
# [10,] "isogroupOFAVBQ_DN219400_c0_g1" "isogroupOFAVBQ_DN328898_c0_g1"
# grey                            magenta                        
# [1,] "isogroupOFAVBQ_DN199552_c1_g2" "isogroupOFAVBQ_DN191751_c1_g1"
# [2,] "isogroupOFAVBQ_DN204434_c4_g2" "isogroupOFAVBQ_DN198838_c4_g1"
# [3,] "isogroupOFAVBQ_DN208787_c0_g1" "isogroupOFAVBQ_DN202544_c6_g3"
# [4,] "isogroupOFAVBQ_DN209674_c1_g1" "isogroupOFAVBQ_DN205781_c1_g4"
# [5,] "isogroupOFAVBQ_DN216078_c0_g2" "isogroupOFAVBQ_DN208579_c2_g3"
# [6,] "isogroupOFAVBQ_DN219778_c0_g1" "isogroupOFAVBQ_DN210734_c2_g1"
# [7,] "isogroupOFAVBQ_DN220244_c3_g2" "isogroupOFAVBQ_DN213809_c5_g1"
# [8,] "isogroupOFAVBQ_DN221390_c0_g1" "isogroupOFAVBQ_DN218885_c2_g7"
# [9,] "isogroupOFAVBQ_DN223590_c0_g1" "isogroupOFAVBQ_DN222846_c0_g2"
# [10,] "isogroupOFAVBQ_DN224683_c1_g1" "isogroupOFAVBQ_DN224689_c3_g1"
# pink                            purple                         
# [1,] "isogroupOFAVBQ_DN196696_c0_g2" "isogroupOFAVBQ_DN17990_c0_g1" 
# [2,] "isogroupOFAVBQ_DN202581_c8_g1" "isogroupOFAVBQ_DN198862_c2_g1"
# [3,] "isogroupOFAVBQ_DN204219_c1_g1" "isogroupOFAVBQ_DN206957_c0_g1"
# [4,] "isogroupOFAVBQ_DN206523_c2_g1" "isogroupOFAVBQ_DN210697_c0_g1"
# [5,] "isogroupOFAVBQ_DN206703_c0_g4" "isogroupOFAVBQ_DN211158_c3_g2"
# [6,] "isogroupOFAVBQ_DN209674_c1_g1" "isogroupOFAVBQ_DN217447_c0_g1"
# [7,] "isogroupOFAVBQ_DN212072_c0_g3" "isogroupOFAVBQ_DN219516_c1_g2"
# [8,] "isogroupOFAVBQ_DN215819_c0_g2" "isogroupOFAVBQ_DN220962_c0_g1"
# [9,] "isogroupOFAVBQ_DN220094_c0_g3" "isogroupOFAVBQ_DN221701_c1_g1"
# [10,] "isogroupOFAVBQ_DN224619_c1_g1" "isogroupOFAVBQ_DN222609_c3_g2"
# red                             salmon                         
# [1,] "isogroupOFAVBQ_DN197511_c2_g3" "isogroupOFAVBQ_DN137595_c0_g1"
# [2,] "isogroupOFAVBQ_DN199946_c5_g1" "isogroupOFAVBQ_DN205820_c4_g1"
# [3,] "isogroupOFAVBQ_DN200871_c0_g1" "isogroupOFAVBQ_DN215351_c1_g1"
# [4,] "isogroupOFAVBQ_DN204496_c3_g1" "isogroupOFAVBQ_DN216608_c4_g1"
# [5,] "isogroupOFAVBQ_DN204827_c6_g4" "isogroupOFAVBQ_DN217266_c0_g1"
# [6,] "isogroupOFAVBQ_DN217563_c1_g1" "isogroupOFAVBQ_DN220450_c0_g1"
# [7,] "isogroupOFAVBQ_DN218804_c4_g8" "isogroupOFAVBQ_DN220982_c3_g2"
# [8,] "isogroupOFAVBQ_DN221071_c0_g1" "isogroupOFAVBQ_DN221306_c7_g5"
# [9,] "isogroupOFAVBQ_DN221385_c0_g2" "isogroupOFAVBQ_DN224860_c3_g1"
# [10,] "isogroupOFAVBQ_DN223832_c4_g1" "isogroupOFAVBQ_DN245573_c0_g1"
# tan                             turquoise                      
# [1,] "isogroupOFAVBQ_DN169978_c0_g1" "isogroupOFAVBQ_DN206162_c0_g1"
# [2,] "isogroupOFAVBQ_DN195745_c0_g1" "isogroupOFAVBQ_DN213863_c1_g1"
# [3,] "isogroupOFAVBQ_DN197166_c0_g1" "isogroupOFAVBQ_DN218184_c1_g1"
# [4,] "isogroupOFAVBQ_DN19927_c0_g1"  "isogroupOFAVBQ_DN218223_c9_g1"
# [5,] "isogroupOFAVBQ_DN201217_c2_g1" "isogroupOFAVBQ_DN218344_c1_g2"
# [6,] "isogroupOFAVBQ_DN204912_c5_g1" "isogroupOFAVBQ_DN224718_c1_g1"
# [7,] "isogroupOFAVBQ_DN241889_c0_g1" "isogroupOFAVBQ_DN225147_c1_g1"
# [8,] "isogroupOFAVBQ_DN256755_c0_g1" "isogroupOFAVBQ_DN225204_c1_g2"
# [9,] "isogroupOFAVBQ_DN28700_c0_g1"  "isogroupOFAVBQ_DN225272_c1_g1"
# [10,] "isogroupOFAVBQ_DN323680_c0_g1" "isogroupOFAVBQ_DN225399_c5_g1"
# yellow                         
# [1,] "isogroupOFAVBQ_DN199373_c3_g1"
# [2,] "isogroupOFAVBQ_DN200852_c4_g2"
# [3,] "isogroupOFAVBQ_DN201751_c1_g1"
# [4,] "isogroupOFAVBQ_DN206790_c4_g2"
# [5,] "isogroupOFAVBQ_DN211653_c2_g1"
# [6,] "isogroupOFAVBQ_DN215093_c1_g1"
# [7,] "isogroupOFAVBQ_DN221306_c7_g5"
# [8,] "isogroupOFAVBQ_DN224809_c0_g1"
# [9,] "isogroupOFAVBQ_DN245573_c0_g1"
# [10,] "isogroupOFAVBQ_DN251696_c0_g1"

save(modules2019, colors2019, ME1_2019, ME2_2021, file = "MEoutput.RData")

#### So far, nothing has been super informative... Gonna try standard techniques used for WGCNA below
# Correlate traits

load("MEoutput.RData")

nSamples2019 = nrow(datExprOut2019)
nSamples2021 = nrow(datExprOut2021)

#Recalculate MEs with color labels
MEs2019 = orderMEs(ME1_2019)
MEs2021 = orderMEs(ME2_2021)
MEs2021 = MEs2021[,colnames(MEs2019)]

names(datTraitsOut2019) = c("CR x CR", "HR x HR", "Hybrid", "Ctrl", "Heat")
moduleTraitCor2019 = cor(MEs2019, datTraitsOut2019, use = "p")

datTraitsOut2021_2C = read.csv("datTraitsOut2021_2Crosses.csv",row.names=1,check.names=FALSE)
names(datTraitsOut2021_2C) = c("HR x HR", "CR x HR", "HR x CR", "Ctrl", "Heat")
moduleTraitCor2021_2C = cor(MEs2021, datTraitsOut2021_2C, use = "p")

moduleTraitPvalue2019 = corPvalueStudent(moduleTraitCor2019, nSamples2019)
moduleTraitPvalue2021_2C = corPvalueStudent(moduleTraitCor2021_2C, nSamples2021)

#Print correlation heatmap between modules and traits
textMatrix2019= paste(signif(moduleTraitCor2019, 2), "\n(", 
                  signif(moduleTraitPvalue2019, 1), ")", sep= "")
dim(textMatrix2019)= dim(moduleTraitCor2019)

textMatrix2021_2C= paste(signif(moduleTraitCor2021_2C, 2), "\n(",
                      signif(moduleTraitPvalue2021_2C, 1), ")", sep= "")
dim(textMatrix2021_2C)= dim(moduleTraitCor2021_2C)

quartz()
par(mfrow=c(1,2))
par(mar = c(6, 8.5, 1.5, 1.5))

labeledHeatmap(Matrix = moduleTraitCor2019,
               xLabels = names(datTraitsOut2019),
               yLabels = names(MEs2019),
               ySymbols = names(MEs2019),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix2019,
               setStdMargins = FALSE,
               # plotLegend = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("2019"))

labeledHeatmap(Matrix = moduleTraitCor2021_2C,
               xLabels = names(datTraitsOut2021_2C),
               yLabels = names(MEs2021),
               ySymbols = names(MEs2021),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix2021_2C,
               setStdMargins = FALSE,
               #plotLegend = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("2021"))

###########plotting massive table of all information - module membership, genes, gene names, etc.
annot=read.table("ofav_iso2gene_uni.tab",sep="\t",quote="")
iso2go=read.table("Updated_iso2go_combined.tab",sep="\t",quote="")

probes=colnames(datExprOut)

probes2go=match(probes,iso2go$V1)
summary(probes2go) # 14731 genes annotated
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#       3    9372   18638   18503   27678   36143   13944 

probes2annot = match(probes,annot$V1)
summary(probes2annot) # 14058 genes annotated
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#       1    8696   17618   17449   26028   34687   14617 

datGS.Traits2019=data.frame(cor(datExprOut2019,datTraitsOut2019,use="p"))
names(datGS.Traits2019)=paste("cor",names(datGS.Traits2019),sep=".")
datME2019=moduleEigengenes(datExprOut2019,modules2019)$eigengenes
datKME2019=signedKME(datExprOut2019, datME2019, outputColumnName="MM.")
datOutput2019=data.frame(ProbeID=names(datExprOut2019),annot[probes2annot,],modules2019,datKME2019,datGS.Traits2019)

datGS.Traits2021=data.frame(cor(datExprOut2021,datTraitsOut2021,use="p"))
names(datGS.Traits2021)=paste("cor",names(datGS.Traits2021),sep=".")
datME2021=moduleEigengenes(datExprOut2021,modules2019)$eigengenes
datKME2021=signedKME(datExprOut2021, datME2021, outputColumnName="MM.")
datOutput2021=data.frame(ProbeID=names(datExprOut2021),annot[probes2annot,],modules2019,datKME2021,datGS.Traits2021)

Modcolors=c("pink", "purple", "magenta", "red", "yellow", "cyan", "black", "salmon", "greenyellow", "turquoise") #change to vector of desired module color names

# Generate categorical data for GO analysis 
# No need to do separate ones for different years bc they have the same module assignment lol
for (col in Modcolors) {
  tab=datOutput2021[,c(1,4)]
  
  tab$modules2019=as.character(tab$modules2019) # Don't change the year as the col name is the same
  #Categorical GeneOntology by Module
  tab$modules2019[tab$modules2019!=col]<-0
  tab$modules2019[tab$modules2019==col]<-1 
  tab$modules2019=as.factor(tab$modules2019) 
  print(col)
  print(summary(tab)) #do counts match table of module colors?
  print(head(tab))
  
  write.csv(tab,file=paste("GO_MM",col,"_categorical.csv", sep=""),quote=F,row.names=F)
}