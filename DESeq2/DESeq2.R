library(DESeq2)
library(pheatmap)
library(vsn)
library(RColorBrewer)
library(gplots) # for venn diagram
library(ggplot2)
library(VennDiagram)

cts=read.csv("AllCountsNH2020.csv",header=TRUE,row.names=1) #Reading in the table of counts per isogroup by sample
head(cts) 
length(cts[,1])  #67239 isogroups

cts=cts[,-32] # Remove the one heat sample from AS

names(cts)<-gsub("X","",names(cts))
names(cts)

######## read in trait data
colData=read.csv('colData.csv',row.names = 1)
colData=colData[-32,] # Remove the one heat sample from AS

head(colData)
colData$DevDay=as.factor(colData$DevDay)
colData$Type=as.factor(colData$Type)
colData$Trmt=as.factor(colData$Trmt)

summary(colData)

#Subset data to only include MS samples
cts_ms=cts[,c(32:77)]
colData_ms=colData[c(32:77),]

cts_ms$low = apply(cts_ms[,1:46],1,function(x){sum(x<=2)})  #making new column counting number of samples with counts <=2 within each isogroup (host) 
cts_ms<-cts_ms[-which(cts_ms$low>41.4),] #41.4 is ~90% of 46 samples - get rid of count less than 2 in more than 90% of samples
nrow(cts_ms) #27894 isogroups

#Check for outliers with too few reads
readsleft=c()
for (column in names(cts)) {
  val=sum(cts[,column])
  readsleft=append(readsleft,val)}

RLtable=data.frame(cbind(names(cts),readsleft)) #Looked fine

save(cts, colData, cts_ms, colData_ms, cts_devday, colData_devday, file="Counts&ColData_LowCountsRemoved.RData")

####check for treatment effect
dds3 <- DESeqDataSetFromMatrix(countData = cts_ms[,1:46],colData = colData_ms,design = ~ Trmt)
dds3 <- DESeq(dds3)

res3 <- results(dds3) #Heat vs. Ctrl
summary(res3) #561 Up, 436 Down

p.adj.cutoff=0.1
# add T/F column if gene is sig
res3$threshold <- as.logical(res3$padj < p.adj.cutoff)

# Sort the results tables
res3_sorted <- res3[order(res3$padj), ]

# Get significant genes
sigTrmt <- row.names(res3_sorted)[which(res3_sorted$threshold)]

##############
# combine the factors of interest into a single factor with all combinations of the original factors
dds4 <- DESeqDataSetFromMatrix(countData = cts_ms[,1:46],colData = colData_ms,design = ~ 1)
dds4$group <- factor(paste0(dds4$Trmt, dds4$Type))

# change the design to include just this factor, e.g. ~ group
design(dds4) <- ~ group
dds4 <- DESeq(dds4)
res4=results(dds4)
summary(res4) # 405 up, 301 down, 7571 low count

res_inshoretrmt <- results(dds4, contrast=c("group","HeatInshore","CtrlInshore"))
summary(res_inshoretrmt) # 133 up, 144 down

res_offshoretrmt <- results(dds4, contrast=c("group","HeatOffshore","CtrlOffshore"))
summary(res_offshoretrmt) # 376 up, 158 down

res_crosstrmt <- results(dds4, contrast=c("group", "HeatCross", "CtrlCross"))
summary(res_crosstrmt) # 12 up, 9 down

res_origin_control <- results(dds4, contrast=c('group', 'CtrlInshore', 'CtrlOffshore'))
summary(res_origin_control) # 991 up, 1023 down

p.adj.cutoff=0.1
# add T/F column if gene is sig
res_inshoretrmt$threshold <- as.logical(res_inshoretrmt$padj < p.adj.cutoff)
res_offshoretrmt$threshold <- as.logical(res_offshoretrmt$padj < p.adj.cutoff)
res_crosstrmt$threshold <- as.logical(res_crosstrmt$padj < p.adj.cutoff)

# Sort the results tables
res_inshoretrmt_sorted <- res_inshoretrmt[order(res_inshoretrmt$padj), ]
res_offshoretrmt_sorted <- res_offshoretrmt[order(res_offshoretrmt$padj), ]
res_crosstrmt_sorted <- res_crosstrmt[order(res_crosstrmt$padj), ]

# Get significant genes
siginshoretrmt <- row.names(res_inshoretrmt_sorted)[which(res_inshoretrmt_sorted$threshold)]
sigoffshoretrmt <- row.names(res_offshoretrmt_sorted)[which(res_offshoretrmt_sorted$threshold)]
sigcrosstrmt <- row.names(res_crosstrmt_sorted)[which(res_crosstrmt_sorted$threshold)]
siginshoretrmtOnly <- siginshoretrmt[!siginshoretrmt %in% c(sigoffshoretrmt,sigcrosstrmt)]
sigoffshoretrmtOnly <- sigoffshoretrmt[!sigoffshoretrmt %in% c(siginshoretrmt,sigcrosstrmt)]
sigcrosstrmtOnly <- sigcrosstrmt[!sigcrosstrmt %in% c(siginshoretrmt,sigoffshoretrmt)]
#sigTrmtStage <- sigTrmt[sigTrmt %in% sigStage]
#sigAll <- sigTrmtStage[sigTrmtStage %in% sigOrigin]

candidates=list('InTrmt'=siginshoretrmt,'OffTrmt'=sigoffshoretrmt,'HybTrmt'=sigcrosstrmt)
quartz()
venn(candidates)
# make pretty venn diagram
venn.diagram(
  x=candidates,
  category.names = c('InTrmt','OffTrmt'),
  filename='~/My Drive/NOAA_USC_Ofav/GeneExp/DESeq2/TypeTrmt_Venn.png',
  col=c("#440154ff", '#21908dff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  cex = 2,
  #cat.col = c("#440154ff", '#21908dff', 'black'),
  cat.cex = 2,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27),
  cat.dist = c(0.08, 0.08)
)

#Try a different package, didn't seem to work
library(venneuler)
v <- venneuler(c(InTrmt=277, OffTrmt=534, HybTrmt=21,"InTrmt&OffTrmt"=22,"InTrmt&HybTrmt"=7,"OffTrmt&HybTrmt"=4,"InTrmt&OffTrmt&HybTrmt"=4))
plot(v)

# transform counts for plotting
vsd=vst(dds4,blind=FALSE)
norm_df=assay(vsd)
means=apply(norm_df,1,mean) # means of rows
norm_df_scaled=norm_df-means #rescale expression data so it's up and down
df=as.data.frame(colData(dds4)[,c("Type","Trmt")])
display.brewer.pal(n = 7, name = 'Set1')
brewer.pal(n = 7, name = "Set1")
ann_colors=list(
  Trmt=c(Ctrl='#377EB8', Heat='#E41A1C'),
  Type=c(Inshore='#FF7F00',Offshore='#4DAF4A',Cross='#984EA3')
)

#Generate GO inputs
res_inshoretrmt$direction=ifelse(res_inshoretrmt$log2FoldChange>0,1,-1) #red=higher expression in heat
res_inshoretrmt$logP<-(-log((res_inshoretrmt$padj+0.0000000001),10))*(res_inshoretrmt$direction)

res_offshoretrmt$direction=ifelse(res_offshoretrmt$log2FoldChange>0,1,-1) #red=higher expression in heat
res_offshoretrmt$logP<-(-log((res_offshoretrmt$padj+0.0000000001),10))*(res_offshoretrmt$direction)

res_crosstrmt$direction=ifelse(res_crosstrmt$log2FoldChange>0,1,-1) #red=higher expression in heat
res_crosstrmt$logP<-(-log((res_crosstrmt$padj+0.0000000001),10))*(res_crosstrmt$direction)

inshoretrmt_out<-as.data.frame(cbind("gene"=row.names(res_inshoretrmt),"logP"=res_inshoretrmt$logP))
offshoretrmt_out<-as.data.frame(cbind("gene"=row.names(res_offshoretrmt),"logP"=res_offshoretrmt$logP))
crosstrmt_out<-as.data.frame(cbind("gene"=row.names(res_crosstrmt),"logP"=res_crosstrmt$logP))

write.csv(inshoretrmt_out,file="GOrank_inshoretrmt.csv",quote=F,row.names=F)
write.csv(offshoretrmt_out,file="GOrank_offshoretrmt.csv",quote=F,row.names=F)
write.csv(crosstrmt_out,file="GOrank_crosstrmt.csv",quote=F,row.names=F)

###############################################
####Front loading/back loading
res_inshoretrmt = as.data.frame(res_inshoretrmt)
res_offshoretrmt = as.data.frame(res_offshoretrmt)
res_origin_control = as.data.frame(res_origin_control)

# save(res_inshoretrmt, res_offshoretrmt, res_origin_control, file="Front_Back_Loading_2019_datasets.RData")

load("Front_Back_Loading_2019_datasets.RData")

iso2gene=read.table("ofav_iso2gene_uni.tab",sep="\t",quote="")

In_trmtGenes=rownames(subset(res_inshoretrmt,padj < 0.1))
Off_trmtGenes=rownames(subset(res_offshoretrmt,padj < 0.1))

# let's check if there are any offshore DEGs being frontloaded / backloaded
# aka look at if any of the treatment responsive genes not responding in offshore larvae are upregulated relative to adult samples in control conditions
IngenesOnlyUP=rownames(subset(res_inshoretrmt,padj < 0.1 & log2FoldChange > 0  & !rownames(res_inshoretrmt) %in% Off_trmtGenes))
OffgenesOnlyUP=rownames(subset(res_offshoretrmt,padj < 0.1 & log2FoldChange > 0  & !rownames(res_offshoretrmt) %in% In_trmtGenes))

IngenesOnlyDown=rownames(subset(res_inshoretrmt,padj < 0.1 & log2FoldChange < 0  & !rownames(res_inshoretrmt) %in% Off_trmtGenes))
OffgenesOnlyDown=rownames(subset(res_offshoretrmt,padj < 0.1 & log2FoldChange < 0  & !rownames(res_offshoretrmt) %in% In_trmtGenes))

# get genes upreg or downreg in offshore relative to inshore
up_control=subset(res_origin_control,padj < 0.1 & log2FoldChange < 0) #upreg in offshore rel to inshore
down_control=subset(res_origin_control,padj < 0.1 & log2FoldChange > 0) #downreg in offshore rel to inshore

origin_control_up=rownames(up_control)
origin_control_down=rownames(down_control)

# frontloading in offshore
candidates=list('inshore heat response up' = IngenesOnlyUP,'offshore vs inshore control up' = origin_control_up)
venn(candidates) # 28 genes frontloaded by offshore, but responding in inshore
off_front=iso2gene[iso2gene$V1 %in% IngenesOnlyUP[IngenesOnlyUP %in% origin_control_up],] # 11 annotated

# backloading in offshore
candidates=list('inshore heat response down' = IngenesOnlyDown,'offshore vs inshore control down' = origin_control_down)
venn(candidates) # 47 genes backloaded by offshore, but responding in inshore
off_back=iso2gene[iso2gene$V1 %in% IngenesOnlyDown[IngenesOnlyDown %in% origin_control_down],] # 31 annotated

####################################
# Now for inshore
up_control=subset(res_origin_control,padj < 0.1 & log2FoldChange > 0)
down_control=subset(res_origin_control,padj < 0.1 & log2FoldChange < 0)

origin_control_up=rownames(up_control)
origin_control_down=rownames(down_control)

# frontloading in inshore
candidates=list('offshore heat response up' = OffgenesOnlyUP,'inshore vs offshore control up' = origin_control_up)
venn(candidates) # 82 genes frontloaded by inshore, but responding in offshore
in_front=iso2gene[iso2gene$V1 %in% OffgenesOnlyUP[OffgenesOnlyUP %in% origin_control_up],] # 42 annotated

# backloading in inshore
candidates=list('offshore heat response down' = OffgenesOnlyDown,'inshore vs offshore control down' = origin_control_down)
venn(candidates) # 23 genes frontloaded by inshore, but responding in offshore
in_back=iso2gene[iso2gene$V1 %in% OffgenesOnlyDown[OffgenesOnlyDown %in% origin_control_down],] # 14 annotated




###############################################
#### Analyze 2021 dataset- only offshore larvae
cts2021=read.csv('AllCountsNH2022.csv',row.names=1,check.names=FALSE)
names(cts2021)<-gsub("-","_",names(cts2021))

colData=read.csv('colData_Combined.csv',row.names = 1)
traits2021 = colData[79:108,c(1,3)]
traits2021$Type=as.factor(traits2021$Type)
traits2021$Trmt=as.factor(traits2021$Trmt)
summary(traits2021)
rownames(traits2021)<-gsub("-","_",rownames(traits2021))
cts2021 = cts2021[,rownames(traits2021)]

table(colnames(cts2021)==rownames(traits2021)) #TRUE 30

#Remove isogroups with low counts from dataset - count less than 2 in more than 90% of samples
cts2021$low = apply(cts2021[,1:30],1,function(x){sum(x<=2)})  #making new column counting number of samples with counts <=2 within each isogroup (host) 
cts2021<-cts2021[-which(cts2021$low>27),] #27 is ~90% of 30 samples - get rid of count less than 2 in more than 90% of samples
nrow(cts2021) #28960 isogroups

#Check for outliers with too few reads
readsleft=c()
for (column in names(cts2021)) {
  val=sum(cts2021[,column])
  readsleft=append(readsleft,val)}

RLtable=data.frame(cbind(names(cts2021),readsleft)) # YZ_6 has half a million (lowest count)

save(cts2021, traits2021, file="Larvae2021_Counts&ColData_LowCountsRemoved.RData")

#### DESeq analysis for 2021 larvae
# combine the factors of interest into a single factor with all combinations of the original factors
dds <- DESeqDataSetFromMatrix(countData = cts2021[,1:30],colData = traits2021,design = ~ 1)
dds$group <- factor(paste0(dds$Trmt, dds$Type))

# change the design to include just this factor, e.g. ~ group
design(dds) <- ~ group
dds <- DESeq(dds)

res=results(dds)
res
summary(res) # 1374 up, 1171 down, 3369 low count

res_offshoretrmt <- results(dds, contrast=c("group","HeatOffshore","CtrlOffshore"))
summary(res_offshoretrmt) # 706 up, 952 down

res_cross1trmt <- results(dds, contrast=c("group","HeatCross1","CtrlCross1"))
summary(res_cross1trmt) # 228 up, 54 down

res_cross2trmt <- results(dds, contrast=c("group","HeatCross2","CtrlCross2"))
summary(res_cross2trmt) # 124 up, 101 down

p.adj.cutoff=0.1

# add T/F column if gene is sig
res_offshoretrmt$threshold <- as.logical(res_offshoretrmt$padj < p.adj.cutoff)
res_cross1trmt$threshold <- as.logical(res_cross1trmt$padj < p.adj.cutoff)
res_cross2trmt$threshold <- as.logical(res_cross2trmt$padj < p.adj.cutoff)

# Sort the results tables
res_offshoretrmt_sorted <- res_offshoretrmt[order(res_offshoretrmt$padj), ]
res_cross1trmt_sorted <- res_cross1trmt[order(res_cross1trmt$padj), ]
res_cross2trmt_sorted <- res_cross2trmt[order(res_cross2trmt$padj), ]

# Get significant genes
sigoffshoretrmt <- row.names(res_offshoretrmt_sorted)[which(res_offshoretrmt_sorted$threshold)]
sigcross1trmt <- row.names(res_cross1trmt_sorted)[which(res_cross1trmt_sorted$threshold)]
sigcross2trmt <- row.names(res_cross2trmt_sorted)[which(res_cross2trmt_sorted$threshold)]

#Generate GO inputs
res_offshoretrmt$direction=ifelse(res_offshoretrmt$log2FoldChange>0,1,-1) #red=higher expression in heat
res_offshoretrmt$logP<-(-log((res_offshoretrmt$padj+0.0000000001),10))*(res_offshoretrmt$direction)

res_cross1trmt$direction=ifelse(res_cross1trmt$log2FoldChange>0,1,-1) #red=higher expression in heat
res_cross1trmt$logP<-(-log((res_cross1trmt$padj+0.0000000001),10))*(res_cross1trmt$direction)

res_cross2trmt$direction=ifelse(res_cross2trmt$log2FoldChange>0,1,-1) #red=higher expression in heat
res_cross2trmt$logP<-(-log((res_cross2trmt$padj+0.0000000001),10))*(res_cross2trmt$direction)

offshoretrmt_out<-as.data.frame(cbind("gene"=row.names(res_offshoretrmt),"logP"=res_offshoretrmt$logP))
cross1trmt_out<-as.data.frame(cbind("gene"=row.names(res_cross1trmt),"logP"=res_cross1trmt$logP))
cross2trmt_out<-as.data.frame(cbind("gene"=row.names(res_cross2trmt),"logP"=res_cross2trmt$logP))

write.csv(offshoretrmt_out,file="GOrank_offshoretrmt_2021.csv",quote=F,row.names=F)
write.csv(cross1trmt_out,file="GOrank_cross1trmt_2021.csv",quote=F,row.names=F)
write.csv(cross2trmt_out,file="GOrank_cross2trmt_2021.csv",quote=F,row.names=F)


###########################
################### DAPC
library('adegenet')
library('dplyr')
library(wesanderson)
library('stringr')
library(ggridges)
library(cowplot)

# For 2019

load("Counts&ColData_LowCountsRemoved.RData")
cts2019=cts_ms
traits2019=colData_ms[, c(1,3)]

dds0 <- DESeqDataSetFromMatrix(countData = cts2019[,1:46],colData = traits2019,design = ~ 1)
design(dds0) <- formula(~ Trmt + Type)

dds0 <- DESeq(dds0)
res0=results(dds0)
res0
summary(res0) # 222 up, 336 down, 3369 low count

vsd0=vst(dds0)
vsd_df0=assay(vsd0)

dim(vsd_df0) #27894    46

traits2019$groupall=as.factor(paste0(traits2019$Type, traits2019$Trmt))
traits2019$groupall=sub('CrossCtrl','Hybrid Ctrl',traits2019$groupall)
traits2019$groupall=sub('CrossHeat','Hybrid Heat',traits2019$groupall)
traits2019$groupall=sub('InshoreCtrl','CR x CR Ctrl',traits2019$groupall)
traits2019$groupall=sub('InshoreHeat','CR x CR Heat',traits2019$groupall)
traits2019$groupall=sub('OffshoreCtrl','HR x HR Ctrl',traits2019$groupall)
traits2019$groupall=sub('OffshoreHeat','HR x HR Heat',traits2019$groupall)

dapc0 <- dapc(t(vsd_df0), traits2019$groupall) # 20 PCs, 4 discriminant functions
temp <- optim.a.score(dapc0, n.sim = 5)
#for the vsd_df, they suggest retaining 8 PCs

dapc <- dapc(t(vsd_df0), traits2019$groupall, n.da=2, n.pca=8)
varexpl <- round((dapc$eig/sum(dapc$eig))[1:2] * 100, 1) # 77.3 14.9

dapc1 <- tibble(sample = rownames(dapc$ind.coord),
                grp = dapc$grp,
                LD1 = dapc$ind.coord[,1],
                LD2 = dapc$ind.coord[,2])
dapc2 <- dapc1 %>%
  group_by(grp) %>%
  summarize(c1 = mean(LD1),
            c2 = mean(LD2)) %>%
  full_join(dapc1)

traits2019$sample = rownames(traits2019)
dapc2.joined = left_join(dapc2, traits2019)
dapc2.joined$Type = factor(dapc2.joined$Type, levels = c("Inshore", "Offshore", "Cross"))
dapc2.joined$Type=sub('Inshore','CR x CR',dapc2.joined$Type)
dapc2.joined$Type=sub('Offshore','HR x HR',dapc2.joined$Type)
dapc2.joined$Type=sub('Cross','Hybrid',dapc2.joined$Type)
dapc2.joined_subset = dapc2.joined[1:34,] # only plotted inshore and offshore
dapc2.joined_subset$Type = factor(dapc2.joined_subset, levels = c("CR x CR", "HR x HR"))

g1 = ggplot(dapc2.joined, aes(x=LD1, y=LD2,colour=Type,shape=Trmt)) +
  geom_point(aes(colour=Type,shape=Trmt),size=3)+
  scale_color_manual(values=c("#FC4E07","#00AFBB","#E7B800"),
                     labels=c("CR x CR", "HR x HR", "Hybrid"))+
  scale_shape_manual(values = c(1,16))+
  labs(x = paste0("LD1 [", varexpl[1],"%]"), y = paste0("LD2 [", varexpl[2],"%]"))+
  #theme_bw()+
  guides(shape = guide_legend(order = 2),col = guide_legend(order = 1))+
  ggtitle("2019")

quartz()

ggplot(dapc2.joined_subset,aes(x=LD2, y=Type, fill=Trmt, height = ..density..))+
  geom_density_ridges(scale = 1, stat = "density") +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_fill_manual(name="Treatment", values =c(alpha("darkgrey",0.2),alpha("darkgrey", 0.9))) +
  theme_ridges()

# Now 2021 larvae

load("Larvae2021_Counts&ColData_LowCountsRemoved.RData")

dds <- DESeqDataSetFromMatrix(countData = cts2021[,1:30],colData = traits2021,design = ~ 1)
design(dds) <- formula(~ Trmt + Type)

dds <- DESeq(dds)
res=results(dds)
res
summary(res) # 1761 up, 1758 down, 3369 low count

vsd=vst(dds)
vsd_df=assay(vsd)

dim(vsd_df) #28960    30

traits2021$groupall=as.factor(paste0(traits2021$Type, traits2021$Trmt))
traits2021$groupall=sub('Cross1Ctrl','CR x HR Ctrl',traits2021$groupall)
traits2021$groupall=sub('Cross1Heat','CR x HR Heat',traits2021$groupall)
traits2021$groupall=sub('Cross2Ctrl','HR x CR Ctrl',traits2021$groupall)
traits2021$groupall=sub('Cross2Heat','HR x CR Heat',traits2021$groupall)
traits2021$groupall=sub('OffshoreCtrl','HR x HR Ctrl',traits2021$groupall)
traits2021$groupall=sub('OffshoreHeat','HR x HR Heat',traits2021$groupall)


dapc0 <- dapc(t(vsd_df), traits2021$groupall)
temp <- optim.a.score(dapc0, n.sim = 5)
#for the vsd_df, they suggest retaining 18 PCs

dapc <- dapc(t(vsd_df), traits2021$groupall, n.da=2, n.pca=4)
varexpl <- round((dapc$eig/sum(dapc$eig))[1:2] * 100, 1) # 92.2 6.8

dapc1 <- tibble(sample = rownames(dapc$ind.coord),
                grp = dapc$grp,
                LD1 = dapc$ind.coord[,1],
                LD2 = dapc$ind.coord[,2])
dapc2 <- dapc1 %>%
  group_by(grp) %>%
  summarize(c1 = mean(LD1),
            c2 = mean(LD2)) %>%
  full_join(dapc1)

traits2021$sample = rownames(traits2021)
dapc2.joined = left_join(dapc2, traits2021)
dapc2.joined$Type=sub('Offshore','HR x HR',dapc2.joined$Type)
dapc2.joined$Type=sub('Cross1','CR x HR',dapc2.joined$Type)
dapc2.joined$Type=sub('Cross2','HR x CR',dapc2.joined$Type)
dapc2.joined$grp = factor(dapc2.joined$grp, levels = c("HR x HR Ctrl", "HR x HR Heat", "CR x HR Ctrl", "CR x HR Heat","HR x CR Ctrl", "HR x CR Heat"))
dapc2.joined$LD2 = dapc2.joined$LD2 * -1
dapc2.joined <- dapc2.joined[order(dapc2.joined$Trmt,dapc2.joined$Type),]

g2 = ggplot(dapc2.joined, aes(x=LD1, y=LD2,colour=Type,shape=Trmt)) +
  geom_point(aes(colour=Type,shape=Trmt),size=3)+
  scale_color_manual(values=c("#00AFBB","#FCD2B2","#FF8D33"),
                     labels=c("HR x HR", "CR x HR", "HR x CR"))+
  scale_shape_manual(values = c(1,16))+
  labs(x = paste0("LD1 [", varexpl[1],"%]"), y = paste0("LD2 [", varexpl[2],"%]"))+
  #theme_bw()+
  guides(shape = guide_legend(order = 2),col = guide_legend(order = 1))+
  ggtitle("2021")

ggplot(dapc2.joined,aes(x=LD2, y=Type, fill=Trmt, height = ..density..))+
  geom_density_ridges(scale = 1, stat = "density") +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_fill_manual(name="Trmt", values =c(alpha("darkgrey",0.2),alpha("darkgrey", 0.9))) +
  theme_ridges()

quartz()
plot_grid(g1, g2, labels = c("a", "b"))
