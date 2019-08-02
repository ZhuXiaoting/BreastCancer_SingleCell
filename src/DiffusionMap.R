library(RColorBrewer)
library(scran)
library(destiny)
library(cowplot)
library(ggplot2)
library(viridis)
library(scatterplot3d)
library(dplyr)

source("../src/functions.R")

#--Load Data-------------------------------------
dataList <- readRDS("../data/Robjects/ExpressionList_QC_norm.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]
# 
# # ---- Prepocessing ----
# 
# Cells and genes 
fD$keep <- rowMeans(m) > 0.01
m <- m[fD$keep,pD$PassAll]
pD <- pD[pD$PassAll,]
m <- m[,pD$SampleID !="4T1"]

pD <- pD[pD$SampleID !="4T1",]
fD <- fD[fD$keep,]

#--Normalize counts by size factor---------------
m.norm <- t(t(m))/pD$sf
rownames(m.norm) <- fD$symbol
fD$gene_short_name <- fD$symbol

m.norm <- m.norm[fD$highVar,]
m.norm <- t(log(m.norm+1))

#--Compute diffusion map-------------------------
set.seed(1000)
dm <- DiffusionMap(m.norm, n_eigs=20, rotate=TRUE)
dms <- eigenvectors(dm)[,1:4]
dms <- data.frame(dms,
                  barcode=pD$barcode)

pD.vp <- left_join(pD,dms,by="barcode")

saveRDS(list(genes=colnames(m.vp),dm=dm),file="../data/Robjects/DiffusionMap.rds")
#-----------------------

t1 <- which.max(dms[,2])
t2 <- which.max(dms[,1])
t3 <- which.min(dms[,1])

#--Compute Pseudotime and branching-------------
set.seed(1000)
dpt <-DPT(dm, tips = c(t1,t2,t3))
root <- which(dpt@tips[,1])[2]
rootdpt <- paste0("DPT",root)

#--Rename branches------------------------------

branch <- dpt@branch[,1]
branch[is.na(branch)]  <- "Intermediate"
branch[branch==1] <- "Root"  #4T1
branch[branch==2] <- "A"  #FF99WT
branch[branch==3] <- "B"  #NEU

#--add branches and pseudotime to pD-----------
dms$branch <- factor(branch,levels=c("Root","Intermediate", "A","B"))
dms$dpt<- dpt[["dpt"]]

write.csv(dms, file="DMS.csv", row.names=FALSE)
saveRDS(list(genes=colnames(m.vp),dm=dm),file="newDiffusionMap.rds")

#------------------
pD <- right_join(pD,dms,by="barcode")
cols <- as.character(pD$Colors)

ggplot(pD.vp, aes(DC1,DC2,color=branch)) +
  geom_point() +
  theme_bw()

scatterplot3d(x=pD[,"DC1"],
              y=-pD[,"DC2"],
              z=pD[,"DC3"],
              color=cols,
              pch=20,
              angle=60,
              cex.symbols=1.0,
              scale.y=0.5,
              mar=c(5,3,-0.1,3)+0.1,
              xlab="Component 1",
              ylab="-Component 2",
              zlab="Component 3",
              box=FALSE,axis=TRUE, tick.marks=TRUE, label.tick.marks=TRUE)

cols <- unique(as.character(pD$Colors))
names(cols) <- unique(as.character(pD$TumorType))

#--Add legend for 3D plot-----------------
clustLeg <- ggplot(pD, aes(x=tSNE1,y=tSNE2,color=TumorType)) +
  geom_point() +
  scale_color_manual(values=cols) +
  theme(legend.direction="horizontal",
        legend.title=element_blank(),
        legend.position="bottom")+
  guides(colour=guide_legend(nrow=2,
                             override.aes=list(size=3)))
clustLeg <- get_legend(clustLeg)

#--Create alternative view inlet---------
p.clust <- ggplot(pD, aes(x=DC1,y=DC2, color=TumorType)) +
  geom_point(size=2, pch=20) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  scale_color_manual(values=cols) +
  guides(colour=FALSE) +
  xlab("Component 1") +
  ylab("Component 2") 

plot_grid(p.clust, clustLeg,nrow=2,rel_heights=c(1,0.1))

#plot.DPT(dpt, dcs= c(2,1),root=1, paths_to=c(2,3), col_by = "branch",col_path = c("red","blue"))
plot(dpt, pal = viridis::magma)

# ---- PlotWithBranchDPT -----#swapped DC1 and DC2 for better looking 
b0 <- pD[pD$branch =="Intermediate",]
pb0 <- ggplot(pD,aes(x=DC1,y=DC2)) +
  geom_point(size=1,color="grey80") +
  geom_point(data=b0,aes(x=DC1,y=DC2,color=dpt),size=1) +
  scale_color_viridis(option="magma",begin=1,end=0) +
  xlab("Component 1") +
  ylab("Component 2") +
  #ggtitle("Secretory lineage") +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.title=element_blank())

b1 <- pD[pD$branch =="Root",]
pb1 <- ggplot(pD,aes(x=DC1,y=DC2)) +
  geom_point(size=1,color="grey80") +
  geom_point(data=b1,aes(x=DC1,y=DC2,color=dpt),size=1) +
  scale_color_viridis(option="magma",begin=1,end=0) +
  xlab("Component 1") +
  ylab("Component 2") +
  #ggtitle("Hormone-sensing lineage") +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.title=element_blank())


b2 <- pD[pD$branch =="A",]
pb2 <- ggplot(pD,aes(x=DC1,y=DC2)) +
  geom_point(size=1,color="grey80") +
  geom_point(data=b2,aes(x=DC1,y=DC2,color=dpt),size=1) +
  scale_color_viridis(option="magma",begin=1,end=0) +
  xlab("Component 1") +
  ylab("Component 2") +
  #ggtitle("Secretory lineage") +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.title=element_blank())

b3 <- pD[pD$branch == "B",]
pb3 <- ggplot(pD,aes(x=DC1,y=DC2)) +
  geom_point(size=1,color="grey80") +
  geom_point(data=b3,aes(x=DC1,y=DC2,color=dpt), size=1) +
  scale_color_viridis(option="magma",begin=1,end=0) +
  xlab("Component 1") +
  ylab("Component 2") +
  #ggtitle("Secretory lineage") +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.title=element_blank())

branches <- plot_grid(pb0,pb1, pb2,pb3,nrow=2,labels=c("Intermediate","Branch1-Root","Branch2","Branch3"))

