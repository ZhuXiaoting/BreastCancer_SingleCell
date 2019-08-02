library(scran)
library(dplyr)
library(Rtsne)
library(plyr)
library(BiocSingular)
library(scater)
library(batchelor)
source("../src/functions.R")

dataList <- readRDS("../data/Robjects/ExpressionList_QC.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]
rm(dataList)

pD$Replicate <- mapvalues(pD$SampleID, c("BRCA1_A","BRCA1_B","NEU_A","NEU_B","FF99WT_A","FF99WT_B","4T1"), 
                          c("A","B","A","B","A","B","A")) %>% factor(., levels = c("A","B"))
pD$Condition <- mapvalues(pD$SampleID, c("BRCA1_A","BRCA1_B","NEU_A","NEU_B","FF99WT_A","FF99WT_B","4T1"),
                          c("BRCA1","BRCA1","NEU","NEU","FF99WT","FF99WT","4T1")) %>% 
  factor(., levels = c("BRCA1","NEU","FF99WT","4T1"))

fD$keep <- rowMeans(m) > 0.01

m <- m[, pD$PassAll]
pD <- pD[pD$PassAll,]
#---------------------------
sce <- SingleCellExperiment(list(counts=as.matrix(m)),
                            colData = DataFrame(pD),
                            rowData = DataFrame(fD))
set.seed(1000)
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
assignments <- cyclone(sce, pairs=mm.pairs)
table(assignments$phases)
pD$phases <- assignments$phases
sce$phases <- assignments$phases

m <- m[fD$keep,]
sce <- sce[fD$keep,]
fD <- fD[fD$keep,]

saveRDS(sce, "../data/Robjects/AllSamples_raw_sce.rds")

#--Normalize-------------------
sce <- sce[,pD$SampleID !="4T1"]

clusters <- quickCluster(sce, method ='igraph', min.mean = 0.1)
table(clusters)

sce <- computeSumFactors(sce, min.mean = 0.1, clusters = clusters)
summary(sizeFactors(sce))
pD$sf <- sizeFactors(sce)

sce <- normalize(sce)

#--High variant genes--------
fit <- trendVar(sce, use.spikes = FALSE)
dec <- decomposeVar(sce, fit)

dec$Symbol <- rowData(sce)$Symbol
dec <- dec[order(dec$bio, decreasing = TRUE), ]
head(dec)
hvg.out <- dec[which(dec$FDR <= 0.1 & dec$bio >=0),]

#----------------------------
fD$highVar <- fD$id %in% rownames(hvg.out)
fD$highVarFDR <- dec$FDR
fD$highVarBiolComp <- dec$bio

#--Dimension reduction-------
set.seed(1000)
osce <- runPCA(sce, feature_set= rownames(hvg.out),BSPARAM=IrlbaParam())
osce <- runTSNE(osce, use_dimred = "PCA")
osce <- runUMAP(osce, use_dimred = "PCA")
plotTSNE(osce, colour_by="SampleID") 
plotUMAP(osce, colour_by = "SampleID")

#--Clustering--------------
snn.gr <- buildSNNGraph(osce, use.dimred="PCA", assay.type="logcounts", k=250)
clusters <- igraph::cluster_walktrap(snn.gr)
osce$Cluster <- factor(clusters$membership)
table(osce$Cluster)

plotTSNE(osce, colour_by = "Cluster")
plotUMAP(osce, colour_by = "Cluster")

#--Plot-------------------

pD <- pD[pD$SampleID != "4T1",]
pD$Cluster <- osce$Cluster

pD$tSNE1 <- reducedDim(osce,"TSNE")[,1]
pD$tSNE2 <- reducedDim(osce,"TSNE")[,2]

pD$UMAP1 <- reducedDim(osce,"UMAP")[,1]
pD$UMAP2 <- reducedDim(osce,"UMAP")[,2]

pD$TumorType <- mapvalues(pD$Condition, c("BRCA1","NEU","FF99WT"), c("BRCA1-null","Neu","PyMT"))%>%
  factor(.,levels = c("BRCA1-null","PyMT","Neu"))

pD$Colors <- mapvalues(pD$TumorType, c("BRCA1-null","PyMT","Neu"),
                       c("#FF0000","#008000","#0101FA")) 

#--UMAP--------------------
P1 <- ggplot(pD, aes(x = UMAP1, y = UMAP2, color = Replicate))+
  geom_point(size=1)+
  theme_void(base_size = 12)+
  xlab("UMAP1")+
  ylab("UMAP2")+
  ggtitle("UMAP")+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  theme_bw()

P2<- ggplot(pD, aes(x = UMAP1, y = UMAP2, color = TumorType))+
  geom_point(size=1)+
  theme_void(base_size = 12)+
  scale_color_manual(values=levels(pD$Colors))+
  xlab("UMAP1")+
  ylab("UMAP2")+
  ggtitle("UMAP")+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  theme_bw()

P3<- ggplot(pD, aes(x = UMAP1, y = UMAP2, color = phases))+
  geom_point(size=1)+
  theme_void(base_size = 12)+
  xlab("UMAP1")+
  ylab("UMAP2")+
  ggtitle("UMAP")+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  theme_bw()

P4<- ggplot(pD, aes(x = UMAP1, y = UMAP2, color = Cluster))+
  geom_point(size=1)+
  theme_void(base_size = 12)+
  xlab("UMAP1")+
  ylab("UMAP2")+
  ggtitle("UMAP")+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  theme_bw()
#--TSNE-----------------------
P5 <- ggplot(pD, aes(x = tSNE1, y = tSNE2, color = Replicate))+
  geom_point(size=1)+
  theme_void(base_size = 12)+
  xlab("tSNE1")+
  ylab("tSNE2")+
  ggtitle("tSNE")+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  theme_bw()
  
P6 <- ggplot(pD, aes(x = tSNE1, y = tSNE2, color = TumorType))+
  geom_point(size=1)+
  theme_void(base_size = 12)+
  scale_color_manual(values=levels(pD$Colors))+
  xlab("tSNE1")+
  ylab("tSNE2")+
  ggtitle("tSNE")+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  theme_bw()

P7 <- ggplot(pD, aes(x = tSNE1, y = tSNE2, color = phases))+
  geom_point(size=1)+
  theme_void(base_size = 12)+
  xlab("tSNE1")+
  ylab("tSNE2")+
  ggtitle("tSNE")+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  theme_bw()

P8 <- ggplot(pD, aes(x = tSNE1, y = tSNE2, color = Cluster))+
  geom_point(size=1)+
  theme_void(base_size = 12)+
  xlab("tSNE1")+
  ylab("tSNE2")+
  ggtitle("tSNE")+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  theme_bw()

saveRDS(osce, "../Data/Robjects/3Tumors_norm_sce.rds")

#---save object--------------
pD.add <- pD[,c("barcode","Replicate","tSNE1","tSNE2","UMAP1","UMAP2", "phases","Cluster","TumorType","Colors","sf")]
fD.add <- fD[,c("id","highVar","highVarFDR","highVarBiolComp")]

dataList <- readRDS("../data/Robjects/ExpressionList_QC.rds")
dataList[["phenoData"]] <- left_join(dataList[["phenoData"]],pD.add)
fD.new <- left_join(dataList[["featureData"]],fD.add)
fD.new$highVar[is.na(fD.new$highVar)] <- FALSE
dataList[["featureData"]] <- fD.new

saveRDS(dataList,file="../data/Robjects/ExpressionList_QC_norm.rds")

#--Find markers------------------
rownames(osce) <- rowData(osce)$symbol
markers <- findMarkers(osce, osce$Cluster, direction="up")

library(openxlsx)
wb <- createWorkbook("../data/Combined_Marker_Genes.xlsx")
for (cluster in names(markers)) {
  addWorksheet(wb,sheetName = cluster)
  writeData(wb, data.frame(markers[[cluster]]), rowNames=TRUE, 
             sheet=cluster)
}
saveWorkbook(wb, "../data/Combined_markerGenes.xlsx")

#--Plot Heatmap of markers----------
library(pheatmap)

top.markers_c <- c()
for (i in 1:6){
  marker.tmp <- markers[[i]]
  top.tmp <- rownames(marker.tmp)[marker.tmp$Top <=10]
  top.markers_c <- c(top.markers_c, top.tmp)
}
top.markers_c <- top.markers_c[!duplicated(top.markers_c)]

top.exprs <- logcounts(osce)[top.markers_c,,drop=FALSE]
heat.vals <- top.exprs - rowMeans(top.exprs)

heat.vals <- heat.vals[,order(factor(pD$Cluster, levels = c(4,5,1,2,3,6)))]

##----------------
forcol <- ggplot_build(P8)
condColors <- unique(dplyr::arrange(forcol$data[[1]],group) %>% .[["colour"]])

col = list(Tumor = levels(pD$Colors), Cluster = condColors)
names(col$Tumor) <- levels(pD$TumorType)
names(col$Cluster) <- levels(pD$Cluster)

#-----------------
pheatmap(heat.vals, 
         cluster_cols=FALSE,
         cluster_rows = TRUE,
         show_colnames = FALSE,
         annotation_colors = col,
         annotation_col=data.frame(Cluster = factor(pD$Cluster),
                                   Tumor=factor(pD$TumorType),
                                   row.names=colnames(osce)),
         gaps_col = c(1608,2309,4268,6483,8942),
         fontsize_row = 5.5)

#------------------------------------------------------------------
markers <- findMarkers(osce, osce$Condition, direction="up")
top.markers_c <- c()
for (i in 1:3){
  marker.tmp <- markers[[i]]
  top.tmp <- rownames(marker.tmp)[marker.tmp$Top <=10]
  top.markers_c <- c(top.markers_c, top.tmp)
}
top.markers_c <- top.markers_c[!duplicated(top.markers_c)]

top.exprs <- logcounts(osce)[top.markers_c,,drop=FALSE]
heat.vals <- top.exprs - rowMeans(top.exprs)

pheatmap(heat.vals, 
         cluster_cols=T,
         cluster_rows = TRUE,
         show_colnames = FALSE,
         annotation_colors = col[1], 
         annotation_col=data.frame(Tumor=factor(pD$TumorType),
                                   row.names=colnames(osce)),
         gaps_col = c(1608,2309,4268,6483,8942),
         fontsize_row = 6.5)

#--Plot Gene expression--------
#sce <- readRDS("../data/Robjects/3Tumors_norm_sce.rds")
#rownames(sce) <- rowData(sce)$symbol

genes <- c("Krt14","Aldh1a3","Csn1s1", "Vim","Ltf","Lalba","Sparc", "Spp1", "Tspan8")
#========================
plt <- list()
for (gene in genes){
  tmp = plotTSNE(sce, colour_by = gene )+scale_fill_gradient2(name = gene, low='grey',high ='red')
  plt[[gene]] <- tmp
}
multiplot(plotlist=plt, cols = 3)

#----------------------------------------------------------------
# c1 <- c("Cited1","Prlr","Esr1","Areg")
# c2 <- c("Rspo1","Atp6v1c2","Fabp3","Thrsp","Wap","Glycam1","Olah")
# c3 <- c("Foxa1","Ly6a","Aldh1a3","Kit","Cd14")
# c5 <- c("Lypd3")
# c6 <- c("1500015O10Rik","Col7a1","Moxd1","Mia","Emid1","Pdpn","Col9a2","Fbln2","Igfbp3","Fst","Il17b")
# c7 <- c("Oxtr","Krt15","Igfbp6","Igfbp2","Tns1")
# c9 <- c("Gng11","Procr","Igfbp7","Notch3","Zeb2")
#=========================
# genes <- c(c1,c3,c5,c2,c6,c7,c9)
# mheat<- logcounts(osce)[genes, ]
# mheat <- mheat/apply(mheat,1,max) 

# annotCol <- data.frame("Cluster" = osce$Cluster,
#                        #"Batch" = osce$Replicate,
#                        "Condition" = osce$Condition)
# rownames(annotCol) <- as.character(osce$barcode)

# Color1<- c("red", "green","blue")
# names(Color1) <- levels(sce$Condition)[1:3]

# Color2 <- c("yellow","brown")
# names(Color2) <- levels(sce$Replicate)
# annotColors <- list("Batch" = Color2)

# mheat
# mh <- mheat[,order(osce$Cluster)]
# p <-  pheatmap(mheat,
#                cluster_rows=FALSE,
#                cluster_cols=FALSE,
#                show_rownames=TRUE,
#                show_colnames=FALSE,
#                annotation_legend=TRUE,
#                annotation_col=annotCol,
#                # gaps_col=c(263,463,663,863,1052,1252,1352),
#                annotation_colors=annotColors,
#                fontsize=9)

# genes <- c("Krt18","Krt8","Krt5", "Krt14","Acta2","Vim","Aldh1a3","Ltf","Wnt5b","Csn1s1","Lalba","Cd36")
# genes <- c("Prlr","Ly6a","Cited1")
# genes <- c("Igfbp2","Wif1","Igfbp4","Id4","Cxcl14","Pdpn")
# genes <- c("Cd14","Bcl3","Osmr","Nfkbia")
# genes <- c("Mki67","Ltf","Aldh1a3", "Csn1s1")





