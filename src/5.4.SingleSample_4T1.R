library(scran)
library(dplyr)
library(Rtsne)
library(plyr)
library(BiocSingular)
library(scater)
library(batchelor)

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
m <- m[fD$keep, pD$PassAll]
pD <- pD[pD$PassAll,]
fD <- fD[fD$keep, ]
rownames(m) <- fD$symbol
rownames(pD) <- pD$barcode
rownames(fD) <- fD$symbol
#-------------------------------------------
BRCA1 <- pD$barcode[pD$Condition %in% c("BRCA1")]
FF99WT <- pD$barcode[pD$Condition %in% c("FF99WT")]
NEU <- pD$barcode[pD$Condition %in% c("NEU")]
Sample_4T1 <- pD$barcode[pD$Condition %in% c("4T1")]

#----------------------------------------------------------------------------------------------------------
#4T1

m_4T1 <- m[, pD$barcode %in% Sample_4T1]
pD_4T1 <- pD[pD$barcode %in% Sample_4T1,]

fD_4T1 <- fD %>% dplyr::mutate(keep = rowMeans(m_4T1)>0.01)
m_4T1 <- m_4T1[fD_4T1$keep, ]
fD_4T1 <- fD_4T1[fD_4T1$keep,]

rownames(pD_4T1) <- pD_4T1$barcode

sce.4T1 <- SingleCellExperiment(list(counts=as.matrix(m_4T1)),
                                colData = DataFrame(pD_4T1),
                                rowData = DataFrame(fD_4T1))

#=========================
clusters <- quickCluster(sce.4T1, method ='igraph',use.ranks=FALSE, min.mean = 0.1)
table(clusters)
sce.4T1 <- computeSumFactors(sce.4T1, min.mean = 0.1, clusters = clusters)
summary(sizeFactors(sce.4T1))
sce.4T1 <- normalize(sce.4T1)
#table(sce.NEU$Replicate)
#batch <- c(rep("1", each=2242), rep("2", each = 1912))
fit <- trendVar(sce.4T1, use.spikes = FALSE)
dec <- decomposeVar(sce.4T1, fit)
dec$Symbol <- rowData(sce.4T1)$symbol
dec_A <- dec[order(dec$bio, decreasing = TRUE), ]
hvg.A <- dec_A[which(dec_A$FDR <= 0.1 & dec_A$bio >=0),]

set.seed(1000)
sce.4T1 <- runPCA(sce.4T1, feature_set= rownames(hvg.A),BSPARAM=IrlbaParam())
sce.4T1  <- runTSNE(sce.4T1 , use_dimred = "PCA")
sce.4T1  <- runUMAP(sce.4T1 , use_dimred = "PCA")

plotTSNE(sce.4T1, colour_by="Cd14")

snn.gr <- buildSNNGraph(sce.4T1, use.dimred="PCA", assay.type="logcounts", k=300)
clusters <- igraph::cluster_louvain(snn.gr)
table(clusters$membership)

sce.4T1$Cluster <- factor(clusters$membership)
plotTSNE(sce.4T1, colour_by="Cluster")
markers <- findMarkers(sce.4T1, sce.4T1$Cluster, direction="up")

library(xlsx)
for (cluster in names(markers)) {
  write.xlsx(data.frame(markers[[cluster]]), row.names=TRUE, 
             file="4T1__MarkerGenes.xlsx",
             sheetName=cluster,  append = TRUE)
  gc()
}