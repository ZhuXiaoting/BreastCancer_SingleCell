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

#--------------------------------------------
set.seed(1000)
m_FF99WT <- m[,pD$barcode %in% FF99WT]
pD_FF99WT <- pD[pD$barcode %in% FF99WT,]

#=======================================================
fD_A <- fD %>% dplyr::mutate(keep = rowMeans(m_FF99WT[,pD_FF99WT$Replicate == "A"])> 0.01)
fD_B <- fD %>% dplyr::mutate(keep = rowMeans(m_FF99WT[,pD_FF99WT$Replicate == "B"])> 0.01)

sce.FF99WT_A <- SingleCellExperiment(list(counts=as.matrix(m_FF99WT[fD_A$keep,pD_FF99WT$Replicate == "A"])),
                                  colData = DataFrame(pD_FF99WT[pD_FF99WT$Replicate == "A",]),
                                  rowData = DataFrame(fD_A[fD_A$keep,]))

sce.FF99WT_B <- SingleCellExperiment(list(counts=as.matrix(m_FF99WT[fD_B$keep,pD_FF99WT$Replicate == "B"])),
                                  colData = DataFrame(pD_FF99WT[pD_FF99WT$Replicate == "B",]),
                                  rowData = DataFrame(fD_B[fD_B$keep,]))
#-------------------------
clusters <- quickCluster(sce.FF99WT_A, method ='igraph',use.ranks=FALSE, min.mean = 0.1)
table(clusters)
sce.FF99WT_A <- computeSumFactors(sce.FF99WT_A, min.mean = 0.1, clusters = clusters)
summary(sizeFactors(sce.FF99WT_A))
sce.FF99WT_A <- normalize(sce.FF99WT_A)
#table(sce.NEU$Replicate)
#batch <- c(rep("1", each=2242), rep("2", each = 1912))
fit <- trendVar(sce.FF99WT_A, use.spikes = FALSE)
dec <- decomposeVar(sce.FF99WT_A, fit)
dec$Symbol <- rowData(sce.FF99WT_A)$symbol
dec_A <- dec[order(dec$bio, decreasing = TRUE), ]
hvg.A <- dec_A[which(dec_A$FDR <= 0.1 & dec_A$bio >=0),]

set.seed(1000)
sce.FF99WT_A <- runPCA(sce.FF99WT_A, feature_set= rownames(hvg.A),BSPARAM=IrlbaParam())
sce.FF99WT_A  <- runTSNE(sce.FF99WT_A , use_dimred = "PCA")
sce.FF99WT_A  <- runUMAP(sce.FF99WT_A , use_dimred = "PCA")

#plotTSNE(sce.FF99WT_A, colour_by="Cd14")

snn.gr <- buildSNNGraph(sce.FF99WT_A, use.dimred="PCA", assay.type="logcounts", k=200)
clusters <- igraph::cluster_louvain(snn.gr)
table(clusters$membership)

sce.FF99WT_A$Cluster <- factor(clusters$membership)
plotTSNE(sce.FF99WT_A, colour_by="Cluster")
markers <- findMarkers(sce.FF99WT_A, sce.FF99WT_A$Cluster, direction="up")

library(xlsx)
for (cluster in names(markers)) {
  write.xlsx(data.frame(markers[[cluster]]), row.names=TRUE, 
             file="FF99WT_A_MarkerGenes.xlsx",
             sheetName=cluster,  append = TRUE)
  gc()
}

#================================================
clusters <- quickCluster(sce.FF99WT_B, method ='igraph',use.ranks=FALSE, min.mean = 0.1)
table(clusters)
sce.FF99WT_B <- computeSumFactors(sce.FF99WT_B, min.mean = 0.1, clusters = clusters)
summary(sizeFactors(sce.FF99WT_B))
sce.FF99WT_B <- normalize(sce.FF99WT_B)
#table(sce.NEU$Replicate)
#batch <- c(rep("1", each=2242), rep("2", each = 1912))
fit <- trendVar(sce.FF99WT_B, use.spikes = FALSE)
dec <- decomposeVar(sce.FF99WT_B, fit)
dec$Symbol <- rowData(sce.FF99WT_B)$symbol
dec_B <- dec[order(dec$bio, decreasing = TRUE), ]

hvg.B <- dec_B[which(dec_B$FDR <= 0.1 & dec_B$bio >=0),]

set.seed(1000)
sce.FF99WT_B <- runPCA(sce.FF99WT_B, feature_set= rownames(hvg.B),BSPARAM=IrlbaParam())
sce.FF99WT_B  <- runTSNE(sce.FF99WT_B, use_dimred = "PCA")
sce.FF99WT_B <- runUMAP(sce.FF99WT_B, use_dimred = "PCA")

snn.gr <- buildSNNGraph(sce.FF99WT_B, use.dimred="PCA", assay.type="logcounts", k=200)
clusters <- igraph::cluster_louvain(snn.gr)
table(clusters$membership)

sce.FF99WT_B$Cluster <- factor(clusters$membership)
plotTSNE(sce.FF99WT_B, colour_by="Cluster")
plotTSNE(sce.FF99WT_B, colour_by="Cd14")

markers <- findMarkers(sce.FF99WT_B, sce.FF99WT_B$Cluster, direction="up")  

for (cluster in names(markers)) {
  write.xlsx(data.frame(markers[[cluster]]), row.names=TRUE, 
             file="FF99WT_B_MarkerGenes.xlsx",
             sheetName=cluster,  append = TRUE)
  gc()
}

#==========================================
genes <- c("Cd14","Bcl3","Osmr","Nfkbia")

plt <- list()
for (gene in genes){
  tmp = plotTSNE(sce.FF99WT_A, colour_by = gene )+scale_fill_gradient2(name = gene, low='grey',high ='red')
  plt[[gene]] <- tmp
}
multiplot(plotlist=plt, cols = 2)
plotTSNE(sce.FF99WT_A, colour_by ="Cluster")

#==================================================================================================

universe <- intersect(rownames(dec_A), rownames(dec_B))
mean.bio <- (dec_A[universe,"bio"] + dec_B[universe,"bio"])/2
chosen <- universe[mean.bio > 0]
length(chosen)

rescaled <- batchelor::multiBatchNorm(
  sce.FF99WT_A[universe,], 
  sce.FF99WT_B[universe,]
)
rescaled.FF99WT_A <- rescaled[[1]]
rescaled.FF99WT_B <- rescaled[[2]]

set.seed(1000) 
unc.FF99WT_A <- logcounts(rescaled.FF99WT_A)[chosen,]
unc.FF99WT_B <- logcounts(rescaled.FF99WT_B)[chosen,]

mnn.out <- batchelor::fastMNN(
  FF99WT_A=unc.FF99WT_A, FF99WT_B=unc.FF99WT_B,
  k=20, d=50, BSPARAM=IrlbaParam(deferred=TRUE)
)
mnn.out

sce.FF99WT <- mnn.out
sce.FF99WT <- runTSNE(sce.FF99WT, use_dimred="corrected")
plotTSNE(sce.FF99WT, colour_by="batch") + ggtitle("Corrected")

sce.FF99WT <- runUMAP(sce.FF99WT, use_dimred = "corrected")

assay(sce.FF99WT, "original") <- cbind(unc.FF99WT_A, unc.FF99WT_B)
osce.FF99WT <- runPCA(sce.FF99WT, exprs_values = "original",ntop = Inf, BSPARAM=IrlbaParam())
osce.FF99WT <- runTSNE(osce.FF99WT, use_dimred = "PCA")
plotTSNE(osce.FF99WT, colour_by = "batch") +ggtitle("Original")

# rowData(sce.NEU)$symbol <- fD$symbol[fD$id %in%rownames(sce.NEU)]
# rownames(sce.NEU) <- rowData(sce.NEU)$symbol
# plotTSNE(sce.NEU, by_exprs_values="reconstructed", colour_by="Cd14")
# genes <- c("Cd14","Bcl3","Osmr","Nfkbia")
# 
# plt <- list()
# for (gene in genes){
#   tmp = plotTSNE(sce.NEU, by_exprs_values="original", colour_by = gene )+scale_fill_gradient2(name = gene, low='grey',high ='red')
#   plt[[gene]] <- tmp
# }
# multiplot(plotlist=plt, cols = 2)

#-----------------------
snn.gr <- buildSNNGraph(sce.FF99WT, use.dimred="corrected", k=40)  # k100 3clusters k50 4clusters
clusters <- igraph::cluster_louvain(snn.gr)
table(clusters$membership, sce.FF99WT$batch)

sce.FF99WT$Cluster <- factor(clusters$membership)

plotTSNE(sce.FF99WT, colour_by="Cluster")

#-----------------------

markers <- findMarkers(sce.FF99WT, sce.FF99WT$Cluster,block = sce.FF99WT$batch,
                       assay.type="original", direction="up")  

for (cluster in names(markers)) {
  write.xlsx(data.frame(markers[[cluster]]), row.names=TRUE, 
             file="FF99WT_combine_5clusters_MarkerGenes.xlsx",
             sheetName=cluster,  append = TRUE)
  gc()
}
