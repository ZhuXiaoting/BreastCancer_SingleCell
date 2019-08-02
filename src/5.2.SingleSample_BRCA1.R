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

set.seed(1000)
m_BRCA <- m[,pD$barcode %in% BRCA1]
pD_BRCA <- pD[pD$barcode %in% BRCA1,]

#=======================================================
fD_A <- fD %>% dplyr::mutate(keep = rowMeans(m_BRCA[,pD_BRCA$Replicate == "A"])> 0.01)
fD_B <- fD %>% dplyr::mutate(keep = rowMeans(m_BRCA[,pD_BRCA$Replicate == "B"])> 0.01)

sce.BRCA_A <- SingleCellExperiment(list(counts=as.matrix(m_BRCA[fD_A$keep,pD_BRCA$Replicate == "A"])),
                                  colData = DataFrame(pD_BRCA[pD_BRCA$Replicate == "A",]),
                                  rowData = DataFrame(fD_A[fD_A$keep,]))

sce.BRCA_B <- SingleCellExperiment(list(counts=as.matrix(m_BRCA[fD_B$keep,pD_BRCA$Replicate == "B"])),
                                  colData = DataFrame(pD_BRCA[pD_BRCA$Replicate == "B",]),
                                  rowData = DataFrame(fD_B[fD_B$keep,]))
#-------------------------
clusters <- quickCluster(sce.BRCA_A, method ='igraph',use.ranks=FALSE, min.mean = 0.1)
table(clusters)
sce.BRCA_A <- computeSumFactors(sce.BRCA_A, min.mean = 0.1, clusters = clusters)
summary(sizeFactors(sce.BRCA_A))
sce.BRCA_A <- normalize(sce.BRCA_A)
#table(sce.NEU$Replicate)
#batch <- c(rep("1", each=2242), rep("2", each = 1912))
fit <- trendVar(sce.BRCA_A, use.spikes = FALSE)
dec <- decomposeVar(sce.BRCA_A, fit)
dec$Symbol <- rowData(sce.BRCA_A)$symbol
dec_A <- dec[order(dec$bio, decreasing = TRUE), ]
hvg.A <- dec_A[which(dec_A$FDR <= 0.1 & dec_A$bio >=0),]

set.seed(1000)
sce.BRCA_A <- runPCA(sce.BRCA_A, feature_set= rownames(hvg.A),BSPARAM=IrlbaParam())
sce.BRCA_A  <- runTSNE(sce.BRCA_A , use_dimred = "PCA")
sce.BRCA_A  <- runUMAP(sce.BRCA_A , use_dimred = "PCA")

plotTSNE(sce.BRCA_A, colour_by="Cd14")

snn.gr <- buildSNNGraph(sce.BRCA_A, use.dimred="PCA", assay.type="logcounts", k= 50)
clusters <- igraph::cluster_louvain(snn.gr)
table(clusters$membership)

sce.BRCA_A$Cluster <- factor(clusters$membership)
plotTSNE(sce.BRCA_A, colour_by="Cluster")
markers <- findMarkers(sce.BRCA_A, sce.BRCA_A$Cluster, direction="up")

library(xlsx)
for (cluster in names(markers)) {
  write.xlsx(data.frame(markers[[cluster]]), row.names=TRUE, 
             file="BRCA_A_MarkerGenes.xlsx",
             sheetName=cluster,  append = TRUE)
  gc()
}

#================================================
clusters <- quickCluster(sce.BRCA_B, method ='igraph',use.ranks=FALSE, min.mean = 0.1)
table(clusters)
sce.BRCA_B <- computeSumFactors(sce.BRCA_B, min.mean = 0.1, clusters = clusters)
summary(sizeFactors(sce.BRCA_B))
sce.BRCA_B <- normalize(sce.BRCA_B)
#table(sce.NEU$Replicate)
#batch <- c(rep("1", each=2242), rep("2", each = 1912))
fit <- trendVar(sce.BRCA_B, use.spikes = FALSE)
dec <- decomposeVar(sce.BRCA_B, fit)
dec$Symbol <- rowData(sce.BRCA_B)$symbol
dec_B <- dec[order(dec$bio, decreasing = TRUE), ]

hvg.B <- dec_B[which(dec_B$FDR <= 0.1 & dec_B$bio >=0),]

set.seed(1000)
sce.BRCA_B <- runPCA(sce.BRCA_B, feature_set= rownames(hvg.B),BSPARAM=IrlbaParam())
sce.BRCA_B  <- runTSNE(sce.BRCA_B, use_dimred = "PCA")
sce.BRCA_B <- runUMAP(sce.BRCA_B, use_dimred = "PCA")

snn.gr <- buildSNNGraph(sce.BRCA_B, use.dimred="PCA", assay.type="logcounts", k=200)
clusters <- igraph::cluster_louvain(snn.gr)
table(clusters$membership)

sce.BRCA_B$Cluster <- factor(clusters$membership)
plotTSNE(sce.BRCA_B, colour_by="Cluster")
plotTSNE(sce.BRCA_B, colour_by="Cd14")

markers <- findMarkers(sce.BRCA_B, sce.BRCA_B$Cluster, direction="up")  

for (cluster in names(markers)) {
  write.xlsx(data.frame(markers[[cluster]]), row.names=TRUE, 
             file="BRCA_B_MarkerGenes.xlsx",
             sheetName=cluster,  append = TRUE)
  gc()
}

#==========================================
genes <- c("Cd14","Bcl3","Osmr","Nfkbia")

plt <- list()
for (gene in genes){
  tmp = plotTSNE(sce.BRCA_A, colour_by = gene )+scale_fill_gradient2(name = gene, low='grey',high ='red')
  plt[[gene]] <- tmp
}
multiplot(plotlist=plt, cols = 2)
plotTSNE(sce.BRCA_A, colour_by ="Cluster")

#==================================================================================================

universe <- intersect(rownames(dec_A), rownames(dec_B))
mean.bio <- (dec_A[universe,"bio"] + dec_B[universe,"bio"])/2
chosen <- universe[mean.bio > 0]
length(chosen)

rescaled <- batchelor::multiBatchNorm(
  sce.BRCA_A[universe,], 
  sce.BRCA_B[universe,]
)
rescaled.BRCA_A <- rescaled[[1]]
rescaled.BRCA_B <- rescaled[[2]]

set.seed(1000) 
unc.BRCA_A <- logcounts(rescaled.BRCA_A)[chosen,]
unc.BRCA_B <- logcounts(rescaled.BRCA_B)[chosen,]

mnn.out <- batchelor::fastMNN(
  BRCA_A=unc.BRCA_A, BRCA_B=unc.BRCA_B,
  k=20, d=50, BSPARAM=IrlbaParam(deferred=TRUE)
)
mnn.out

sce.BRCA <- mnn.out
sce.BRCA <- runTSNE(sce.BRCA, use_dimred="corrected")
plotTSNE(sce.BRCA, colour_by="batch") + ggtitle("Corrected")

sce.BRCA <- runUMAP(sce.BRCA, use_dimred = "corrected")

assay(sce.BRCA, "original") <- cbind(unc.BRCA_A, unc.BRCA_B)
osce.BRCA <- runPCA(sce.BRCA, exprs_values = "original",ntop = Inf, BSPARAM=IrlbaParam())
osce.BRCA <- runTSNE(osce.BRCA, use_dimred = "PCA")
plotTSNE(osce.BRCA, colour_by = "batch") +ggtitle("Original")

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
snn.gr <- buildSNNGraph(sce.BRCA, use.dimred="corrected", k=50) 
clusters <- igraph::cluster_louvain(snn.gr)
table(clusters$membership, sce.BRCA$batch)

sce.BRCA$Cluster <- factor(clusters$membership)

plotTSNE(sce.BRCA, colour_by="Cluster")

# plotTSNE(sce.NEU, by_exprs_values="original", colour_by="Cd14")+
#   scale_fill_gradient2(name = gene, low='grey',high ='red')

#-----------------------

markers <- findMarkers(sce.BRCA, sce.BRCA$Cluster,block = sce.BRCA$batch,
                       assay.type="original", direction="up")  

for (cluster in names(markers)) {
  write.xlsx(data.frame(markers[[cluster]]), row.names=TRUE, 
             file="BRCA_combine_5clusters_MarkerGenes.xlsx",
             sheetName=cluster,  append = TRUE)
  gc()
}

sceList <- list("BRCA_A" = sce.BRCA_A, "BRCA1_B" = sce.BRCA_B, "combine" = sce.BRCA)
saveRDS(sceList, "../data/Robjects/BRCA1_sceList.rds")


#==============================================================================
fD <- fD %>% dplyr::mutate(keep = rowMeans(m_BRCA)> 0.01)

sce.BRCA <- SingleCellExperiment(list(counts=as.matrix(m_BRCA[fD$keep,])),
                                   colData = DataFrame(pD_BRCA),
                                   rowData = DataFrame(fD[fD$keep,]))

clusters <- quickCluster(sce.BRCA, method ='igraph',use.ranks=FALSE, min.mean = 0.1)
table(clusters)
sce.BRCA<- computeSumFactors(sce.BRCA, min.mean = 0.1, clusters = clusters)
summary(sizeFactors(sce.BRCA))
sce.BRCA <- normalize(sce.BRCA)
#table(sce.NEU$Replicate)
batch <- c(rep("1", each=240), rep("2", each =2044))
fit <- trendVar(sce.BRCA,block = batch, use.spikes = FALSE)
dec <- decomposeVar(sce.BRCA, fit)
dec$Symbol <- rowData(sce.BRCA)$symbol
dec <- dec[order(dec$bio, decreasing = TRUE), ]

hvg <- dec[which(dec$FDR <= 0.1 & dec$bio >=0),]

set.seed(1000)
sce.BRCA <- runPCA(sce.BRCA, feature_set= rownames(hvg),BSPARAM=IrlbaParam())
sce.BRCA  <- runTSNE(sce.BRCA, use_dimred = "PCA")
sce.BRCA <- runUMAP(sce.BRCA, use_dimred = "PCA")

plotTSNE(sce.BRCA, colour_by="Replicate")

snn.gr <- buildSNNGraph(sce.BRCA, use.dimred="PCA", assay.type="logcounts", k=200)
clusters <- igraph::cluster_louvain(snn.gr)
table(clusters$membership)

sce.BRCA$Cluster <- factor(clusters$membership)
plotTSNE(sce.BRCA, colour_by="Cluster")

