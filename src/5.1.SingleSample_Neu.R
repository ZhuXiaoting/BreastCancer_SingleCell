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
# NEU
set.seed(1000)
m_NEU <- m[,pD$barcode %in% NEU]
pD_NEU <- pD[pD$barcode %in% NEU,]

#=======================================================
fD_A <- fD %>% dplyr::mutate(keep = rowMeans(m_NEU[,pD_NEU$Replicate == "A"])> 0.01)
fD_B <- fD %>% dplyr::mutate(keep = rowMeans(m_NEU[,pD_NEU$Replicate == "B"])> 0.01)
  
sce.NEU_A <- SingleCellExperiment(list(counts=as.matrix(m_NEU[fD_A$keep,pD_NEU$Replicate == "A"])),
                                colData = DataFrame(pD_NEU[pD_NEU$Replicate == "A",]),
                                rowData = DataFrame(fD_A[fD_A$keep,]))

sce.NEU_B <- SingleCellExperiment(list(counts=as.matrix(m_NEU[fD_B$keep,pD_NEU$Replicate == "B"])),
                                  colData = DataFrame(pD_NEU[pD_NEU$Replicate == "B",]),
                                  rowData = DataFrame(fD_B[fD_B$keep,]))
#-------------------------
clusters <- quickCluster(sce.NEU_A, method ='igraph',use.ranks=FALSE, min.mean = 0.1)
table(clusters)
sce.NEU_A <- computeSumFactors(sce.NEU_A, min.mean = 0.1, clusters = clusters)
summary(sizeFactors(sce.NEU_A))
sce.NEU_A <- normalize(sce.NEU_A)
#table(sce.NEU$Replicate)
#batch <- c(rep("1", each=2242), rep("2", each = 1912))
fit <- trendVar(sce.NEU_A, use.spikes = FALSE)
dec <- decomposeVar(sce.NEU_A, fit)
dec$Symbol <- rowData(sce.NEU_A)$symbol
dec_A <- dec[order(dec$bio, decreasing = TRUE), ]
hvg.A <- dec_A[which(dec_A$FDR <= 0.1 & dec_A$bio >=0),]

set.seed(1000)
sce.NEU_A <- runPCA(sce.NEU_A, feature_set= rownames(hvg.A),BSPARAM=IrlbaParam())
sce.NEU_A  <- runTSNE(sce.NEU_A , use_dimred = "PCA")
sce.NEU_A  <- runUMAP(sce.NEU_A , use_dimred = "PCA")
#rownames(sce.NEU_A) <- rowData(sce.NEU_A)$symbol
plotTSNE(sce.NEU_A, colour_by="Cd14")

snn.gr <- buildSNNGraph(sce.NEU_A, use.dimred="PCA", assay.type="logcounts", k=500)
clusters <- igraph::cluster_louvain(snn.gr)
table(clusters$membership)

sce.NEU_A$Cluster <- factor(clusters$membership)
plotTSNE(sce.NEU_A, colour_by="Cluster")
markers <- findMarkers(sce.NEU_A, sce.NEU_A$Cluster, direction="up")

library(xlsx)
for (cluster in names(markers)) {
  write.xlsx(data.frame(markers[[cluster]]), row.names=TRUE, 
             file="NEU_A__MarkerGenes.xlsx",
             sheetName=cluster,  append = TRUE)
  gc()
}

#----------------------------
clusters <- quickCluster(sce.NEU_B, method ='igraph',use.ranks=FALSE, min.mean = 0.1)
table(clusters)
sce.NEU_B <- computeSumFactors(sce.NEU_B, min.mean = 0.1, clusters = clusters)
summary(sizeFactors(sce.NEU_B))
sce.NEU_B <- normalize(sce.NEU_B)
#table(sce.NEU$Replicate)
#batch <- c(rep("1", each=2242), rep("2", each = 1912))
fit <- trendVar(sce.NEU_B, use.spikes = FALSE)
dec <- decomposeVar(sce.NEU_B, fit)
dec$Symbol <- rowData(sce.NEU_B)$symbol
dec_B <- dec[order(dec$bio, decreasing = TRUE), ]

hvg.B <- dec_B[which(dec_B$FDR <= 0.1 & dec_B$bio >=0.1),]

sce.NEU_B <- runPCA(sce.NEU_B, feature_set= rownames(hvg.B),BSPARAM=IrlbaParam())
sce.NEU_B  <- runTSNE(sce.NEU_B, use_dimred = "PCA")
sce.NEU_B <- runUMAP(sce.NEU_B, use_dimred = "PCA")
rownames(sce.NEU_B) <- rowData(sce.NEU_B)$symbol

snn.gr <- buildSNNGraph(sce.NEU_B, use.dimred="PCA", assay.type="logcounts", k=300)
clusters <- igraph::cluster_walktrap(snn.gr)
table(clusters$membership)

sce.NEU_B$Cluster <- factor(clusters$membership)
plotTSNE(sce.NEU_B, colour_by="Cluster")
plotTSNE(sce.NEU_B, colour_by="Cd14")

markers <- findMarkers(sce.NEU_B, sce.NEU_B$Cluster, direction="up")  

for (cluster in names(markers)) {
  write.xlsx(data.frame(markers[[cluster]]), row.names=TRUE, 
             file="NEU_B_MarkerGenes.xlsx",
             sheetName=cluster,  append = TRUE)
  gc()
}

#--Plot Genes expression--------------
genes <- c("Cd14","Bcl3","Osmr","Nfkbia")

plt <- list()
for (gene in genes){
  tmp = plotTSNE(sce.NEU_A, colour_by = gene )+
  scale_fill_gradient2(name = gene, low='grey',high ='red')+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  plt[[gene]] <- tmp
}
multiplot(plotlist=plt, cols = 2)
plotTSNE(sce.NEU_A, colour_by ="Cluster")

#==CombineDataset=============================================================

universe <- intersect(rownames(dec_A), rownames(dec_B))
mean.bio <- (dec_A[universe,"bio"] + dec_B[universe,"bio"])/2
chosen <- universe[mean.bio > 0]
length(chosen)

rescaled <- batchelor::multiBatchNorm(
  sce.NEU_A[universe,], 
  sce.NEU_B[universe,]
)
rescaled.NEU_A <- rescaled[[1]]
rescaled.NEU_B <- rescaled[[2]]

set.seed(1000) 
unc.NEU_A <- logcounts(rescaled.NEU_A)[chosen,]
unc.NEU_B <- logcounts(rescaled.NEU_B)[chosen,]

mnn.out <- batchelor::fastMNN(
  NEU_A=unc.NEU_A, NEU_B=unc.NEU_B,
  k=20, d=50, BSPARAM=IrlbaParam(deferred=TRUE)
)
mnn.out

sce.NEU <- mnn.out
sce.NEU <- runTSNE(sce.NEU, use_dimred="corrected")
plotTSNE(sce.NEU, colour_by="batch") + ggtitle("Corrected")

assay(sce.NEU, "original") <- cbind(unc.NEU_A, unc.NEU_B)
osce.NEU <- runPCA(sce.NEU, exprs_values = "original",ntop = Inf, BSPARAM=IrlbaParam())
osce.NEU <- runTSNE(osce.NEU, use_dimred = "PCA")
plotTSNE(osce.NEU, colour_by = "batch") +ggtitle("Original")
#osce.NEU <- runUMAP(osce.NEU, use_dimred = "PCA")
#plotUMAP(osce.NEU, colour_by = "batch")

sceList <- list("sce_A" = sce.NEU_A, "sce_B" = sce.NEU_B, "combine"=sce.NEU)
saveRDS(sceList, "../data/Robjects/Neu_sceList.rds")

#--Plot Gene expression-----------------
genes <- c("Cd14","Bcl3","Osmr","Nfkbia","Aqp5","Kcnn4","Col9a1","Apoc1")

plt <- list()
for (gene in genes){
  tmp = plotTSNE(sce.NEU, by_exprs_values="original", colour_by = gene )+scale_fill_gradient2(name = gene, low='grey',high ='red')
  plt[[gene]] <- tmp
}
multiplot(plotlist=plt, cols = 2)

#--Clustering------------------------------------
snn.gr <- buildSNNGraph(sce.NEU, use.dimred="corrected", k=100) #k=100 3 cluster k=50: 4 cluster k=40 : 5clusters
clusters <- igraph::cluster_louvain(snn.gr)
table(clusters$membership, sce.NEU$batch)

sce.NEU$Cluster <- factor(clusters$membership)
plotTSNE(sce.NEU, colour_by="Cluster")

#-----------------------
markers <- findMarkers(sce.NEU, sce.NEU$Cluster,block = sce.NEU$batch,
                       assay.type="original", direction="up")  

for (cluster in names(markers)) {
  write.xlsx(data.frame(markers[[cluster]]), row.names=TRUE, 
             file="NEU_combine_3clusters_MarkerGenes.xlsx",
             sheetName=cluster,  append = TRUE)
  gc()
}

#===============================================================================
# NEU_A <- CreateSeuratObject(counts =m_NEU[,pD_NEU$Replicate == "A"] , project = "A", min.cells = 5)
# NEU_A$Rep <- "A"
# NEU_A <- NormalizeData(NEU_A, verbose = FALSE)
# NEU_A <- FindVariableFeatures(NEU_A, selection.method = "vst", nfeatures = 2000)

# # Set up stimulated object
# NEU_B <- CreateSeuratObject(counts = m_NEU[,pD_NEU$Replicate == "B"], project = "B", min.cells = 5)
# NEU_B$Rep <- "B"
# NEU_B <- NormalizeData(NEU_B, verbose = FALSE)
# NEU_B <- FindVariableFeatures(NEU_B, selection.method = "vst", nfeatures = 2000)

# anchors <- FindIntegrationAnchors(object.list = list(NEU_A, NEU_B), dims = 1:20)
# Neu <- IntegrateData(anchorset = anchors, dims = 1:20)

# DefaultAssay(Neu) <- "integrated"

# Neu <- ScaleData(Neu, verbose = FALSE)
# Neu <- RunPCA(Neu, npcs = 30, verbose=FALSE)
# Neu <- RunUMAP(Neu, reduction = "pca", dims = 1:20)
# Neu <- RunTSNE(Neu, reduction = "pca",dims = 1:20)

# DimPlot(Neu, reduction = "tsne", group.by = "Rep")

# Neu <- FindNeighbors(Neu, reduction ="pca", dims = 1:20)
# Neu <- FindClusters(Neu, resolution = 0.2)

# DimPlot(Neu, reduction = "tsne", label=T)

# DefaultAssay(Neu) <- "RNA"
# markers <- FindAllMarkers(Neu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)



# FeaturePlot(Neu, features = c("Cd14")) 


# head(dec)
# hvg.out_A <- dec[which(dec$FDR <= 0.1 & dec$bio >=0),]

# sce.NEU_A <- runPCA(sce.NEU,feature_set= rownames(hvg.out), BSPARAM=IrlbaParam())
# sce.NEU <- runTSNE(sce.NEU, use_dimred = "PCA")
# plotTSNE(sce.NEU, colour_by = "SampleID")
# sce.NEU <- runUMAP(sce.NEU, use_dimred ="PCA")
# plotUMAP(sce.NEU, colour_by = "SampleID")

# snn.gr <- buildSNNGraph(sce.NEU, use.dimred="PCA")
# clusters <- igraph::cluster_fast_greedy(snn.gr)
# sce.NEU$Cluster <- factor(clusters$membership)
# plotTSNE(sce.NEU, colour_by = "Cluster")

# rownames(sce) <- rowData(sce)$symbol
# sceList[[name]] <- sce

# saveRDS(sceList, "SingleSamples_norm_sce.rds")


# genes <- c("Tslp", "Ctla2a", 'H2-Ab1', "Il24", "Pdgfrl", 'Ly6a', "H2-Aa", "Slpi", "H2-Eb1")

# plotTSNE(sceList[["NEU"]], colour_by ="Cluster")
# plotUMAP(sce, colour_by ="Cluster")

# rownames(sce) <- rowData(sce)$symbol


# genes <- c("Mgp","Kctd1","Matn4","Cp","Mt1","Steap4","Selenop","Mt2")

# plt <- list()
# for (gene in genes){
#   tmp = plotTSNE(sce.BRCA1, colour_by = gene )+scale_fill_gradient2(name = gene, low='grey',high ='red')
#   plt[[gene]] <- tmp
# }
# multiplot(plotlist=plt, cols = 4)
# #=======================================================================

# fD_NEU <- fD %>% dplyr::mutate(keep = rowMeans(m_NEU) > 0.01)
# sce.NEU <- SingleCellExperiment(list(counts=as.matrix(m_NEU[fD_NEU$keep,])),
#                                   colData = DataFrame(pD_NEU),
#                                   rowData = DataFrame(fD_NEU[fD_NEU$keep,]))

# clusters <- quickCluster(sce.NEU, method ='igraph',use.ranks=FALSE, min.mean = 0.1)
# table(clusters)
# sce.NEU <- computeSumFactors(sce.NEU, min.mean = 0.1, clusters = clusters)
# summary(sizeFactors(sce.NEU))
# sce.NEU <- normalize(sce.NEU)
# table(sce.NEU$Replicate)
# batch <- c(rep("1", each=2242), rep("2", each = 1912))
# fit <- trendVar(sce.NEU, use.spikes = FALSE, block = batch)
# dec <- decomposeVar(sce.NEU, fit)
# dec$Symbol <- rowData(sce.NEU)$symbol
# dec <- dec[order(dec$bio, decreasing = TRUE), ]
# hvg <- dec[which(dec_A$FDR <= 0.1 & dec_A$bio >=0),]

# sce.NEU <- runPCA(sce.NEU, BSPARAM=IrlbaParam())
# sce.NEU <- runTSNE(sce.NEU, use_dimred = "PCA")
# plotTSNE(sce.NEU, colour_by = "SampleID")


# plotTSNE(sce.BRCA1, colour_by="Cluster")
# plotUMAP(sce.4T1, colour_by="Cluster")


# snn.gr <- buildSNNGraph(sce.BRCA1, use.dimred="PCA", k = 50)
# clusters <- igraph::cluster_fast_greedy(snn.gr)
# table(clusters$membership)
# sce.BRCA1$Cluster <- factor(clusters$membership)

# jgc <- function()
# {
#   .jcall("java/lang/System", method = "gc")
# } 
  
# options(java.parameters = "-Xmx8g")
# markers <- findMarkers(sce.BRCA1, sce.BRCA1$Cluster, direction="up")  
# library(xlsx)
# for (cluster in names(markers)) {
#   write.xlsx(data.frame(markers[[cluster]]), row.names=TRUE, 
#              file="BRCA1_MarkerGenes.xlsx",
#              sheetName=cluster,  append = TRUE)
#   jgc()
# }


# sce.4T1 <- sceList[["4T1"]]
# plotTSNE(sce.4T1, colour_by="Cluster")
# plotUMAP(sce.4T1, colour_by="Cluster")

# markers <- findMarkers(sce.4T1, sce.4T1$Cluster, direction="up")  

# for (cluster in names(markers)) {
#   write.xlsx(data.frame(markers[[cluster]]), row.names=TRUE, 
#              file="4T1_MarkerGenes.xlsx",
#              sheetName=cluster,  append = TRUE)
#   gc()
# }


# sce.NEU <- sceList[["NEU"]]
# plotTSNE(sce.NEU, colour_by="Cluster")
# plotUMAP(sce.NEU, colour_by="Cluster")

# markers <- findMarkers(sce.NEU, sce.NEU$Cluster, direction="up")  

# for (cluster in names(markers)) {
#   write.xlsx(data.frame(markers[[cluster]]), row.names=TRUE, 
#              file="Neu_MarkerGenes.xlsx",
#              sheetName=cluster,  append = TRUE)
#   gc()
# }


# sce.FF99WT <- sceList[["FF99WT"]]
# plotTSNE(sce.FF99WT, colour_by="Cluster")

# snn.gr <- buildSNNGraph(sce.FF99WT, use.dimred="PCA", k = 50)
# clusters <- igraph::cluster_fast_greedy(snn.gr)
# table(clusters$membership)
# sce.FF99WT$Cluster <- factor(clusters$membership)

# plotTSNE(sce.FF99WT, colour_by="Cluster")
# plotUMAP(sce.FF99WT, colour_by="Cluster")

# markers <- findMarkers(sce.FF99WT, sce.FF99WT$Cluster, direction="up")  

# for (cluster in names(markers)) {
#   write.xlsx(data.frame(markers[[cluster]]), row.names=TRUE, 
#              file="FF99WT_MarkerGenes.xlsx",
#              sheetName=cluster,  append = TRUE)
#   gc()
# }



# #====================================================================


# top.markers_c <- c()
# for (i in 1:3){
#   marker.tmp <- markers[[i]]
#   top.tmp <- rownames(marker.tmp)[marker.tmp$Top <=10]
#   top.markers_c <- c(top.markers_c, top.tmp)
# }
# top.markers_c <- top.markers_c[!duplicated(top.markers_c)]

# top.exprs <- logcounts(sce.BRCA1)[top.markers_c,,drop=FALSE]
# heat.vals <- top.exprs - rowMeans(top.exprs)
# pheatmap(heat.vals, cluster_cols=TRUE,
#          show_colnames = FALSE,
#          annotation_col=data.frame(Cluster = factor(sce.BRCA1$Cluster),
#                                    row.names=colnames(sce.BRCA1)),
#          fontsize_row = 7)

  
# #-------------------------------------------
# universe <- intersect(rownames(decList[[1]]), rownames(decList[[2]]))

# mean.bio <- (decList[[1]][universe, "bio"] + decList[[2]][universe, "bio"])/2
# chosen <- universe[mean.bio >0] 
# length(chosen)

# rescaled <- batchelor::multiBatchNorm(sceList[[1]][universe, ], sceList[[2]][universe,])

# #-------------------------------------------

# repA <- logcounts(rescaled[[1]])[chosen,]
# repB <- logcounts(rescaled[[2]])[chosen,]


# mnn.out <- batchelor::fastMNN(repA = repA, repB=repB, k = 20, d = 50, BSPARAM=IrlbaParam(deferred=TRUE))
# dim(reducedDim(mnn.out,"corrected"))

# Rle(mnn.out$batch)

# metadata(mnn.out)$merge.info$pairs[[1]]
# #-------------------------------------
# sce <- mnn.out

# assay(sce,"original") <- cbind(repA, repB)
# sce <- SingleCellExpriment(list(logcounts = omat))
# reducedDim(sce, "MNN") <- mnn.out$corrected
# sce$Batch <- as.character(mnn.out$batch)
# sce
# #---------------------------------------------
# osce <- runPCA(sce, exprs_values="original",ntop = Inf, BSPARAM=IrlbaParam())
# osce <- runTSNE(osce, use_dimred = "PCA")
# ot <- plotTSNE(osce, colour_by="batch") + ggtitle("Original")

# csce <- runTSNE(sce, use_dimred = "corrected")
# ct <- plotTSNE(csce, colour_by = "batch") + ggtitle("Corrected")

# multiplot(ot, ct, cols = 2)


# metadata(mnn.out)$merge.info$lost.var

