library(Seurat)
library(cowplot)
library(stringr)
library(plyr)
library(dplyr)
source("../src/functions.R")

#--Load Normal ref data----------------- 
Bach<- readRDS("../data/refData/Bach_ExpressionList_QC_norm_clustered_clean.rds")
set.seed(1000)
Bach <- subSample(Bach, cell.number=12000)

m.r1 <- Bach[["counts"]]
pD.r1 <- Bach[["phenoData"]]
fD.r1 <- Bach[["featureData"]]
rm(Bach)

m.r1 <- m.r1[fD.r1$keep, pD.r1$PassAll]
pD.r1 <- pD.r1[pD.r1$PassAll,]
fD.r1 <- fD.r1[fD.r1$keep,]
rownames(m.r1) <- fD.r1$symbol

colnames(m.r1) <- str_replace_all(colnames(m.r1),"-","-Bach-")
pD.r1$barcode <- str_replace_all(pD.r1$barcode, "-","-Bach-")
rownames(pD.r1)<- pD.r1$barcode

#Setup Normal object 
Bach <- CreateSeuratObject(counts = m.r1 , meta.data = pD.r1, project = "Bach", min.cells = 5)
Bach$type <- "Normal"
Bach <- NormalizeData(Bach, verbose = FALSE)
Bach <- FindVariableFeatures(Bach, selection.method = "vst", nfeatures = 2000)

#--Load Tumor Data-----------------
dataList <- readRDS("../data/Robjects/ExpressionList_QC.rds")
m <- dataList[[1]]
pD <- dataList[[2]]
fD <- dataList[[3]]
rm(dataList)

pD$Replicate <- mapvalues(pD$SampleID, c("BRCA1_A","BRCA1_B","NEU_A","NEU_B","FF99WT_A","FF99WT_B","4T1"), 
                          c("A","B","A","B","A","B","A")) %>% 
                          factor(., levels = c("A","B"))
pD$Condition <- pD$SubClusterNumbers <- mapvalues(pD$SampleID, 
                          c("BRCA1_A","BRCA1_B","NEU_A","NEU_B","FF99WT_A","FF99WT_B","4T1"),
                          c("BRCA1_null","BRCA1_null","Neu","Neu","PyMT","PyMT","4T1")) %>% 
                          factor(., levels = c("BRCA1_null","Neu","PyMT","4T1"))
pD$Colors <- mapvalues(pD$SubClusterNumbers, c("BRCA1_null","Neu","PyMT","4T1"),
                       c("#FF00EE","#FF0000","#008000","#0101FA"))

m <- m[, pD$SampleID != "4T1"]
pD <- pD[pD$SampleID !="4T1",]

fD$keep <- rowMeans(m) > 0.01
m <- m[fD$keep, pD$PassAll]
fD <- fD[fD$keep,]
pD <- pD[pD$PassAll,]
rownames(m) <- fD$symbol
rownames(pD) <- pD$barcode

#Set up Tumor object 
Tumors <- CreateSeuratObject(counts = m, meta.data = pD, project = "Tumors", min.cells = 5)
Tumors$type = "Tumor"
Tumors <- NormalizeData(Tumors, verbose = FALSE)
Tumors <- FindVariableFeatures(Tumors, selection.method = "vst", nfeatures = 2000)

#---------------------------
rm(m)
rm(m.r1)
gc()

anchors <- FindIntegrationAnchors(object.list = list(Tumors, Bach), dims = 1:20)
combined <- IntegrateData(anchorset = anchors, dim = 1:20)

#--Run intergrated analysis on all cells------------
DefaultAssay(combined) <- "integrated"

#Run the standard workflow for visualization and clustering 
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)

#--t-SNE and clustering------------
combined <- RunUMAP(combined, reduction = "pca", dims = 1:20)
combined <- RunTSNE(combined, reduction = "pca", dims = 1:20)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)

#--Visualization----------------
p1 <- DimPlot(combined, reduction ="tsne", group.by = "type")
p2 <- DimPlot(combined, reduction = "tsne", label = TRUE)
plot_grid(p1, p2)

DimPlot(combined, reduction = "tsne", split.by = "type",label=TRUE)

combined$SubClusterNumbers <- factor(combined$SubClusterNumbers,
                                     levels = c("BRCA1_null","PyMT","Neu",
                                                "C1","C2","C3","C4","C5",
                                                "C6","C7","C8","C9","C10",
                                                "C11","C12","C13","C14","C15",
                                                "C16","C17","C18","C19","C20"))

Colors <- c("#FF0000","#008000","#0101FA",
            "#9DAFFF","#81C57A","#2A4BD7","#AD2323","#1D6914",
            "#FF9233","#29D0D0","#8126C0","#FFCDF3","#814A19",
            "#FFEE15","#E9DEBB","#A0A0A0","#575757","#000000",
            "#CFFF00","#aaff00","#57FF00","#FF6300","#00A3FE"
            )
names(Colors) <- levels(combined$SubClusterNumbers)

#--Plots-------------------
Colors_Brca <- c("#FF0000",rep("#D3D3D3", each=22))
names(Colors_Brca) <- levels(combined$SubClusterNumbers)
DimPlot(combined, reduction ="tsne", group = "SubClusterNumbers", cols= Colors_Brca,label = F, label.size = 3.5)

Colors_Neu <- c("#D3D3D3","#D3D3D3","#0101FA",rep("#D3D3D3", each=20))
names(Colors_Neu) <- levels(combined$SubClusterNumbers)
DimPlot(combined, reduction ="tsne", group = "SubClusterNumbers", cols= Colors_Neu,label = F, label.size = 5)

Colors_Pymt <- c("#D3D3D3","#008000", rep("#D3D3D3", each=21))
names(Colors_Pymt) <- levels(combined$SubClusterNumbers)
DimPlot(combined, reduction ="tsne", group = "SubClusterNumbers", cols= Colors_Pymt,label = F, label.size = 5)

#--Combined Plots--------
DimPlot(combined, reduction ="tsne", group = "SubClusterNumbers", cols= Colors,label = F, label.size = 5)

genes <- c("Col9a1","Aldh1a3","Mfge8","Spp1","Col4a1","Col1a2","Icam1","Sparcl1","Csn1s2a","Aldoc")
FeaturePlot(combined, features = genes, min.cutoff="q9", reduction="tsne",cols = c("lightgrey","red"))

#---Identify conserved celltype markers-------
DefaultAssay(combined) <- "RNA"
markers <- FindConservedMarkers(combined, ident.1 = 0, grouping.var = "type", verbose = FALSE)
head(markers)
FeaturePlot(combined, features =rownames(markers)[1:6], min.cutoff="q9", reduction="tsne")
FeaturePlot(combined, features = c("Mki67","Ltf","Aldh1a3","Csn1s1"), min.cutoff="q9", reduction="tsne")

#---------------------------------------------
all.genes <- rownames(combined)
combined <- ScaleData(combined, features = all.genes)
all.markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.1, logfc.threshold=1)
all.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

top5 <- all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(combined, features = top5$gene, size = 4) + NoLegend() +theme(axis.text.y = element_text(size = 6))

saveRDS(combined, "../data/Robjects/combinedDataset.rds")
