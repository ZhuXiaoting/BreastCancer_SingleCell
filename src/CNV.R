# install.packages("beanplot")
# install.packages("mixtools")
# install.packages("pheatmap")
# install.packages("zoo")
# install.packages("squash")
# 
# devtools::install_github("diazlab/CONICS/CONICSmat", dep = FALSE)

library(CONICSmat)
library(dplyr)
library(plyr)

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

pD$Condition <- pD$SubClusterNumbers <- mapvalues(pD$SampleID, 
                                                  c("BRCA1_A","BRCA1_B","NEU_A","NEU_B","FF99WT_A","FF99WT_B","4T1"),
                                                  c("BRCA1_null","BRCA1_null","Neu","Neu","PyMT","PyMT","4T1")) %>% 
  factor(., levels = c("BRCA1_null","Neu","PyMT","4T1"))
pD$Colors <- mapvalues(pD$SubClusterNumbers, c("BRCA1_null","Neu","PyMT","4T1"),
                       c("#FF00EE","#FF0000","#008000","#0101FA"))

fD$keep <- rowMeans(m) > 0.01
m <- m[fD$keep, pD$PassAll]
pD <- pD[pD$PassAll,]
fD <- fD[fD$keep, ]
rownames(m) <- fD$symbol
rownames(pD) <- pD$barcode
rownames(fD) <- fD$symbol

calc_cpm <- function(expr_mat) {
  norm_factor <- colSums(expr_mat)
  return (t(t(expr_mat)/norm_factor)*10^6)
}

cpm <- calc_cpm(m)

matrix <- log2(cpm/10 + 1)

regions <- read.table("../chromosome_full_positions_mm10.txt", 
                      sep="\t",row.names = 1,header = T)

gene_pos=getGenePositions(rownames(matrix),species ="mouse")

m=filterMatrix(matrix,gene_pos[,"mgi_symbol"],minCells=5)
normFactor=calcNormFactors(m)

condition = pD$Condition
names(condition) = pD$barcode
sample = pD$SampleID
names(sample) = pD$barcode

l=plotAll(m,normFactor,regions,gene_pos,"test")

par(mar=c(5.1,4.1,4.1,2.1))
hi=plotHistogram(l,m,clusters=2,zscoreThreshold=4,celltypes=sample)

vg=detectVarGenes(m,500)
ts=calculateTsne(m,vg)
plotTsneGene(ts,m,c("MBP","CSF1R","ALDOC","OLIG1"))

lrbic = read.table("CNV_BIC_LR.txt", sep="\t", header = T, row.names = 1, check.names = F)
colnames(lrbic)
candRegions = rownames(lrbic)[which(lrbic[,"BIC difference"]>2000 & lrbic[,"LRT adj. p-val"]<0.01)]

 hi=plotHistogram(l[,candRegions],m,clusters=4,zscoreThreshold=4,sample)

normal= which(hi==1)
tumor=which(hi!=1)
pD$Type <- ""
pD[pD$barcode %in% names(tumor),"Type"] <- "Tumor"
pD[pD$barcode %in% names(normal),"Type"] <- "Normal"

table(pD$SampleID, pD$Type)
table(pD$Condition, pD$Type)

library(ggplot2)

toPlot <- pD %>% group_by(SampleID,Type) %>%
  dplyr::summarise(count = dplyr::n()) %>% 
  mutate(percentage = count/sum(count))

ggplot(toPlot, aes(fill = Type,y= percentage ,x=SampleID))+
  geom_bar(position = "fill", stat= "identity")


redu=plotAll(m,normFactor,regions[candRegions,],gene_pos,"CNVs/Filter1/CNVs_with_info_filter1.pdf",normal=normal,tumor=tumor)

bin_mat=binarizeMatrix(redu,normal,tumor,0.8)
bin_mat[is.na(bin_mat)] <- 0
plotBinaryMat(bin_mat,condition,normal,tumor,"Neu")
plotBinaryMat(bin_mat,condition,normal,tumor,"BRCA1_null")

#-----------
plotBM <- function (mati, sample, normal,  brca, tumor, patient = NULL, k = 3) 
{
  celltypes = rep("Normal", nrow(mati))
  celltypes[tumor] = "Tumor"
  names(celltypes) = rownames(mati)
  patientcolors = data.frame(celltypes = celltypes, BRCA1cluster = brca)
  patientcolors = cbind(patientcolors,sample)
  rownames(patientcolors) = names(celltypes)
  if (!is.null(patient)) {
    p = pheatmap::pheatmap(t(mati[which(sample == patient), 
                                  ]), cluster_cols = T, cutree_cols = 3, annotation = patientcolors, 
                           col = c("lightgrey", "black"), border_color = "grey60", 
                           show_colnames = F, clustering_distance_cols = "euclidean")
  }
  else {
    p = pheatmap::pheatmap(t(mati), cluster_cols = T, cutree_cols = k, 
                           annotation = patientcolors, col = c("lightgrey", 
                                                               "black"), border_color = "grey60", show_colnames = F, 
                           clustering_distance_cols = "euclidean")
  }
  ord = unique(cutree(p$tree_col, k = k)[p$tree_col[["order"]]])
  numb = table(cutree(p$tree_col, k = k))[ord]
  n = length(numb)
  grid::grid.text(expression(bold("Cluster ID \n(left to right)")), 
                  x = rep(0.92), y = c(n * 0.03 + 0.03), 
                  gp = grid::gpar(fontsize = 8, col = "grey"))
  grid::grid.text(ord, x = rep(0.92, length(numb)), 
                  y = seq(n * 0.03, 0.03, -0.03), 
                  gp = grid::gpar(fontsize = 8, col = "grey"))
  return(cutree(p$tree_col, k = k))
}

pD$BRCA1 <- brca@colData[pD$barcode, "Cluster"]
pD$BRCA1 <- mapvalues(pD$BRCA1, c(1,2,3,4,5),c("BRCA1-1","BRCA1-2","BRCA1-3","BRCA1-4","BRCA1-5"))
BRCA1 <- pD$BRCA1
BRCA1[BRCA1 != "BRCA1-3"] <- NA
names(BRCA1) <- rownames(pD)
plotBM(bin_mat,condition,normal,BRCA1, tumor,"BRCA1_null")
#------------


detectBreakPoints (m,normal,tumor,windowsize=1000,gene_pos=gene_pos,chr=1,patients=condition,breakpoints=regions)
detectBreakPoints (m,normal,tumor,windowsize=1000,gene_pos=gene_pos,chr=2,patients=condition,breakpoints=regions)

plotAllChromosomes (mat = m ,normal = normal,tumor = tumor,windowsize = 101,
                    gene_pos = gene_pos,fname = " ",patients = condition,breakpoints = regions)

r=generatePvalMat(m,regions[candRegions,],normFactor,normal,tumor,gene_pos,threshold=0.8)
binr=ifelse(r>0.1,0,1)
boxplot(r)

par(mfrow=c(1,1))
plotChromosomeHeatmap(m,normal = normal, plotcells = which(condition=="BRCA1_null"), 
                      gene_pos = gene_pos, windowsize = 10000, chr=T, expThresh=0.2, thresh = 1)
plotChromosomeHeatmap(m,normal = normal, plotcells = which(condition=="BRCA1_null"), 
                      gene_pos = gene_pos, windowsize = 121, chr=T, expThresh=0.2, thresh = 1)
plotChromosomeHeatmap(m,normal = normal, plotcells = which(condition=="BRCA1_null"), 
                      gene_pos = gene_pos, windowsize = 121, chr=T, expThresh=0.2, thresh = 1)
plotChromosomeHeatmap(m,normal = normal, plotcells = which(condition=="BRCA1_null"), 
                      gene_pos = gene_pos, windowsize = 121, chr=T, expThresh=0.2, thresh = 1)
#--------------------------------

library(Seurat)
library(CaSpER)

library(dplyr)

data <- Read10X(data.dir ='../data/outs/filtered_feature_bc_matrix/')
# counts  <- read.delim("GSE110499_GEO_processed_MM_10X_raw_UMI_count_martix.txt", stringsAsFactor=F, header=T)
# rownames(counts) <- counts[, 1] 
# counts <- counts[, -1]

sc <- CreateSeuratObject(counts = data, project = "bc", min.cells = 3, min.features = 200)

sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
sc <- subset(sc, subset = nFeature_RNA > 500 & percent.mt < 20)

sc <- NormalizeData(sc , scale.factor = 1e6, normalization.method = "RC")
sc <- FindVariableFeatures(sc, do.plot = T, nfeatures = 1000)
sc <- ScaleData(sc)

sc <- RunPCA(sc, features = VariableFeatures(object = sc),npcs = 100)
sc <- RunTSNE(sc, dims.use = 1:10)

sc <- FindNeighbors(sc, dims = 1:10)
sc <- FindClusters(sc, resolution = 0.5)

log.ge <- as.matrix(sc@assays$RNA@data)
control <- names(Idents(sc) )[Idents(sc) %in% c(2,7)]
mm <- names(Idents(sc) )[Idents(sc) %in% c(0, 1, 3, 4)]

genes <- rownames(log.ge)
annotation <- generateAnnotation(id_type="hgnc_symbol", genes=genes, centromere=centromere, ishg19 = T)
log.ge <- log.ge[match( annotation$Gene,rownames(log.ge)) , ]
rownames(log.ge) <- annotation$Gene
log.ge <- log2(log.ge +1)