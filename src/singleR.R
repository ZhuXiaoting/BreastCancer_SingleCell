library(Seurat)
library(scran)
library(dplyr)
library(plyr)
library(SingleR)
source("../src/functions.R")

#======load data ================
dataList <- readRDS("../data/Robjects/ExpressionList_QC_pairNorm.rds")
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

fD$keep <- rowMeans(m) > 0.01
m <- m[fD$keep, pD$PassAll]
fD <- fD[fD$keep,]
pD <- pD[pD$PassAll,]

#======make reference=============
# Ref1
Bach <- readRDS("../data/Bach_ExpressionList_QC_norm_clustered_clean.rds")
set.seed(1000)
Bach <- subSample(Bach, cell.number=12000)

m.r1 <- Bach[["counts"]]
pD.r1 <- Bach[["phenoData"]]
fD.r1 <- Bach[["featureData"]]
rm(Bach)

m.r1 <- m.r1[fD.r1$keep, pD.r1$PassAll]
pD.r1 <- pD.r1[pD.r1$PassAll,]
fD.r1 <- fD.r1[fD.r1$keep,]

m.r1 <- t(t(m.r1)/pD.r1$sf)
m.r1 <- log2(m.r1+1)

rownames(m.r1) <- fD.r1$symbol

name = "ref1_Bach"
expr = as.matrix(m.r1)
types = as.character(pD.r1$SubClusterNumbers)
main_type = as.character(pD.r1$SubCluster)
ref = list(name = name, data = expr, types= types, main_types = main_type)
save(ref,file='../data/ref1_Bach_v2.RData')

#=================
#Ref2 

#=================
#Ref3 
gunzip("GSE111113_Table_S1_FilterNormal10xExpMatrix.txt.gz")
Giraddi <- read.table("GSE111113_Table_S1_FilterNormal10xExpMatrix.txt", sep="\t", header = T)

library(stringr)
p_G <-  as.character(colnames(Giraddi)[-c(1:3)]) %>%
  str_replace(., "__","_") %>%
  str_split(., "\\_|\\__", simplify = TRUE) %>%
  as.data.frame()%>%
  dplyr::rename(Sample = V1, barcode=V2)
 
f_G <- Giraddi[,1:3]%>% 
  dplyr::rename(gene_name = X, transcript_id = transcript_id.s.)

m_G <- Giraddi[,-c(1:3)]
rownames(m_G) <- make.names(Giraddi$X, unique=TRUE)

Datalist <- list("phenoData" = p_G, "featureData" = f_G, "counts" = m_G)
saveRDS(Datalist, file = "../data/Giraddi_Expression.rds")

m_G <- log2(m_G+1)
name = "ref3_Giraddi"
expr = as.matrix(m_G)
types = as.character(p_G$Sample)
main_types = mapvalues(p_G$Sample, c("Adu1","Adu2","AduBas", "E16", "E18","P4"), 
                      c("Adu","Adu","Adu","E16","E18","P4"))
ref3 = list(name = name, data = expr, types= types, main_types = main_types)
save(ref3,file='../data/ref3_Giraddi.RData')

#=======================================================
# calculate correlation with singleR

rownames(m) <- fD$symbol
annot <- pD$SampleID
names(annot) <- pD$barcode

load("ref.RData")
load("../data/ref2_Pal_AdultEpithelial.RData")

singler = CreateSinglerSeuratObject(m, annot =annot, project.name = "Ref-TCGA", 
                                    min.genes = 200, technology = "10X", species="Mouse",
                                    citation="TCGA-BRCA", ref.list = list(ref),
                                    normalize.gene.length=F,
                                    variable.genes = "de", fine.tune = F, reduce.file.size = T, 
                                    do.signatures = F, min.cells = 2, npca = 10, 
                                    regress.out = "nUMI", do.main.types = T, 
                                    reduce.seurat.object = T, numCores = SingleR.numCores)
#Reference Bach and default reference

singler$meta.data$orig.ID <- pD$Condition

# out = SingleR.PlotTsne(singler$singler[[1]]$SingleR.single,
#                                                singler$meta.data$xy, do.label = FALSE, do.letters = T,
#                                                labels=singler$meta.data$orig.ID,label.size = 6, 
#                                                dot.size = 1, alpha = 0.5)
# out$p


SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single, top.n = Inf,
                                         clusters = singler$meta.data$orig.ID)

SingleR.DrawHeatmap(singler3$singler[[1]]$SingleR.single, top.n = Inf,
                                          clusters = singler3$meta.data$orig.ID)


library(knitr)
kable(table(singler$singler[[1]]$SingleR.single$labels, singler$meta.data$orig.ID))

save(singler2,file="SingleR_combine_ref1_Bach_v2.RData")
