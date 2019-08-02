library(plyr)
library(dplyr)
library(ggplot2)
library(monocle)
library(RColorBrewer)
library(scran)
source("../src/functions.R")

# Load Data
dataList <- readRDS("../data/Robjects/ExpressionList_QC_norm.rds")
 m <- dataList[[1]]
 pD <- dataList[[2]]
 fD <- dataList[[3]]
 rm(dataList)
 gc()

#-- Prepocessing --------------------
fD$keep <- rowMeans(m) > 0.01
m <- m[fD$keep,pD$PassAll]
pD <- pD[pD$PassAll,]
m <- m[,pD$SampleID !="4T1"]
pD <- pD[pD$SampleID !="4T1",]
fD <- fD[fD$keep,]

#--Normalize counts by size factor---
m.norm <- t(t(m))/pD$sf
rownames(m.norm) <- fD$symbol
fD$gene_short_name <- fD$symbol

#------------------
pD.af <- new("AnnotatedDataFrame", data = pD)
sampleNames(pD.af) <- pD$barcode
fD.af <- new("AnnotatedDataFrame", data=fD)
sampleNames(fD.af) <- fD$symbol

cds <- newCellDataSet(as.matrix(m.norm),
                      phenoData=pD.af,
                      featureData=fD.af,
                      lowerDetectionLimit=1,
                      expressionFamily=negbinomial.size())

#--Trajectory inference according to vignette
phenoData(cds)$Size_Factor <- pD$sf
cds <- estimateDispersions(cds)

#--Feature Selection---------------------
genes <- fD$id[fD$highVar]
cds <- setOrderingFilter(cds,genes)

#--Dim Reduction-------------------------
cds <- reduceDimension(cds, max_components=2,norm_method="log")
cds <- orderCells(cds,reverse=FALSE)

#--Plot trajectory colored by clusters----

p0 <- plot_cell_trajectory(cds,x=1,y=2, color_by="Condition", cell_size=1,show_branch_points=F)+
  theme(legend.position = "right")+
  scale_color_manual(values=levels(pD$Colors))+
  facet_wrap(~Condition, nrow=1)

p1 <- plot_cell_trajectory(cds,x=1,y=2, color_by="Pseudotime", cell_size=1,show_branch_points=F)+
  theme(legend.position = "right")+
  #scale_color_manual(values=levels(pD$Colors))+
  facet_wrap(~Condition, nrow=1)

#--Plot genes expression-------------------
blast_genes <- row.names(subset(fData(cds), symbol %in% c("Cd14","Krt4")))
cds_subset <- cds[blast_genes,]

plot_genes_branched_pseudotime(cds_subset, color_by = "Condition", branch_point = 1)+
  scale_color_manual(values = levels(pD$Colors))

#--save------------------------------------
pD.monoc <- pData(cds)
monoc <- list("pD"=pD.monoc,
              "plot"=p0)
saveRDS(monoc,"../data/Robjects/ExpressionList_Monocle.rds")

#----------------------------------------------
# return significance score for each gene. Genes score significant are said to be branch-dependent in their expression 
# BEAMres <- BEAM(cds, branch_point = 1, cores = 1)

# BEAMres <- BEAMres[order(BEAMres$qval),]
# BEAMres <- BEAMres[,c("symbol", "pval", "qval")]
# plot_genes_branched_heatmap(cds[row.names(subset(BEAMres,qval < 1e-4)),],
#                             branch_point = 1,
#                             num_clusters = 4,
#                             cores = 1,
#                             use_gene_short_name = T,
#                             show_rownames = T)

# # Plot Genes
# plot_genes <- row.names(subset(fData(cds),
#                                gene_short_name %in% c("Ccnd2", "Sftpb", "Pdpn")))
# plot_genes_branched_pseudotime(cds[plot_genes,],
#                                branch_point = 1,
#                                color_by = "Condition",
#                                ncol = 1)
