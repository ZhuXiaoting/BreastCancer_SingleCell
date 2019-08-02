# Script to prepare cellranger data for downstream analysis

library(Matrix)
library(dplyr)
library(plyr)
# ---- ReadData ----

# Read in output from cell ranger using cellrangerRkit
# The latest version of loading cellranger matrices into R

matrix_dir = "../data/outs/filtered_feature_bc_matrix/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1
cDat<-as.matrix(mat)
pDat<-data.frame(colnames(mat))
colnames(pDat)<-c("barcode")
rownames(pDat)<-pDat$barcode
fDat<-data.frame(feature.names[,1:2])
colnames(fDat)<-c("id","symbol")
rownames(fDat)<-fDat$id

#reduce size of matrix
keep <- rowSums(cDat) > 1
cDat <- cDat[keep,]
fDat <- data.frame(fDat[rownames(cDat),])

# ---- Formatting ----

# Add more info to phenotype Data
pDat <- mutate(pDat, SeqID=substr(barcode,18,19)) %>%
	mutate(barcode=as.character(barcode)) %>%
        mutate(SampleID=mapvalues(SeqID,c("1","2","3","4","5","6","7"),
				 c("BRCA1_A","BRCA1_B","NEU_A","NEU_B","FF99WT_A","FF99WT_B","4T1")))
# pDat$Condition<-as.factor(pDat$Condition)


# Add more info to the feature Data
mitoGenes <- read.table("../data/refData/MitoGenes.txt")
tfCheck <- read.table("../data/refData/TFcheckpoint_WithENSID.tsv",
		header=TRUE, sep="\t")

fDat$Mitochondrial <- fDat$id %in% mitoGenes$V1
fDat$TranscriptionFactor <- fDat$id %in% tfCheck$ensembl_gene_id


# Save data
stopifnot(identical(rownames(fDat),rownames(cDat)) & identical(colnames(cDat),pDat$barcode))
DataList <- list("phenoData"=pDat, "featureData"=fDat, "counts"=cDat)
saveRDS(DataList,file="../data/Robjects/ExpressionList.rds")
