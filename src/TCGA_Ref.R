#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)

query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  experimental.strategy = "RNA-Seq",
                  platform = "Illumina HiSeq",
                  file.type = "results",
                  legacy = TRUE)
GDCdownload(query, method = "api")
BRCARnaseq_assay <- GDCprepare(query)
saveRDS(BRCARnaseq_assay, "TCGA-BRCA-RNASeq.rds")

table(BRCARnaseq_assay$subtype_BRCA_Subtype_PAM50)

BRCAMatrix <- assay(BRCARnaseq_assay,"scaled_estimate")
BRCAMatrix_TPM <- BRCAMatrix * 1000000

symbol <- rownames(BRCAMatrix_TPM) %>% strsplit("[|]") 
symbol <- data.frame(matrix(unlist(symbol), nrow=length(symbol), byrow=T),stringsAsFactors=FALSE)
  
rownames(BRCAMatrix_TPM) <- symbol$X1
BRCAMatrix_TPM <- BRCAMatrix_TPM[, !is.na(BRCARnaseq_assay$subtype_BRCA_Subtype_PAM50)]

name = "ref_TCGA_PAM50"
expr = as.matrix(log2(BRCAMatrix_TPM+1))
types = as.character(BRCARnaseq_assay$subtype_BRCA_Subtype_PAM50[!is.na(BRCARnaseq_assay$subtype_BRCA_Subtype_PAM50)])
main_types = types 
ref = list(name = name, data = expr, types= types, main_types = main_types)
save(ref,file='../data/ref_TCGA_PAM50.RData')

#Homolog
homolog <- read.delim("../data/homologene.data.txt", header=F)
names(homolog) <- c("HID","Taxon", "EntrezID","Symbol","ProteinGI","ProteinAcc")

Human <- homolog[homolog$Taxon =="9606",]
Mouse <- homolog[homolog$Taxon == "10090",]
homoGene <- merge(Human, Mouse, by.x="HID",by.y="HID",suffixes = c(".human", ".mouse"))

match <- match(rownames(BRCAMatrix_TPM), homoGene$Symbol.human)
rownames(BRCAMatrix_TPM) <- homoGene$Symbol.mouse[match]

BRCAMatrix_TPM <- BRCAMatrix_TPM[!is.na(rownames(BRCAMatrix_TPM)),]

name = "ref_TCGA_PAM50_homolog"
expr = as.matrix(log2(BRCAMatrix_TPM+1))
types = as.character(BRCARnaseq_assay$subtype_BRCA_Subtype_PAM50[!is.na(BRCARnaseq_assay$subtype_BRCA_Subtype_PAM50)])
main_types = types 
ref = list(name = name, data = expr, types= types, main_types = main_types)
save(ref,file='../data/ref_TCGA_PAM50_homolog.RData')

