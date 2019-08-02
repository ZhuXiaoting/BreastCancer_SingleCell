library(scran)
library(dplyr)
library(plyr)


#--get homolog genes between human and mouse------
homolog <- read.delim("../data/homologene.data.txt", header=F)
names(homolog) <- c("HID","Taxon", "EntrezID","Symbol","ProteinGI","ProteinAcc")

Human <- homolog[homolog$Taxon =="9606",]
Mouse <- homolog[homolog$Taxon == "10090",]
homoGene <- merge(Human, Mouse, by.x="HID",by.y="HID",suffixes = c(".human", ".mouse"))
#write.table(homoGene, "results/Human_Mouse_homolog.txt", sep="\t")

#-------------------------------------------------- 
library(genefu)
library(xtable)
library(rmeta)
library(Biobase)
library(caret)

set.seed(1000)
sce <- readRDS("../data/Robjects/3Tumors_norm_sce.rds")
scdata <- counts(sce)
rownames(scdata) <- rowData(sce)$symbol

sf <- sizeFactors(sce)
scdata <- t(t(scdata)/sf)

match <- match(rownames(scdata), homoGene$Symbol.mouse)
rownames(scdata) <- homoGene$EntrezID.human[match]

scdata <- scdata[!is.na(rownames(scdata)),]
scdata <- scdata[!duplicated(rownames(scdata)),]

annot <- as.data.frame(rownames(scdata))
names(annot) <- "EntrezGene.ID"
annot$probe <- annot$EntrezGene.ID
rownames(annot) <- as.character(annot$EntrezGene.ID)

#-------------------------------
claudinLowPreds<-molecular.subtyping(sbt.model = "claudinLow",
                                     data = t(scdata),
                                     annot = annot,do.mapping = T)
#[1] "Number of genes used: 649"

table(claudinLowPreds$subtype)

#-------------------------------
PAM50Preds<-molecular.subtyping(sbt.model = "pam50",data=t(scdata),
                                annot=annot, do.mapping = T)
table(PAM50Preds$subtype) 

#------------------------------
pD <- colData(sce)
pD$Claudin <- claudinLowPreds$subtype
pD$pam50 <- as.character(PAM50Preds$subtype)
index <- pD$Claudin =="Claudin"
pD$pam50[index] <- "ClaudinLow"
pD$pam50 <- factor(pD$pam50)

table(pD$pam50,pD$SampleID)

#-----------------------------------------------------------------------
sceList <- readRDS("SingleSamples_norm_sce.rds")
sce.4T1 <- sceList[["4T1"]]
rm(sceList)

sf <- sizeFactors(sce.4T1)
m.4T1 <- counts(sce.4T1)
m.4T1 <- t(t(m.4T1)/sf)
rownames(m.4T1) <- rowData(sce.4T1)$symbol

scdata.4T1 <- m.4T1

match <- match(rownames(scdata.4T1), homoGene$Symbol.mouse)
rownames(scdata.4T1) <- homoGene$EntrezID.human[match]

scdata.4T1 <- scdata.4T1[!is.na(rownames(scdata.4T1)),]
scdata.4T1 <- scdata.4T1[!duplicated(rownames(scdata.4T1)),]

annot.4T1 <- as.data.frame(rownames(scdata.4T1))
names(annot.4T1) <- "EntrezGene.ID"
annot.4T1$probe <- annot.4T1$EntrezGene.ID
rownames(annot.4T1) <- as.character(annot.4T1$EntrezGene.ID)

#------------------------------------
PAM50Preds.4T1<-molecular.subtyping(sbt.model = "pam50",data=t(scdata.4T1),
                                annot=annot.4T1, do.mapping = T)
table(PAM50Preds.4T1$subtype) 

claudinLowPreds.4T1<-molecular.subtyping(sbt.model = "claudinLow",data=t(scdata.4T1),
                                    annot=annot.4T1, do.mapping = T)
table(claudinLowPreds.4T1$subtype)

#---------------------------------------------------------------------
res <- data.frame(pam50=pD$pam50, Condition=pD$Condition, ClaudinLow = pD$claudinLow)
res.4T1 <- data.frame(pam50 = PAM50Preds.4T1$subtype, Condition = sce.4T1$Condition, 
                      ClaudinLow = claudinLowPreds.4T1$subtype)

res <- rbind(res, res.4T1)

res$pam50 <- as.character(res$pam50)
res[res$ClaudinLow =="Claudin",]$pam50 <- "ClaudinLow"
res$pam50 <- factor(res$pam50, level = c("ClaudinLow","Basal","Normal",
                                             "LumA","LumB","Her2"))

ggplot(res, mapping = aes(x = Condition, group= pam50))+
  geom_bar(mapping = aes(y = ..prop..,fill = factor(..x..)), stat = "count")+
  scale_y_continuous(labels = scales::percent)+
  labs(y = "% of cells", x = "TumorType", fill ="Tumor Type")+
  scale_fill_discrete(labels = c("BRCA1-null","Neu","PyMT","4T1"))+
  scale_x_discrete(labels =  c("BRCA1-null","Neu","PyMT","4T1"))+
  theme_classic()+
  facet_wrap(~pam50)


ggplot(res, mapping = aes(x = variable, y = value,group = pam50,fill = variable))+
  geom_bar(stat = "identity")+
  scale_y_continuous(labels = scales::percent)+
  scale_fill_discrete(labels = c("BRCA1-null","Neu","PyMT","4T1"))+
  scale_x_discrete(labels =  c("BRCA1-null","Neu","PyMT","4T1"))+
  labs(y = "% of cells", x = "TumorType", fill ="Tumor Type")+
  theme_classic()+
  facet_wrap(~pam50)  
  
saveRDS(res, "CellType_summary.rds")
