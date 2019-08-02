# QC-Analysis

library(scran)
library(dplyr)
library(knitr)
library(ggplot2)
library(Rtsne)
library(cowplot)
source("../src/functions.R")

# Load Data
DataList <- readRDS("../data/Robjects/ExpressionList.rds")
m <- DataList[["counts"]]
pD <- DataList[["phenoData"]]
fD <- DataList[["featureData"]]
rm(DataList)

# ---- QCOverview ----

# Sequencing Depth and Genes detected
pD$UmiSums<- colSums(m)
pD$GenesDetected <- colSums(m!=0)
genesDetected <- ggplot(pD, aes(x=SampleID,y=GenesDetected)) +
    geom_violin(draw_quantiles=0.5)+
    scale_y_log10() +
    ylab("Total number of genes detected") +
    theme_bw()
genesDetected
LibrarySize <- ggplot(pD, aes(x=SampleID,y=UmiSums)) +
    geom_violin(draw_quantiles=0.5)+
    scale_y_log10() +
    ylab("Total number of molecules") +
    theme_bw()
LibrarySize
# Genewise QC
ntop <- 50 
mRel <- t(t(m)/colSums(m))
rownames(mRel)  <- fD$symbol
topExpressed <- rowMedians(mRel)
names(topExpressed) <- rownames(mRel)
topExpressed <- topExpressed %>% sort(.,decreasing=TRUE) %>% names
plotData <- t(mRel)[,topExpressed[1:ntop]] %>%
    reshape2::melt() 
colnames(plotData)<-c("Cell","Gene","RelativeExpression")

topGenes <- ggplot(plotData, aes(x=Gene, y=RelativeExpression)) +
    geom_boxplot() +
    coord_flip() +
    theme_bw()
freqOfExp <- m!=0
rownames(freqOfExp) <- fD$symbol
freqOfExp <- sort(rowSums(freqOfExp)/ncol(freqOfExp),decreasing=TRUE)
plotData <- data.frame("Gene"=names(freqOfExp),"Frequency"=freqOfExp)
topFreq <- ggplot(plotData[1:ntop,], aes(x=factor(Gene,levels=Gene),y=Frequency)) +
    geom_bar(stat="identity") +
    coord_flip() +
    xlab("Gene") +
    theme_bw()

# Cell Viability
mMito <- m[fD$Mitochondrial,]
idtop <- fD[fD$symbol %in% names(freqOfExp)[1:ntop],"id"]
mTop <- m[idtop,]!=0
pD$prcntTop <- colSums(mTop)/ntop
pD$prcntMito <- colSums(mMito)/colSums(m)
cellViability <- ggplot(pD, aes(x=prcntMito, y=GenesDetected,color=SampleID))+
    geom_point() 
cellViability
prcntTop <- ggplot(pD, aes(x=prcntTop, y=GenesDetected, color=SampleID))+
    geom_point() 
prcntTop
rm(mMito)
rm(mTop)

plot_grid(genesDetected,LibrarySize,cellViability,
	  topGenes)

# ---- Index-Swapping ----

barcodes <- as.character(pD$barcode)
spt <- strsplit(barcodes, split = "-", fixed = T)
pD$sample <- sapply(spt, function(x) x[2])
pD$bcs <- sapply(spt, function(x) x[1])
pD.add <- data.frame(bcs=names(table(pD$bcs)),
                     bc.obs=(table(pD$bcs)))
pD.add <- pD.add[,-2]
pD <- dplyr::left_join(pD,pD.add)

# P1
index.p1 <- ggplot(pD, aes(x=SampleID, fill=as.factor(bc.obs.Freq))) +
  geom_bar() +
  #     ggtitle("Samples contain many shared barcodes") +
  scale_fill_discrete(name="Times barcode observed") +
  theme(legend.position="bottom",
        legend.direction="horizontal")

compare <- function(barcodes, samples) {
  out <- NULL
  ids <- levels(samples)
  combs <- combn(ids,m=2, simplify=FALSE)
  for (i in seq_along(combs)) {
    comb <- combs[[i]]
    s1 <- comb[1]
    s2 <- comb[2]
    bc1 <- as.character(barcodes[samples==s1])
    bc2 <- as.character(barcodes[samples==s2])
    x <- length(intersect(bc1,bc2))
    m <- length(bc1)
    n <- 750000
    k <- length(bc2)
    p.val <- phyper(q=x-1,m=m,n=n,k=k,lower.tail=FALSE)
    tmp <- data.frame(s1=s1,
                      s2=s2,
                      n1=m,
                      n2=k,
                      shared=x,
                      p.val=p.val)
    out <- rbind(out,tmp)
  }
  return(out)
}

# P2
pD$SampleID <- factor(pD$SampleID, levels = c("BRCA1_A",  "BRCA1_B",  "NEU_A",    "NEU_B",    "FF99WT_A", "FF99WT_B", "4T1" ))
compDf <- compare(pD$bcs, pD$SampleID)
index.p2 <- ggplot(compDf, aes(x=shared, y=-log10(p.val))) +
  geom_point() +
  xlab("# shared Barcodes") +
  ylab("-log10(P)") +
  geom_hline(yintercept=2,lty="dashed",color="red") 

# ---- QCThresholding ----

#Fixed threshold
MitoCutOff <- 0.2

#Thresholding per Condition
pD$Condition<-as.factor(pD$SampleID)
pD <- mutate(pD,
	     ThresholdViability = 0,
	     ThresholdGenesDet = 0,
	     ThresholdLibSize = 0)


smryByGroup <- pD %>% dplyr::group_by(Condition) %>%
  dplyr::summarise(threshold_PrcntMito=MitoCutOff,
             mGenesDetected=median(log10(GenesDetected)),
             madGenesDetected=mad(log10(GenesDetected)),
             threshold_GenesDetected=max(mGenesDetected-3*madGenesDetected,log10(500)),
             mUmiSums=median(log10(UmiSums)),
             madUmiSums=mad(log10(UmiSums)),
             threshold_UmiSums=max(mUmiSums-3*madUmiSums,log10(1000))) %>%
            dplyr::select(Condition,starts_with("threshold")) 
kable(smryByGroup)

grps <- as.character(unique(pD$Condition))
grps
#"BRCA1_A"  "BRCA1_B"  "NEU_A"    "NEU_B"    "FF99WT_A" "FF99WT_B" "4T1"    

for (grp in grps) {
  thrs <- filter(smryByGroup, Condition == grp) %>%
    select(-Condition) %>% t() %>%
    as.vector()
  names(thrs) <- filter(smryByGroup, Condition == grp) %>%
    select(-Condition) %>% t() %>%
    rownames()
  pD <- mutate(pD,
               ThresholdViability= ifelse(Condition==grp, MitoCutOff, ThresholdViability),
               ThresholdGenesDet= ifelse(Condition==grp, 10^thrs["threshold_GenesDetected"],ThresholdGenesDet),
               ThresholdLibSize= ifelse(Condition==grp, 10^thrs["threshold_UmiSums"],ThresholdLibSize))
}

pD <- mutate(pD, 
             PassViability=prcntMito < ThresholdViability,
             PassGenesDet=GenesDetected > ThresholdGenesDet,
             PassLibSize=UmiSums > ThresholdLibSize,
             PassAll= PassViability & PassGenesDet & PassLibSize & bc.obs.Freq==1)

# Illustrate thresholds
gdHist <- ggplot(pD, aes(x=GenesDetected,y=..density..)) +
    geom_histogram(fill="white",color="black",bins=100) +
    geom_vline(data=smryByGroup,aes(xintercept=200),color="red",lty="longdash") +
    scale_x_log10() +
    xlab("Total number of genes detected") +
    facet_wrap(~Condition) 

libSizeHist <- ggplot(pD, aes(x=UmiSums,y=..density..)) +
    geom_histogram(fill="white",color="black",bins=100) +
    geom_vline(data=smryByGroup,aes(xintercept=10^threshold_UmiSums),color="red",lty="longdash") +
    geom_vline(data=smryByGroup,aes(xintercept=10^2.825158),color="red",lty="longdash") + 
    scale_x_log10() +
    facet_wrap(~Condition) +
    xlab("Total number of unique molecules") 
mitoprHist <- ggplot(pD, aes(x=prcntMito,y=..density..)) +
  geom_histogram(fill="white",color="black",bins=100) +
  geom_vline(data=smryByGroup,aes(xintercept=10^threshold_UmiSums),color="red",lty="longdash") +
  geom_vline(data=smryByGroup,aes(xintercept=0.2),color="red",lty="longdash") + 
  scale_x_log10() +
  facet_wrap(~Condition) +
  xlab("Total number of unique molecules")

cellViability <- cellViability %+% pD
cellViability <- cellViability + 
    annotate("rect",ymin=-Inf, ymax=Inf, xmax=Inf, xmin=0.2,
	     fill="grey", alpha=0.3) 


# Overview over cells removed
table(pD$Condition,pD$PassGenesDet)
table(pD$Condition,pD$PassLibSize)
table(pD$Condition,pD$PassViability)
table(pD$Condition,pD$PassAll)

plot_grid(gdHist,libSizeHist,cellViability)

# Summary table---------
sumry <- dplyr::group_by(pD, SampleID) %>%
  dplyr::summarise("Number of cells"=n(),
            "Total molecules"=median(UmiSums),
            "Genes Detected"=median(GenesDetected),
            "After QC" = sum(PassAll==TRUE))
tmp =NULL
stats=NULL
for (grp in grps){
  tmp <- read.table(paste0("../../",grp,"/outs/metrics_summary.csv"), sep=",", header=TRUE, check.names = FALSE)
  stat <- data.frame(SampleID=grp, NumberOfReads = tmp$`Number of Reads`, Saturation=tmp$`Sequencing Saturation`)
  stats <- rbind(stats, stat)
}

sumry <- left_join(sumry, stats[,c("SampleID","NumberOfReads","Saturation")])

library(gridExtra)
p1 <- tableGrob(sumry,rows=NULL)

# fD
#fD$keep <- rowMeans(m) > 0.01

# Save Data
stopifnot(identical(as.character(pD$barcode),colnames(m)))
stopifnot(identical(as.character(fD$id),rownames(m)))
out <- list()
out[["counts"]] <- m
out[["phenoData"]] <- pD
out[["featureData"]] <- fD
saveRDS(out,file="../data/Robjects/ExpressionList_QC.rds")
