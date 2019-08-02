library(ggplot2)
library(cowplot)
library(topGO)
library(org.Mm.eg.db)

library(xlsx)
output <- data.frame()

for (i in (1:3)){
  Input <- read.xlsx("../data/FF99WT_combine_3clusters_MarkerGenes.xlsx",sheetIndex = i)
  
  deGenes <- Input[1:50,1]
  univrs <- Input[,1]
  
  alG <- factor(as.numeric(univrs %in% deGenes))
  names(alG) <- univrs
  
  Go.data <- new("topGOdata", description = "Go", ontology = "BP",
                 allGenes = alG, annot = annFUN.org, mapping = "org.Mm.eg.db",
                 nodeSize=20, ID = "symbol")
  result.classic <- runTest(Go.data, statistic = "fisher")
  tmp <- GenTable(Go.data, Fisher.classic =result.classic, 
                  orderBy="topgoFisher", topNodes =50, numChar = 300)
  tmp$Term <- factor(tmp$Term, levels=unique(rev(tmp$Term)))
  tmp$Clusters <- i
  
  output <- rbind(tmp, output)
  
}

output <- output[output$Fisher.classic < 0.01,]

#toPlot<- output %>% group_by(Clusters) %>% dplyr::top_n(-3, Fisher.classic)

ggplot(output, aes(x=Term, y=-log10(as.numeric(Fisher.classic)))) +
  geom_bar(stat="identity",color="black",fill="white") +
  coord_flip() +
  ylab("-log10(P)") +
  geom_hline(yintercept=2,lty="dashed") +
  xlab("GO-Term [BP]") +
  theme(axis.text.y=element_text(size=7)) +
  facet_wrap(Clusters~.,scales="free")

write.csv(output, "PyMT_3clusters_GO.csv")
