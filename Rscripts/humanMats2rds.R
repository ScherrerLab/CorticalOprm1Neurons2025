#human ssv4 count data to seurat objects
library(Seurat)
input <- 'human.matrix.csv'
metafile <- 'humna.metadata.csv'
s <- read.table(input,sep=',',row.names=1,header=T)
s <- t(s)
meta <- read.csv(metafile,header=T,row.names=1)
ssmeta <- meta[which(meta$region_label=='S1lm' | meta$region_label =='S1ul'),]
acmeta <- meta[which(meta$region_label == 'A1C'),]
s <- CreateSeuratObject(counts=s,meta.data = meta)

aca <- subset(s,cells=rownames(acmeta))
ss <- subset(s,cells=rownames(ssmeta))
Idents(aca) <- aca$subclass_label
Idents(ss) <- ss$subclass_label
clusters <- levels(aca)
debris <- c('','Microglia','Oligodendrocyte','Astrocyte','OPC','Endothelial','VLMC','Pericyte')
bad <- clusters[clusters %in% debris]
aca <- subset(aca,idents=bad,invert=T)

clusters <- levels(ss)
bad <- clusters[clusters %in% debris]
ss <- subset(ss, idents = bad, invert=T)

saveRDS(aca,'ACA_AIBSsubclass.rds')
saveRDS(ss, 'SS_AIBSsubclass.rds')
