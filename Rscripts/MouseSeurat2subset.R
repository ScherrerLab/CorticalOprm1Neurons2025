#mouse ssv4 full seurat object to ACA/SSp specific seurat objects
library(Seurat)
load('Seurat.ss.rda')
aca <- subset(ss.seurat,subset=region_label=='ACA')
ss <- subset(ss.seurat,subset=region_label=='SSp')

Idents(aca) <- aca$subclass_label
Idents(ss) <- ss$subclass_label
clusters <- levels(aca)
debris <- c('','Micro','Micro-PVM','Oligo','Astro','OPC','Endo','VLMC','SMC-Peri')
bad <- clusters[clusters %in% debris]
aca <- subset(aca,idents=bad,invert=T)

clusters <- levels(ss)
bad <- clusters[clusters %in% debris]
ss <- subset(ss, idents = bad, invert=T)

saveRDS(aca,'ACA_AIBSsubclass.rds')
saveRDS(ss, 'SS_AIBSsubclass.rds')
