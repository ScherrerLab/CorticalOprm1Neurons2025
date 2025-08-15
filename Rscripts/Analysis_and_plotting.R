#compare mouse human and generate plots
library(biomaRt)
library(dplyr)
library(Seurat)
library(ggplot2)
library(gtools)
library(gplots)
library(matrixStats)
library(dendextend)
library(scCustomize)
library(dittoSeq)
library(patchwork)
#preprocess data for comparable UMAPs and isolate z-scored log1p count data (scale)
huaca=readRDS('Human/ACA_AIBSsubclass.rds')
huss <- readRDS('Human/SS_AIBSsubclass.rds')
huaca=FindVariableFeatures(huaca,nfeatures = 3000) %>% NormalizeData() %>% ScaleData() %>% RunPCA(npcs=10) %>% RunUMAP(dims=1:10)
huss=FindVariableFeatures(huss,nfeatures = 3000) %>% NormalizeData() %>% ScaleData() %>% RunPCA(npcs=10) %>% RunUMAP(dims=1:10)
huacascaled=as.matrix(huaca[['RNA']]$scale.data)
hussscale <- as.matrix(huss[['RNA']]$scale.data)

#ms=readRDS('ACA/ACA_clustered_all.rds')
msaca <- readRDS('ACA_ms_AIBSsubclass.rds')
msss <- readRDS('SS_ms_AIBSsubclass.rds')
meta <- read.csv('/proj/gs25/projects/Allen_SmartSeq/U19.SS.metadat.allcol.csv',header=T,row.names=1)
metaaca <- meta[rownames(meta) %in% colnames(msaca),]
metass <- meta[rownames(meta) %in% colnames(msss),]
msaca$full_genotype <- metaaca$full_genotype
msss$full_genotype <- metass$full_genotype
msaca <- FindVariableFeatures(msaca,nfeatures = 3000) %>% NormalizeData() %>% ScaleData() %>% RunPCA(npcs=10) %>% RunUMAP(dims=1:10)
msss <- FindVariableFeatures(msss,nfeatures = 3000) %>% NormalizeData() %>% ScaleData() %>% RunPCA(npcs=10) %>% RunUMAP(dims=1:10)
#plot umaps
p1 <- DimPlot_scCustom(msaca,pt.size = .1,label=F) + NoAxes() + ggtitle(label = 'Mouse ACA') + theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot_scCustom(msss,pt.size = .1,label=F) + NoAxes() + ggtitle(label = 'Mouse SSp') + theme(plot.title = element_text(hjust = 0.5))
h1 <- DimPlot_scCustom(huaca,pt.size = .1,label=F) + NoAxes() + ggtitle(label = 'Human ACA') + theme(plot.title = element_text(hjust = 0.5))
h2 <- DimPlot_scCustom(msss,pt.size = .1,label=F) + NoAxes() + ggtitle(label = 'Human SSp') + theme(plot.title = element_text(hjust = 0.5))
pdf('fourUMAPs.pdf',width=14,height=7)
print(p1 + p2 + h1 + h2 + plot_layout(ncol=2))
dev.off()
#of all Excitatory (e.g. Slc17a7+, Slc17a6, anud Snap25 + Slc17a7+ neurons that express Oprm1, what percent are IT/PT/CT/
#msacaExcit <- subset

#Get centroids
hgacaavg=AverageExpression(huaca,slot = 'scale.data')
hgacaavg=hgacaavg$RNA
hgssavg=AverageExpression(huss,slot = 'scale.data')
hgssavg=hgssavg$RNA
msacaavg=AverageExpression(msaca,slot = 'scale.data')
msacaavg=msacaavg$RNA
msssavg=AverageExpression(msss,slot = 'scale.data')
msssavg=msssavg$RNA
#gene ms/hu gene keys
mouse_human_genes_all=read.table('http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt',header=T,sep='\t')
mouse_human_genes=mouse_human_genes_all[,c(1,2,4)]
mouse_genes=mouse_human_genes[which(mouse_human_genes$Common.Organism.Name=='mouse, laboratory'),]
human_genes=mouse_human_genes[which(mouse_human_genes$Common.Organism.Name=='human'),]
#get list of DB.Class.Keys present in both; subset by overlap
mouse_dedup=mouse_genes[-which(duplicated(mouse_genes$Symbol)),]
human_dedup=human_genes[-which(duplicated(human_genes$Symbol)),]


mskeys=mouse_dedup$DB.Class.Key
hgkeys=human_dedup$DB.Class.Key
interkey=intersect(mskeys,hgkeys)

ms.conserved=mouse_dedup[mskeys %in% interkey,]
hg.conserved=human_dedup[hgkeys %in% interkey,]
if (length(which(duplicated(ms.conserved$DB.Class.Key)))>0) {
  ms.fin=ms.conserved[-which(duplicated(ms.conserved$DB.Class.Key)),]
} else {
  ms.fin=ms.conserved
}
if (length(which(duplicated(hg.conserved$DB.Class.Key)))>0) {
  hg.fin=hg.conserved[-which(duplicated(hg.conserved$DB.Class.Key)),]
} else {
  hg.fin=hg.conserved
}
#subset avg datasets to genes that have orthologs
hgacasub=as.matrix(hgacaavg[rownames(hgacaavg) %in%hg.fin$Symbol,])
hgsssub=as.matrix(hgssavg[rownames(hgssavg) %in%hg.fin$Symbol,])
msacasub=as.matrix(msacaavg[rownames(msacaavg) %in%ms.fin$Symbol,])
mssssub=as.matrix(msssavg[rownames(msssavg) %in%ms.fin$Symbol,])
ms.aca.fin.sub <- ms.fin[ms.fin$Symbol %in% rownames(msacasub),]
ms.ss.fin.sub <- ms.fin[ms.fin$Symbol %in% rownames(mssssub),]
hg.aca.fin.sub <- hg.fin[hg.fin$Symbol %in% rownames(hgacasub),]
hg.ss.fin.sub <- hg.fin[hg.fin$Symbol %in% rownames(hgsssub),]
#make ms genes into human genes based on DB#s

ma=msacasub[match(ms.aca.fin.sub$Symbol,rownames(msacasub)),]#reorder genes to match key
ma=ma[complete.cases(ma),]
ms=mssssub[match(ms.ss.fin.sub$Symbol,rownames(mssssub)),]#reorder genes to match key
ms=ms[complete.cases(ms),]
ha=hgacasub[match(hg.aca.fin.sub$Symbol,rownames(hgacasub)),]
ha=ha[complete.cases(ha),]
hs=hgsssub[match(hg.ss.fin.sub$Symbol,rownames(hgsssub)),]
hs=hs[complete.cases(hs),]
#add corresponding human to mouse and vice versa
ms.aca.fin.sub$human=hg.fin$Symbol[hg.fin$DB.Class.Key %in% ms.aca.fin.sub$DB.Class.Key]
rownames(ms.aca.fin.sub)=ms.aca.fin.sub$Symbol
hg.aca.fin.sub$mouse=ms.fin$Symbol[ms.fin$DB.Class.Key %in% hg.aca.fin.sub$DB.Class.Key]
rownames(hg.aca.fin.sub)=hg.aca.fin.sub$Symbol
ms.ss.fin.sub$human=hg.fin$Symbol[hg.fin$DB.Class.Key %in% ms.ss.fin.sub$DB.Class.Key]
rownames(ms.ss.fin.sub)=ms.ss.fin.sub$Symbol
hg.ss.fin.sub$mouse=ms.fin$Symbol[ms.fin$DB.Class.Key %in% hg.ss.fin.sub$DB.Class.Key]
rownames(hg.ss.fin.sub)=hg.ss.fin.sub$Symbol
#subset human key by human df
#subset mouse key by human genes and vice verssa
ms.aca.fin.sub=ms.aca.fin.sub[rownames(ms.aca.fin.sub) %in% hg.aca.fin.sub$mouse,]
hg.aca.fin.sub=hg.aca.fin.sub[rownames(hg.aca.fin.sub) %in% ms.aca.fin.sub$human,]
ms.ss.fin.sub=ms.ss.fin.sub[rownames(ms.ss.fin.sub) %in% hg.ss.fin.sub$mouse,]
hg.ss.fin.sub=hg.ss.fin.sub[rownames(hg.ss.fin.sub) %in% ms.ss.fin.sub$human,]
#subset mouse df by human genes
ma.fin=ma[rownames(ma) %in% rownames(ms.aca.fin.sub),]
ms.fin=ms[rownames(ms) %in% rownames(ms.ss.fin.sub),]
ha.fin <- ha[rownames(ha) %in% rownames(hg.aca.fin.sub),]
hs.fin <- hs[rownames(hs) %in% rownames(hg.ss.fin.sub),]
#make gene names equivalent; find most variable genes b/w each;
library(matrixStats)
#subset human df by mouse genes
aa=ha[rownames(ha) %in% rownames(hg.aca.fin.sub),]
as=hs[rownames(hs) %in% rownames(hg.ss.fin.sub),]
ba <- ma[rownames(ma) %in% rownames(ms.aca.fin.sub),]
bs <- ms[rownames(ms) %in% rownames(ms.ss.fin.sub),]
#renaming based on layer spec markers and AIBS seq atlas
#rename seurat objects

hc1<-hclust(as.dist(1-cor(aa)))
hc2<-hclust(as.dist(1-cor(ba)))
d1=as.dendrogram(hc1)
d2=as.dendrogram(hc2)
plot(d1)
plot(d2)

cormat1=cor(aa,ba,method = 'pearson') #l6bmouse and CTGF human match; larger correlation range

hc3<-hclust(as.dist(1-cor(as)))
hc4<-hclust(as.dist(1-cor(bs)))
d3=as.dendrogram(hc3)
d4=as.dendrogram(hc4)
plot(d3)
plot(d4)

cormat2=cor(as,bs,method = 'pearson') #l6bmouse and CTGF human match; larger correlation range

msoa <- Percent_Expressing(msaca,features='Oprm1',assay='RNA')
msos <- Percent_Expressing(msss,features='Oprm1',assay='RNA')

hsoa <- Percent_Expressing(huaca,features='OPRM1',assay='RNA')
hsos <- Percent_Expressing(huss,features='OPRM1',assay='RNA')

hsoa_bold <- hsoa >= 50  # Logical vector for rows
msoa_bold <- msoa >= 50  # Logical vector for columns
hsos_bold <- hsos >= 50  # Logical vector for rows
msos_bold <- msos >= 50  # Logical vector for columns
# Define row and column highlights
row_fonta <- ifelse(hsoa_bold, 2, 1)
col_fonta <- ifelse(msoa_bold, 2, 1)
row_fonts <- ifelse(hsos_bold, 2, 1)
col_fonts <- ifelse(msos_bold, 2, 1)

#condensed gradient based on quantiles
quan <- c(min(cormat1),0,max(cormat1)/2,max(cormat1))
colors=unique(c(seq(quan[1],quan[2],length=50),seq(quan[2],quan[3],length=50),seq(quan[3],quan[4],length=50)))
#linear gradient 
#m=diff(range(cormat1))/3
#quan=c(min(cormat1),min(cormat1)+m,min(cormat1)+m+m,max(cormat1))
#colors=unique(c(seq(quan[1],quan[2],length=50),seq(quan[2],quan[3],length=50),seq(quan[3],quan[4],length=50)))
my_palette<-colorRampPalette(c('darkblue','white','firebrick','darkred'))(n=147)
#par(mar=c(7,4,4,2)+0.1) 
pdf('ACA_MsHuHeatmap.pdf',width=5.8,height=4.5)
print(
heatmap.2(as.matrix(cormat1),
          scale='none',trace='none',density.info = 'none',
          margins = c(12,14),
          cexRow = .75,cexCol = .75,
          col=my_palette,breaks=colors,
          key.par = list(cex=0.5),
          main=paste0('HuMs ACA Correlation: \n',nrow(aa),' variable genes'),
          Rowv=1:12,
          labRow=rownames(cormat1),
          labCol=colnames(cormat1),
          colRow=row_fonta,
          colCol = col_fonta,)
)
dev.off()


#condensed gradient based on quantiles
quan <- c(min(cormat2),0,max(cormat2)/2,max(cormat2))
colors=unique(c(seq(quan[1],quan[2],length=50),seq(quan[2],quan[3],length=50),seq(quan[3],quan[4],length=50)))
#linear gradient 
#m=diff(range(cormat1))/3
#quan=c(min(cormat1),min(cormat1)+m,min(cormat1)+m+m,max(cormat1))
#colors=unique(c(seq(quan[1],quan[2],length=50),seq(quan[2],quan[3],length=50),seq(quan[3],quan[4],length=50)))
my_palette<-colorRampPalette(c('darkblue','white','firebrick','darkred'))(n=147)
#par(mar=c(7,4,4,2)+0.1) 
pdf('SSp_MsHuHeatmap.pdf',width=5.8,height=4.5)
print( 
heatmap.2(as.matrix(cormat2),
          scale='none',trace='none',density.info = 'none',
          margins = c(12,14),
          cexRow = .75,cexCol = .75,
          col=my_palette,breaks=colors,
          key.par = list(cex=0.5),
          main=paste0('HuMs SSp Correlation: \n',nrow(aa),' variable genes'),
          Rowv=1:12,
          labRow=rownames(cormat2),
          labCol=colnames(cormat2),
          colRow=row_fonts,
          colCol = col_fonts)
)
dev.off()


#MS: L5-2, L4-14, & Sst 12
#HU: L24-2, L5-4, SST-7
pdf('MsACAdotplot.pdf',width=4.5,height=3)
print(DotPlot(msaca,features=c('Oprm1','Sst','Cplx3','Ctgf','Foxp2','Sulf1'),cols = c('blue','red')) + theme(axis.text.x=element_text(angle=45,hjust=1)))
dev.off()
pdf('MsSSpdotplot.pdf',width=4.5,height=3)
print(DotPlot(msss,features=c('Oprm1','Sst','Cplx3','Ctgf','Foxp2','Sulf1'),cols = c('blue','red')) + theme(axis.text.x=element_text(angle=45,hjust=1)))
dev.off()
pdf('HuACAdotplot.pdf',width=4.5,height=3)
print(DotPlot(huaca,features=toupper(c('Oprm1','Sst','Cplx3','Ctgf','Foxp2','Sulf1')),cols = c('blue','red')) + theme(axis.text.x=element_text(angle=45,hjust=1)))
dev.off()
pdf('HuSSpdotplot.pdf',width=4.5,height=3)
print(DotPlot(huss,features=toupper(c('Oprm1','Sst','Cplx3','Ctgf','Foxp2','Sulf1')),cols = c('blue','red')) + theme(axis.text.x=element_text(angle=45,hjust=1)))
dev.off()

#Oprm1+ cell enrichment; excit vs inhib?
#L5 and L6 excitatory, not 6b
msacae <- subset(msaca,idents=levels(msaca)[c(1,3,5,6,7,8)])
mssse <- subset(msss,idents=levels(msss)[c(2,3,6,7,10,11)])
huacae <- subset(huaca,idents=levels(huaca)[c(3,4,9,12)])
husse <- subset(huss,idents=levels(huss)[c(7,8,9,10)])
msacae$Oprm1Expressing <- ifelse(msacae[['RNA']]$counts['Oprm1',]>0,'Oprm1+','Oprm1-')
mssse$Oprm1Expressing <- ifelse(mssse[['RNA']]$counts['Oprm1',]>0,'Oprm1+','Oprm1-')
huacae$Oprm1Expressing <- ifelse(huacae[['RNA']]$counts['OPRM1',]>0,'OPRM1+','OPRM1-')
husse$Oprm1Expressing <- ifelse(husse[['RNA']]$counts['OPRM1',]>0,'OPRM1+','OPRM1-')

Idents(msacae) <- msacae$Oprm1Expressing
Idents(mssse) <- mssse$Oprm1Expressing
msacae@misc$OpMark <- FindMarkers(msacae,ident.1 = 'Oprm1+',only.pos=TRUE,logfc.threshold=1)
mssse@misc$OpMark <- FindMarkers(mssse,ident.1 = 'Oprm1+',only.pos=TRUE,logfc.threshold=1)
Idents(msacae) <- msacae$subclass_label
Idents(mssse) <- mssse$subclass_label

DotPlot(msacae,features=unique(head(rownames(msacae@misc$OpMark),20)),cols=c('blue','red')) + theme(axis.text.x=element_text(angle=45,hjust=1))
DotPlot(mssse,features=unique(head(rownames(mssse@misc$OpMark),20)),cols=c('blue','red')) + theme(axis.text.x=element_text(angle=45,hjust=1))

#Of all Excitatory, Oprm1+ cells, what proportion are IT/PT/CT/NP
table(msacae$full_genotype)
table(mssse$full_genotype)
msacae$Projection <- substr(Idents(msacae),nchar(as.character(Idents(msacae)))-5,nchar(as.character(Idents(msacae)))-4)
mssse$Projection <- substr(Idents(mssse),nchar(as.character(Idents(mssse)))-5,nchar(as.character(Idents(mssse)))-4)
abar <- dittoBarPlot(msacae,var='Oprm1Expressing',group.by='Projection',color.panel = c('grey','lightblue'),data.out = T)
sbar <- dittoBarPlot(mssse,var='Oprm1Expressing',group.by='Projection',color.panel = c('grey','lightblue'),data.out = T)
pdf('MsACA_stackedByPercOprm1.pdf',width=2.6,height=2.25)
print(abar$p + geom_text(aes(label=round(percent*100,digits=0)),size=3,position=position_stack(vjust=0.5)))
dev.off()
pdf('MsSSp_stackedByPercOprm1.pdf',width=2.6,height=2.25)
print(sbar$p + geom_text(aes(label=round(percent*100,digits=0)),size=3,position=position_stack(vjust=0.5)))
dev.off()
pdf('MsACA_stackedByProj.pdf',width=2,height=3.5)
habar <- dittoBarPlot(msacae,var='Projection',group.by='Oprm1Expressing',data.out = T,color.panel = c('lavender','lightgreen','lightblue','pink'))
print(habar$p + geom_text(aes(label=round(percent*100,digits = 0)),size=3, position = position_stack(vjust = 0.5)))
dev.off()
pdf('MsSSp_stackedByProj.pdf',width=2,height=3.5)
hsbar <- dittoBarPlot(mssse,var='Projection',group.by='Oprm1Expressing',data.out = T,color.panel = c('lavender','lightgreen','lightblue','pink'))
print(hsbar$p + geom_text(aes(label=round(percent*100,digits = 0)),size=3, position = position_stack(vjust = 0.5)))
dev.off()


dittoBarPlot(msacae,var='Oprm1Expressing',group.by='Projection')
dittoBarPlot(msacae,var='Projection',group.by='Oprm1Expressing')

dittoBarPlot(mssse,var='Oprm1Expressing',group.by='Projection')
dittoBarPlot(mssse,var='Projection',group.by='Oprm1Expressing')

huacae$Projection=0
huacae$Projection[which(huacae$subclass_label=='L6 CT')] <- 'CT'
huacae$Projection[which(huacae$subclass_label=='L5/6 IT Car3')] <- 'IT'
huacae$Projection[which(huacae$subclass_label=='L5/6 NP')] <- 'NP'
huacae$Projection[which(huacae$subclass_label == 'L5 ET')] <- 'PT'
husse$Projection <- 0
husse$Projection[which(husse$subclass_label=='L6 CT')] <- 'CT'
husse$Projection[which(husse$subclass_label=='L5/6 IT Car3')] <- 'IT'
husse$Projection[which(husse$subclass_label=='L5/6 NP')] <- 'NP'
husse$Projection[which(husse$subclass_label == 'L5 ET')] <- 'PT'

abar <- dittoBarPlot(huacae,var='Oprm1Expressing',group.by='Projection',color.panel = c('grey','lightblue'),data.out = T)
sbar <- dittoBarPlot(husse,var='Oprm1Expressing',group.by='Projection',color.panel = c('grey','lightblue'),data.out = T)

pdf('HuACA_stackedByPercOprm1.pdf',width=2.6,height=2.25)
print(abar$p + geom_text(aes(label=round(percent*100,digits=0)),size=3,position=position_stack(vjust=0.5)))
dev.off()
pdf('HuSSp_stackedByPercOprm1.pdf',width=2.6,height=2.25)
print(sbar$p + geom_text(aes(label=round(percent*100,digits=0)),size=3,position=position_stack(vjust=0.5)))
dev.off()

habar <- dittoBarPlot(msacae,var='Projection',group.by='Oprm1Expressing',data.out = T,color.panel = c('lavender','lightgreen','lightblue','pink'))
hsbar <- dittoBarPlot(mssse,var='Projection',group.by='Oprm1Expressing',data.out = T,color.panel = c('lavender','lightgreen','lightblue','pink'))
pdf('HuACA_stackedByProj.pdf',width=2,height=3.5)
print(habar$p + geom_text(aes(label=round(percent*100,digits = 0)),size=3, position = position_stack(vjust = 0.5)))
dev.off()
pdf('HuSSp_stackedByProj.pdf',width=2,height=3.5)
print(hsbar$p + geom_text(aes(label=round(percent*100,digits = 0)),size=3, position = position_stack(vjust = 0.5)))
dev.off()


dittoBarPlot(huacae,var='Oprm1Expressing',group.by='Projection')
dittoBarPlot(huacae,var='Projection',group.by='Oprm1Expressing')

dittoBarPlot(husse,var='Oprm1Expressing',group.by='Projection')
dittoBarPlot(husse,var='Projection',group.by='Oprm1Expressing')

pl1 <- c()
pl2 <- c()
Idents(msacae) <- msacae$Projection
Idents(mssse) <- mssse$Projection
Idents(huacae) <- huacae$Projection
Idents(husse) <- husse$Projection
proj <- c('IT','CT','PT','NP')
for (i in 1:4){
  pl1[[i]] <- euler_plot(msacae,metacol='Projection',metavar=proj[i],genes='Oprm1',title=proj[i])
  pl2[[i]] <- euler_plot(mssse,metacol='Projection',metavar=proj[i],genes='Oprm1',title=proj[i])
  pl1[[4+i]] <- euler_plot(huacae,metacol='Projection',metavar=proj[i],genes='OPRM1',title=proj[i])
  pl2[[4+i]] <- euler_plot(husse,metacol='Projection',metavar=proj[i],genes='OPRM1',title=proj[i])
}
library(gridExtra)
pdf('projEulers_ACA.pdf',width=8,height=5)
print(grid.arrange(grobs = pl1, nrow = 2, ncol = 4))
dev.off()
pdf('projEulers_SSp.pdf',width=8,height=5)
print(grid.arrange(grobs = pl2, nrow = 2, ncol = 4))
dev.off()

#get percentages of L5/6 IT/PT/CT cells expressing Oprm1
#IT
#First subset by deep laminar neurons classified by IT, PT, or CT
#for mouse, subset by CTX (removes inhibitory) then remove anything L6b or L2/3
good <- levels(msaca)[-grep('L6b',levels(msaca))]
good <- good[grepl('L4/5',good) | grepl('L5',good) | grepl('L6',good)]
msaca <- subset(msaca, idents = good)

good <- levels(msss)[-grep('L6b',levels(msss))]
good <- good[grepl('L4/5',good) | grepl('L5',good) | grepl('L6',good)]
msss <- subset(msss, idents = good)

#for human, remove L6b. then subset by L5 or L6-containing
good <- levels(huaca)[-grep('L6b',levels(huaca))]
good <- good[grepl('L5',good) | grepl('L6',good)]
huaca <- subset(huaca,idents=good)
good <- levels(huss)[-grep('L6b',levels(huss))]
good <- good[grepl('L5',good) | grepl('L6',good)]
huss <- subset(huss,idents=good)
#then subset by expected projection target
it <- levels(msaca)[grepl('IT',levels(msaca))]
msacait <- subset(msaca,idents=it)
it <- levels(msss)[grepl('IT',levels(msss))]
msssit <- subset(msss,idents=it)
it <- levels(huaca)[grepl('IT',levels(huaca))]
huacait <- subset(huaca,idents=it)
it <- levels(huss)[grepl('IT',levels(huss))]
hussit <- subset(huss,idents=it)
#CT
CT <- levels(msaca)[grepl('CT ',levels(msaca))]
msacact <- subset(msaca,idents=CT)
CT <- levels(msss)[grepl('CT ',levels(msss))]
msssct <- subset(msss,idents=CT)
CT <- levels(huaca)[grepl('CT',levels(huaca))]
huacact <- subset(huaca,idents=CT)
CT <- levels(huss)[grepl('CT',levels(huss))]
hussct <- subset(huss,idents=CT)
#ET/PT
pt <- levels(msaca)[grepl('PT',levels(msaca)) | grepl('ET',levels(msaca))]
msacapt <- subset(msaca,idents=pt)
pt <- levels(msss)[grepl('PT',levels(msss)) | grepl('ET',levels(msss))]
mssspt <- subset(msss,idents=pt)
pt <- levels(huaca)[grepl('PT',levels(huaca)) | grepl('ET',levels(huaca))]
huacapt <- subset(huaca,idents=pt)
pt <- levels(huss)[grepl('PT',levels(huss)) | grepl('ET',levels(huss))]
husspt <- subset(huss,idents=pt)
#NP
NP <- levels(msaca)[grepl('NP ',levels(msaca))]
msacanp <- subset(msaca,idents=NP)
NP <- levels(msss)[grepl('NP ',levels(msss))]
msssnp <- subset(msss,idents=NP)
NP <- levels(huaca)[grepl('NP',levels(huaca))]
huacanp <- subset(huaca,idents=NP)
NP <- levels(huss)[grepl('NP',levels(huss))]
hussnp <- subset(huss,idents=NP)

pl1 <- c()
pl2 <- c()
mslist <- c(msacait,msacact,msacapt,msacanp)
hulist <- c(huacait,huacact,huacapt,huacanp)
mssslist <- c(msssit,msssct,mssspt,msssnp)
husslist <- c(hussit,hussct,husspt,hussnp)
proj <- c('IT','CT','PT/ET','NP')
#out of each excitatory projection-specific cell type, what percentage is Oprm1+?
for (i in 1:4){
  pl1[[i]] <- euler_plot(mslist[[i]],'Oprm1',title = proj[i])
  pl1[[i+4]] <- euler_plot(hulist[[i]],'OPRM1')
  pl2[[i]] <- euler_plot(mssslist[[i]],'Oprm1',title=proj[i])
  pl2[[i+4]] <- euler_plot(husslist[[i]],'OPRM1')
}

library(gridExtra)
pdf('projEulers_ACA.pdf',width=8,height=5)
print(grid.arrange(grobs = pl1, nrow = 2, ncol = 4))
dev.off()
pdf('projEulers_SSp.pdf',width=8,height=5)
print(grid.arrange(grobs = pl2, nrow = 2, ncol = 4))
dev.off()
