rm(list=ls())
library(reshape2)
library(edgeR)
library(DESeq2)

setwd("G:/mRNA/DEG")
a=read.table('hisat2_mm10_htseq.txt',stringsAsFactors = F)
######################################################################
#ESCTSA01.geneCounts	Nek1	2790
#ESCTSA01.geneCounts	Nek10	18
#ESCTSA01.geneCounts	Nek11	2
#ESCTSA01.geneCounts	Nek2	4945
######################################################################
colnames(a)=c('sample','gene','reads')
exprSet=dcast(a,gene~sample)
head(exprSet)

# write.table(exprSet[grep("^__",exprSet$gene),],'hisat2.stats.txt',quote=F,sep='\t')
# exprSet=exprSet[!grepl("^__",exprSet$gene),] 

geneLists=exprSet[,1]
exprSet=exprSet[,-1]
head(exprSet)

rownames(exprSet)=geneLists
colnames(exprSet)=do.call(rbind,strsplit(colnames(exprSet),'\\.'))[,1]

write.csv(exprSet,'raw_reads_matrix.csv') 

keepGene=rowSums(cpm(exprSet)>0) >=2
table(keepGene);dim(exprSet)
dim(exprSet[keepGene,])
exprSet=exprSet[keepGene,]
rownames(exprSet)=geneLists[keepGene]

str(exprSet)

group_list=c('control','control','treat_12','treat_12','treat_2','treat_2')

write.csv(exprSet,'filter_reads_matrix.csv' )
 


######################################################################
###################      Firstly for DEseq2      #####################
######################################################################
if(T){
  
  suppressMessages(library(DESeq2)) 
  (colData <- data.frame(row.names=colnames(exprSet), group_list=group_list) )
  dds <- DESeqDataSetFromMatrix(countData = exprSet,
                                colData = colData,
                                design = ~ group_list)
  dds <- DESeq(dds)
  png("qc_dispersions.png", 1000, 1000, pointsize=20)
  plotDispEsts(dds, main="Dispersion plot")
  dev.off()
  
  
  rld <- rlogTransformation(dds)
  exprMatrix_rlog=assay(rld) 
  write.csv(exprMatrix_rlog,'exprMatrix.rlog.csv' )
  
  normalizedCounts1 <- t( t(counts(dds)) / sizeFactors(dds) )
  # normalizedCounts2 <- counts(dds, normalized=T) # it's the same for the tpm value
  # we also can try cpm or rpkm from edgeR pacage
  exprMatrix_rpm=as.data.frame(normalizedCounts1) 
  head(exprMatrix_rpm)
  write.csv(exprMatrix_rpm,'exprMatrix.rpm.csv' )
  
  png("DEseq_RAWvsNORM.png",height = 800,width = 800)
  par(cex = 0.7)
  n.sample=ncol(exprSet)
  if(n.sample>40) par(cex = 0.5)
  cols <- rainbow(n.sample*1.2)
  par(mfrow=c(2,2))
  boxplot(exprSet, col = cols,main="expression value",las=2)
  boxplot(exprMatrix_rlog, col = cols,main="expression value",las=2)
  hist(as.matrix(exprSet))
  hist(exprMatrix_rlog)
  dev.off()
  
  library(RColorBrewer)
  (mycols <- brewer.pal(8, "Dark2")[1:length(unique(group_list))])
  cor(as.matrix(exprSet))
  # Sample distance heatmap
  sampleDists <- as.matrix(dist(t(exprMatrix_rlog)))
  #install.packages("gplots",repos = "http://cran.us.r-project.org")
  library(gplots)
  png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
  heatmap.2(as.matrix(sampleDists), key=F, trace="none",
            col=colorpanel(100, "black", "white"),
            ColSideColors=mycols[group_list], RowSideColors=mycols[group_list],
            margin=c(10, 10), main="Sample Distance Matrix")
  dev.off()
  
  cor(exprMatrix_rlog) 
  
  
  res <- results(dds, contrast=c("group_list","treat_2","control"))
  resOrdered <- res[order(res$padj),]
  head(resOrdered)
  DEG_treat_2=as.data.frame(resOrdered)
  write.csv(DEG_treat_2,"DEG_treat_2_deseq2.results.csv")
  
  res <- results(dds, contrast=c("group_list","treat_12","control"))
  resOrdered <- res[order(res$padj),]
  head(resOrdered)
  DEG_treat_12=as.data.frame(resOrdered)
  write.csv(DEG_treat_12,"DEG_treat_12_deseq2.results.csv")
  
  
  
}

######################################################################
###################      Then  for edgeR        #####################
######################################################################
if(T){
  
  library(edgeR)
  d <- DGEList(counts=exprSet,group=factor(group_list))
  d$samples$lib.size <- colSums(d$counts)
  d <- calcNormFactors(d)
  d$samples
  
  ## The calcNormFactors function normalizes for RNA composition by finding a set of scaling
  ## factors for the library sizes that minimize the log-fold changes between the samples for most
  ## genes. The default method for computing these scale factors uses a trimmed mean of Mvalues
  ## (TMM) between each pair of samples
  
  png('edgeR_MDS.png')
  plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
  dev.off()
  
  # The glm approach to multiple groups is similar to the classic approach, but permits more general comparisons to be made
  
  dge=d
 
  design <- model.matrix(~0+factor(group_list))
  rownames(design)<-colnames(dge)
  colnames(design)<-levels(factor(group_list))
  
  dge <- estimateGLMCommonDisp(dge,design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  
  fit <- glmFit(dge, design)
  
  lrt <- glmLRT(fit,  contrast=c(-1,1,0))
  nrDEG=topTags(lrt, n=nrow(exprSet))
  nrDEG=as.data.frame(nrDEG)
  head(nrDEG)
  write.csv(nrDEG,"DEG_treat_12_edgeR.csv",quote = F)
  
  lrt <- glmLRT(fit, contrast=c(-1,0,1) )
  nrDEG=topTags(lrt, n=nrow(exprSet))
  nrDEG=as.data.frame(nrDEG)
  head(nrDEG)
  write.csv(nrDEG,"DEG_treat_2_edgeR.csv",quote = F)
  summary(decideTests(lrt))
  plotMD(lrt)
  abline(h=c(-1, 1), col="blue")
}

######################################################################
###################      Then  for limma/voom        #################
######################################################################


suppressMessages(library(limma))
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprSet)

dge <- DGEList(counts=exprSet)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log=TRUE, prior.count=3)

v <- voom(dge,design,plot=TRUE, normalize="quantile")
fit <- lmFit(v, design)

group_list
cont.matrix=makeContrasts(contrasts=c('treat_12-control','treat_2-control'),levels = design)
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)
 
tempOutput = topTable(fit2, coef='treat_12-control', n=Inf)
DEG_treat_12_limma_voom = na.omit(tempOutput)
write.csv(DEG_treat_12_limma_voom,"DEG_treat_12_limma_voom.csv",quote = F)

tempOutput = topTable(fit2, coef='treat_2-control', n=Inf)
DEG_treat_2_limma_voom = na.omit(tempOutput)
write.csv(DEG_treat_2_limma_voom,"DEG_treat_2_limma_voom.csv",quote = F)



png("limma_voom_RAWvsNORM.png",height = 600,width = 600)
exprSet_new=v$E
par(cex = 0.7)
n.sample=ncol(exprSet)
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
par(mfrow=c(2,2))
boxplot(exprSet, col = cols,main="expression value",las=2)
boxplot(exprSet_new, col = cols,main="expression value",las=2)
hist(as.matrix(exprSet))
hist(exprSet_new)
dev.off()






