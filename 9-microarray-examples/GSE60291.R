## paper results: Messenger RNA expression analysis identified 731 probe sets with significant differential expression
## mRNA expression array - GSE60291  (Affymetrix Human Genome U133 Plus 2.0 Array)

## Cell-line	GPL570	HG-U133_Plus_2
## different treatment :  ET1 control
rm(list=ls())
library(GEOquery)
library(limma) 
GSE60291 <- getGEO('GSE60291', destdir=".",getGPL = F)
exprSet=exprs(GSE60291[[1]])
library("annotate")
GSE60291[[1]]

pdata=pData(GSE60291[[1]])
treatment=factor(unlist(lapply(pdata$title,function(x) strsplit(as.character(x),"-")[[1]][1]))) 
#treatment=relevel(treatment,'control')

platformDB='hgu133plus2.db'
library(platformDB, character.only=TRUE)
probeset <- featureNames(GSE60291[[1]])
#EGID <- as.numeric(lookUp(probeset, platformDB, "ENTREZID"))
SYMBOL <-  lookUp(probeset, platformDB, "SYMBOL")

a=cbind(SYMBOL,exprSet)
## remove the duplicated probeset 
rmDupID <-function(a=matrix(c(1,1:5,2,2:6,2,3:7),ncol=6)){
  exprSet=a[,-1]
  rowMeans=apply(exprSet,1,function(x) mean(as.numeric(x),na.rm=T))
  a=a[order(rowMeans,decreasing=T),]
  exprSet=a[!duplicated(a[,1]),]
  #
  exprSet=exprSet[!is.na(exprSet[,1]),]
  rownames(exprSet)=exprSet[,1]
  exprSet=exprSet[,-1]
  return(exprSet)
}
exprSet=rmDupID(a)
rn=rownames(exprSet)
exprSet=apply(exprSet,2,as.numeric)
rownames(exprSet)=rn
exprSet[1:4,1:4]
#exprSet=log(exprSet) ## based on e 
boxplot(exprSet,las=2)

design=model.matrix(~ treatment)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
#vennDiagram(decideTests(fit))
DEG=topTable(fit,coef=2,n=Inf,adjust='BH')
dim(DEG[abs(DEG[,1])>1.2 & DEG[,5]<0.05,])  ## 806 genes
write.csv(DEG,"ET1-normal.DEG.csv")
