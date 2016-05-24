## CLL cells	GPL570	HG-U133_Plus_2
## different treatment : non-del11q del11q
rm(list=ls())
library(GEOquery)
library(limma) 
GSE26526 <- getGEO('GSE26526', destdir=".",getGPL = F)
exprSet=exprs(GSE26526[[1]])
library("annotate")
GSE26526[[1]]

pdata=pData(GSE26526[[1]])
treatment=factor(unlist(lapply(pdata$title,function(x) strsplit(as.character(x),"[_ ]")[[1]][2]))) 
treatment=relevel(treatment,'non-del11q')
platformDB='hgu133plus2.db'
library(platformDB, character.only=TRUE)
probeset <- featureNames(GSE26526[[1]])
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

write.csv(topTable(fit,coef=2,n=Inf,adjust='BH'),"del11q-normal.DEG.csv")
 
