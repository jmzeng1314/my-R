## kidney	GPL96	HG-U133A
## different subtype tissue :Glomerulus Tubulus
## glomerulus apparently overrepresentation  compared with Tubulus
rm(list=ls())
library(GEOquery)
library(limma) 
GSE21785 <- getGEO('GSE21785', destdir=".",getGPL = F)
exprSet=exprs(GSE21785[[1]])
library("annotate")
GSE21785[[1]]

pdata=pData(GSE21785[[1]])
tissue=factor(unlist(lapply(pdata$title,function(x) strsplit(as.character(x),"_")[[1]][1]))) 

platformDB='hgu133a.db'
library(platformDB, character.only=TRUE)
probeset <- featureNames(GSE21785[[1]])
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
## very interesting here,one group is apparently bigger than other.

design=model.matrix(~ tissue)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
#vennDiagram(decideTests(fit))

write.csv(topTable(fit,coef='tissueTubulus',n=Inf,adjust='BH'),"Tubulus-Glomerulus.DEG.csv")
