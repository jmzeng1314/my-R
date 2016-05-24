## obesity adipose	GPL570	HG-U133_Plus_2
## different subtype tissue :omental subcutaneous
## different reactive for insulin : resistant sensitive
rm(list=ls())
library(GEOquery)
library(limma) 
GSE20950 <- getGEO('GSE20950', destdir=".",getGPL = F)
exprSet=exprs(GSE20950[[1]])
library("annotate")
GSE20950[[1]]

pdata=pData(GSE20950[[1]])
tissue=factor(unlist(lapply(pdata$title,function(x) strsplit(as.character(x)," ")[[1]][6])))
insulin=factor(unlist(lapply(pdata$title,function(x) strsplit(as.character(x)," ")[[1]][2])))

platformDB='hgu133plus2.db'
library(platformDB, character.only=TRUE)
probeset <- featureNames(GSE20950[[1]])
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

design=model.matrix(~ tissue+insulin)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
vennDiagram(fit,include=c("up","down"))

write.csv(topTable(fit,coef='tissuesubcutaneous',n=Inf,adjust='BH'),"subcutaneous-omental.DEG.csv")
write.csv(topTable(fit,coef='insulinsensitive',n=Inf,adjust='BH')  ,"sensitive-resistant.DEG.csv") 
 
