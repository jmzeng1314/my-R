## GPL96   HG-U133A 
library(GEOquery)
library(limma) 
GSE474 <- getGEO('GSE474', destdir=".",getGPL = F)
exprSet=exprs(GSE474[[1]])
library("annotate")
GSE474[[1]]
pdata=pData(GSE474[[1]])
status=factor(unlist(lapply(pdata$title,function(x) strsplit(as.character(x),"-")[[1]][3])))
status=relevel(status,'NonObese') 

platformDB='hgu133a.db'
library(platformDB, character.only=TRUE)
probeset <- featureNames(GSE474[[1]])
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
exprSet=log(exprSet) ## based on e 
boxplot(exprSet,las=2)

design=model.matrix(~ status)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
vennDiagram(fit,include=c("up","down"))

write.csv(topTable(fit,coef='statusMObese',n=Inf,adjust='BH'),"MObese-NonObese.DEG.csv")
write.csv(topTable(fit,coef='statusObese',n=Inf,adjust='BH')  ,"Obese-NonObese.DEG.csv") 
 
