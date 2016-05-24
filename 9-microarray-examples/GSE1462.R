
library(GEOquery)
library(limma) 
GSE1462 <- getGEO('GSE1462', destdir=".",getGPL = F)
exprSet=exprs(GSE1462[[1]])
library("annotate")
GSE1462[[1]]
pdata=pData(GSE1462[[1]])
group_list=c(rep('normal',3),rep('MELAS',4),rep('PEO',4),rep('mtDNA_deletion',4))
group_list=factor(group_list)
group_list <- relevel(group_list, ref="normal")

platformDB='hgu133a.db'
library(platformDB, character.only=TRUE)
probeset <- featureNames(GSE1462[[1]])
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
exprSet=log(exprSet)
boxplot(exprSet,las=2)

design=model.matrix(~ group_list)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
vennDiagram(fit,include=c("up","down"))

write.csv(topTable(fit,coef='group_listMELAS',n=Inf,adjust='BH'),"MELAS-normal.DEG.csv")
write.csv(topTable(fit,coef='group_listmtDNA_deletion',n=Inf,adjust='BH'),"mtDNA_deletion-normal.DEG.csv")
write.csv(topTable(fit,coef='group_listPEO',n=Inf,adjust='BH'),"PEO-normal.DEG.csv")
