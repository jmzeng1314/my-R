## Skeletal Muscle	GPL9486	Hs133P_Hs_ENST
## different subtype tissue :Glomerulus Tubulus
rm(list=ls())
library(GEOquery)
library(limma) 
GSE18732 <- getGEO('GSE18732', destdir=".",getGPL = F)
exprSet=exprs(GSE18732[[1]])
library("annotate")
GSE18732[[1]]

pdata=pData(GSE18732[[1]])
treatment=factor(unlist(lapply(pdata$title,function(x) strsplit(as.character(x),"[_ ]")[[1]][2]))) 
treatment=relevel(treatment,'normal')
table(treatment)
#treatment
#           normal          diabetic glucoseIntolerant 
#               47                45                26 
platformDB='hgu133plus2.db'
library(platformDB, character.only=TRUE)
probeset <- featureNames(GSE18732[[1]])
ENST=sub("_at","",probeset)
library(org.Hs.eg.db)
ENST2ENG=toTable(org.Hs.egENSEMBLTRANS)
EGID=ENST2ENG[match(ENST,ENST2ENG[,2]),1]

a=cbind(EGID,exprSet)
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
 
design=model.matrix(~ treatment)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
vennDiagram(decideTests(fit))

write.csv(topTable(fit,coef='treatmentdiabetic',n=Inf,adjust='BH'),"diabetic-normal.DEG.csv")
write.csv(topTable(fit,coef='treatmentglucoseIntolerant',n=Inf,adjust='BH'),"glucoseIntolerant-normal.DEG.csv") 
