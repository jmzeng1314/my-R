## NASH	liver	GPL15207	Affymetrix Human Gene Expression Array
##  Patients were grouped according histology: 
##  group I (<5% steatosis), group II (NAFLD, 30-50% steatosis) and group III (NASH).
rm(list=ls())
library(GEOquery)
library(limma) 
GSE59045 <- getGEO('GSE59045', destdir=".",getGPL = F)
exprSet=exprs(GSE59045[[1]])
library("annotate")
GSE59045[[1]]

pdata=pData(GSE59045[[1]])
treatment=factor(c(rep('steatosis',6),rep('NAFLD',4),rep('NASH',5)) ) 
treatment=relevel(treatment,'steatosis')

#platformDB='hgu133plus2.db'
#library(platformDB, character.only=TRUE)
probeset <- featureNames(GSE59045[[1]])
library(Biobase)
library(GEOquery)
#Download GPL file, put it in the current directory, and load it:
gpl <- getGEO('GPL15207', destdir=".")
colnames(Table(gpl)) ## [1] 41108    17
head(Table(gpl)[,c(1,15,19)])  ## you need to check this , which column do you need 
probe2symbol=Table(gpl)[,c(1,15)]
 
symbol=probe2symbol[match(probeset,probe2symbol[,1]),2]
symbol= as.character(symbol)
exprSet=exprs(GSE59045[[1]])
a=cbind(symbol,exprSet)
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
vennDiagram(decideTests(fit))
NASH_steatosis=topTable(fit,coef='treatmentNASH',n=Inf,adjust='BH')
NAFLD_steatosis=topTable(fit,coef='treatmentNAFLD',n=Inf,adjust='BH')
write.csv(topTable(fit,coef='treatmentNASH',n=Inf,adjust='BH'),"NASH-steatosis.DEG.csv")
write.csv(topTable(fit,coef='treatmentNAFLD',n=Inf,adjust='BH'),"NAFLD-steatosis.DEG.csv") 
