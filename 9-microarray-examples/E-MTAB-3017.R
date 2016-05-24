## GPL6480//Agilent//ArrayExpress
getwd()
# http://www.bioconductor.org/packages/release/bioc/vignettes/ArrayExpress/inst/doc/ArrayExpress.pdf
# http://www.bioconductor.org/packages/release/bioc/manuals/Biobase/man/Biobase.pdf
# http://bioinformatics.oxfordjournals.org/content/25/16/2092.full.pdf
# https://www.bioconductor.org/packages/3.3/bioc/vignettes/gCMAP/inst/doc/diffExprAnalysis.pdf
library(ArrayExpress)
## getAE download all raw data for this study
## ArrayExpress just donwload the expressionSet 
rawset = ArrayExpress("E-MTAB-3017")    
# http://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-3017/
## it's not a AffyBatch object , but a NChannelSet  for this study .
#exprSet=exprs(rawset)
#Investigation description       E-MTAB-3017.idf.txt
#Sample and data relationship    E-MTAB-3017.sdrf.txt
#Raw data (1)                    E-MTAB-3017.raw.1.zip
#Array design                    A-AGIL-28.adf.txt
## A-AGIL-28 - Agilent Whole Human Genome Microarray 4x44K 014850 G4112F (85 cols x 532 rows)
#R ExpressionSet                 E-MTAB-3017.eSet.r

##  for this platform the element names: E, Eb 
exprSet_Eb=assayDataElement(rawset,"Eb")
exprSet_E=assayDataElement(rawset,"E")

library(Biobase)
library(GEOquery)
#Download GPL file, put it in the current directory, and load it:
gpl <- getGEO('GPL6480', destdir=".")
colnames(Table(gpl)) ## [1] 41108    17
head(Table(gpl)[,c(1,6,7)])  ## you need to check this , which column do you need 
write.csv(Table(gpl)[,c(1,6,7)],"GPL6400.csv")


pdata=pData(rawset)
group_list=pdata$Factor.Value.disease.

##I don't know how to create a proper dat , so just stop here !!!

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
exprSet=rmDupID(dat)

rn=rownames(exprSet)
exprSet=apply(exprSet,2,as.numeric)
rownames(exprSet)=rn
#exprSet=log(exprSet)
#boxplot(exprSet,las=2)
library(limma) 
design=model.matrix(~factor(group_list))
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
topTable(fit,coef=2,n=20)
write.csv(topTable(fit,coef=2,n=Inf,adjust='BH'),"obesityVSnormal.DEG.csv")
