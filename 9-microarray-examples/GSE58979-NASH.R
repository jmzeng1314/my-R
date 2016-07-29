## 	Homo sapiens/NASH/expression array // 
## http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58979
## GPL15207	[PrimeView] Affymetrix Human Gene Expression Array
## the number of 'visceral fat samples' for stages I-II-III-IV is 10-7-8-5.
## Patients were grouped according histology: group I (<5% steatosis), group II (NAFLD, 30-50% steatosis), group III (NASH) and group IV (NASH + fibrosis F2-F3).
## steatosis//NAFLD//NASH//NASHfibrosis
library(GEOquery)
library(limma) 
GSE58979 <- getGEO('GSE58979', destdir=".",getGPL = F)
exprSet=exprs(GSE58979[[1]])
library("annotate")
GSE58979[[1]]

# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 49395 features, 53 samples 
# element names: exprs 
# protocolData: none
# phenoData
# sampleNames: GSM1423267 GSM1423268 ... GSM1423319 (53 total)
# varLabels: title geo_accession ... data_row_count (38 total)
# varMetadata: labelDescription
# featureData: none
# experimentData: use 'experimentData(object)'
# Annotation: GPL15207 

probeset <- featureNames(GSE58979[[1]])
library(Biobase)
library(GEOquery)
#Download GPL file, put it in the current directory, and load it:
gpl <- getGEO('GPL15207', destdir=".")
colnames(Table(gpl)) ## [1] 49395    24
head(Table(gpl)[,c(1,15,19)])  ## you need to check this , which column do you need 
probe2symbol=Table(gpl)[,c(1,15)]

symbol=probe2symbol[match(probeset,probe2symbol[,1]),2]
symbol= as.character(symbol)
exprSet=exprs(GSE58979[[1]])

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

#exprSet=log(exprSet)
#boxplot(exprSet,las=2)

pdata=pData(GSE58979[[1]])
## check the metadat and find the correct group information
# we must tell R that group should be interpreted as a factor !
HistologyClass=factor( unlist(lapply(pdata$source_name_ch1,function(x) strsplit(as.character(x)," ")[[1]][5]))  )
 
design=model.matrix(~HistologyClass)
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
topTable(fit,coef=2,n=20)
options(digits = 4)
vennDiagram(fit,include=c("up","down"))
topTable(fit,coef=2,adjust='BH',n=20)     
topTable(fit,coef=3,adjust='BH',n=20)  
topTable(fit,coef=4,adjust='BH',n=20)  
## steatosis//NAFLD//NASH//NASHfibrosis
write.csv(topTable(fit,coef=2,adjust='BH',n=Inf),'NAFLD-steatosis.csv') 
write.csv(topTable(fit,coef=3,adjust='BH',n=Inf),'NASH-steatosis.csv') 
write.csv(topTable(fit,coef=4,adjust='BH',n=Inf),'NASHfibrosis-steatosis.csv') 
