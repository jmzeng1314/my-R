rm(list=ls())
library(lumi)
studyID = 'GSE30669'
setwd('G:/array/GSE30669')
# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE30nnn/GSE30669/suppl/GSE30669_HEK_Sample_Probe_Profile.txt.gz
fileName <- 'GSE30669_HEK_Sample_Probe_Profile.txt'  
x.lumi <- lumiR.batch(fileName) ##, sampleInfoFile='sampleInfo.txt')
## it make me so sad that I can't find the sampleInfo.txt
pData(phenoData(x.lumi))  

## Do all the default preprocessing in one step
lumi.N.Q <- lumiExpresso(x.lumi) 
### retrieve normalized data
dataMatrix <- exprs(lumi.N.Q)
##　To speed up the processing and reduce false positives, remove the unexpressed　and un-annotated genes
presentCount <- detectionCall(x.lumi)
selectDataMatrix <- dataMatrix[presentCount > 0,]
probe2target=pData(featureData(x.lumi)) 
## estimate the detect call (percentage of expressed genes) of each sample
temp <- detectionCall(x.lumi, type='sample')
print(temp)
## estimate the present count of each gene (probe)
temp <- detectionCall(x.lumi, type='probe')
hist(temp)

## QC, sampleRelation.png is especially important.
QClumi <- function(example.lumi){
  width <- 800;height <- 800;
  summary(example.lumi, 'QC')
  png('1_density.png',width = width, height = height)
  plot(example.lumi, what='density')	## plot the density
  dev.off()
  png('2_cdf.png',width = width, height = height)
  plotCDF(example.lumi) 
  dev.off()
  png('3_pairwise.png',width = width*2, height = height*2)
  plot(example.lumi, what='pair')		## pairwise plots
  #plot(example.lumi, what='pair', smoothScatter=T)
  dev.off()
  png('4_pairwiseMAplot.png',width = width*2, height = height*2)
  plot(example.lumi, what='MAplot')	
  #plot(example.lumi, what='MAplot', smoothScatter=T)	
  dev.off()
  png('5_density_plot_of_coefficient_of_varience.png',width = width, height = height)
  plot(example.lumi, what='cv')
  dev.off()
  png('6_sampleRelation.png',width = width, height = height)
  plot(example.lumi, what='sampleRelation')
  dev.off()
  png('7_MDS.png',width = width, height = height)
  plot(example.lumi, what='sampleRelation', method='mds' )
  dev.off()
}
QClumi(x.lumi)


# library(GEOquery)
# library(limma)
# GSE30669 <- getGEO('GSE30669', destdir='.',getGPL = F)
# exprSet=exprs(GSE30669[[1]])
# GSE30669[[1]]
# pdata=pData(GSE30669[[1]])
# exprSet=exprs(GSE30669[[1]])


###################################################
### code chunk number 31: Identify differentially expressed genes
###################################################
## Specify the sample type
# It's wrorng to just arrange them in order. 
# sampleType <- factor(c(rep('PMN',3),rep('PDK1',3),rep('MYC',3),rep('E545K',3)))
# limmaArg='PDK1-PMN,MYC-PMN,E545K-PMN'
sampleType <- factor( c('G1','G4','G3','G2','G1','G4','G3','G2','G1','G2','G3','G4' )  )
limmaArg='G2-G1,G3-G1,G4-G1,G4-G2,G3-G2,G4-G3'
## we should check the sampleRelation.png to group them properly, but we could just give them a anonymous label.
if (require(limma)) { 
  design <- model.matrix(~0+  sampleType )
  colnames(design) <-  levels(sampleType)
  fit <- lmFit(selectDataMatrix, design) 

  contrastsCommand=unlist(strsplit(limmaArg, split=","))
  cont.matrix <- makeContrasts(contrasts=contrastsCommand, levels=design)
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)
  options(digits = 4) 
  for(i in 1:length(contrastsCommand)){
    tempOutFile <- paste(studyID,".diffexp.", contrastsCommand[i],".csv", sep="")
    tempOutput = topTable(fit2, coef=i, n=Inf)
    tempOutput = na.omit(tempOutput)
    #### change probeID to geneSymbol according to the probe2target  !!!
    tempOutput$geneSymbol=probe2target[match(rownames(tempOutput),probe2target[,1]),2] 
    write.csv(tempOutput,tempOutFile,quote=FALSE,row.names = F)
  }
} 

## Significant analysis of microarray identified 1,750, 1,080, and 297 differentially expressed genes 
## in these transformed cells when compared with nontransformed control cells, respectively 
##  false discovery rate < 0.05; P < 0.01; 

for(i in 1:length(contrastsCommand)){ 
  tempOutput = topTable(fit2, coef=i, n=Inf) 
  print(nrow(tempOutput[tempOutput$adj.P.Val<0.05,]))
}

 