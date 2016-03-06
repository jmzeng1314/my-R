#Introduction for Differential Expression Analysis(DEG):
    Differential analysis, also known as marker selection, is the search for genes that are differentially expressed in distinct phenotypes.
#about this APP
    Run an Differential analysis by preprocessing gene expression data and visualizing the resulting data !
#We offer 5 bioconductor packages to do DEG
    *lima(trend/notrend/voom)
    *deseq/deseq2
    *edgeR(classic/gml)
    *bayseq
#we also offer 3 examples which come from bioconductor
    *CLL(sCLLex)  expressionset object  exprs(eset)
    *pasilla(pasillaGenes)  CountDataSet object
    *airway(airway)
#Then you can upload the data you need analysis by yourself !
So far,we just accept 2 files in certain formats which are GCT and CLS !

    *GCT file   Gene expression data 
    *CLS file   The class of each sample 

#you can look through the result from 4 aspects
    *results table()
    *sample clustering 
    *MA plot
    *volcano plot