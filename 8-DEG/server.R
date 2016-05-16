library(shiny)
suppressPackageStartupMessages(library(DESeq))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(baySeq))
suppressPackageStartupMessages(library(pasilla))
suppressPackageStartupMessages(library(pasilla))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(CLL))

#library(airway)
#runApp("8-DEG/")

rmDupID <-function(a=matrix(c(1,1:5,2,2:6,2,3:7),ncol=6)){
  exprSet=a[,-1]
  rowMeans=apply(exprSet,1,mean)
  a=a[order(rowMeans,decreasing=T),]
  return(a[!duplicated(a[,1]),])
}

DEG_edger_classic <- function(exprSet=exprSet,group_list=group_list){
  d <- DGEList(counts=exprSet,group=factor(group_list))
  d.full <- d # keep the old one in case we mess up
  #apply(d$counts, 2, sum) # total gene counts per samplekeep <- rowSums(cpm(d)>100) >= 2
  keep <- rowSums(cpm(d)>100) >= 2
  d <- d[keep,]
  d$samples$lib.size <- colSums(d$counts)
  d <- calcNormFactors(d)
  
  png("MDS.png")
  plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
  legend("bottomleft", as.character(unique(d$samples$group)), col=1:3, pch=20)
  dev.off()
  
  d1 <- estimateCommonDisp(d, verbose=T)
  d1 <- estimateTagwiseDisp(d1)
  
  png("BCV.png")
  plotBCV(d1)
  dev.off()
  
  et12 <- exactTest(d1) 
  
  png("MA.png")
  de1 <- decideTestsDGE(et12, adjust.method="BH", p.value=0.05)
  de1tags12 <- rownames(d1)[as.logical(de1)] 
  plotSmear(et12, de.tags=de1tags12)
  dev.off()
  
  nrDEG=topTags(et12, n=nrow(exprSet))
  nrDEG=as.data.frame(nrDEG)
  #write.csv(nrDEG,"edger_classic.results.csv",quote = F)
  
}

DEG_edger_glm <- function(exprSet=exprSet,group_list=group_list){
  d <- DGEList(counts=exprSet,group=factor(group_list))
  d.full <- d # keep the old one in case we mess up
  #apply(d$counts, 2, sum) # total gene counts per samplekeep <- rowSums(cpm(d)>100) >= 2
  keep <- rowSums(cpm(d)>100) >= 2
  d <- d[keep,]
  d$samples$lib.size <- colSums(d$counts)
  d <- calcNormFactors(d)   
  
  png("MDS.png")
  plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
  legend("bottomleft", as.character(unique(d$samples$group)), col=1:3, pch=20)
  dev.off()
  
  design.mat <- model.matrix(~ 0 + d$samples$group)
  colnames(design.mat) <- levels(d$samples$group)
  d2 <- estimateGLMCommonDisp(d,design.mat)
  d2 <- estimateGLMTrendedDisp(d2,design.mat, method="auto")
  # You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
  # The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
  d2 <- estimateGLMTagwiseDisp(d2,design.mat)
  
  png("BCV.png")
  plotBCV(d2)
  dev.off()
  
  fit <- glmFit(d2, design.mat)
  # compare (group 1 - group 2) to 0:
  # this is equivalent to comparing group 1 to group 2
  lrt12 <- glmLRT(fit)
  
  png("MA.png")
  de2 <- decideTestsDGE(lrt12, adjust.method="BH", p.value = 0.05)
  de2tags12 <- rownames(d2)[as.logical(de2)]
  plotSmear(lrt12, de.tags=de2tags12)
  abline(h = c(-2, 2), col = "blue")
  dev.off()
  
  nrDEG=topTags(lrt12, n=nrow(exprSet))
  nrDEG=as.data.frame(nrDEG)
  #write.csv(nrDEG,"edger_glm.results.csv",quote = F)
}

DEG_limma <- function(exprSet=exprSet,group_list=group_list){
	suppressMessages(library(limma))
	design <- model.matrix(~0+factor(group_list))
	colnames(design)=levels(factor(group_list))
	rownames(design)=colnames(exprSet)
	fit <- lmFit(exprSet,design)
	png("MA.png")
	limma::plotMA(fit)
	abline(0,0,col="blue")
	dev.off()
	fit2 <- eBayes(fit)  ## default no trend !!!
	##eBayes() with trend=TRUE
	tempOutput = topTable(fit2, coef=1, n=Inf)
	nrDEG2 = na.omit(tempOutput) 
	#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
	nrDEG2
}

DEG_limma_trend <- function(exprSet=exprSet,group_list=group_list){
	suppressMessages(library(limma))
	design <- model.matrix(~0+factor(group_list))
	colnames(design)=levels(factor(group_list))
	rownames(design)=colnames(exprSet)
	fit <- lmFit(exprSet,design)
	png("MA.png")
	limma::plotMA(fit)
	abline(0,0,col="blue")
	dev.off()
	fit2 <- eBayes(fit,trend=TRUE)  ## default no trend !!!
	tempOutput = topTable(fit2, coef=1, n=Inf)
	nrDEG2 = na.omit(tempOutput) 
	write.csv(nrDEG2,"limma_trend.results.csv",quote = F)
}

DEG_voom <- function(exprSet=exprSet,group_list=group_list){
	suppressMessages(library(limma))
	design <- model.matrix(~0+factor(group_list))
	colnames(design)=levels(factor(group_list))
	rownames(design)=colnames(exprSet)
	v <- voom(exprSet,design,normalize="quantile")
	png("RAWvsNORM.png")
	exprSet_new=v$E
	par(cex = 0.7)
	n.sample=ncol(exprSet)
	if(n.sample>40) par(cex = 0.5)
	cols <- rainbow(n.sample*1.2)
	par(mfrow=c(2,2))
	boxplot(exprSet, col = cols,main="expression value",las=2)
	boxplot(exprSet_new, col = cols,main="expression value",las=2)
	hist(exprSet)
	hist(exprSet_new)
	dev.off()
	
	png("RAWvsNORM.png")
	rld <- rlogTransformation(dds)
	exprSet_new=assay(rld)
	par(cex = 0.7)
	n.sample=ncol(exprSet)
	if(n.sample>40) par(cex = 0.5)
	cols <- rainbow(n.sample*1.2)
	par(mfrow=c(2,2))
	boxplot(exprSet, col = cols,main="expression value",las=2)
	boxplot(exprSet_new, col = cols,main="expression value",las=2)
	hist(exprSet)
	hist(exprSet_new)
	dev.off()
	
	png("MDS.png")
	plotMDS(v, labels=1:ncol(exprSet),col=rainbow(ncol(exprSet)))
	legend("topright", legend=unique(group_list), col=1:2, pch=15)
	dev.off()
	fit <- lmFit(v,design)
	png("MA.png")
	limma::plotMA(fit)
	#abline(0,0,col="blue")
	dev.off()

	fit2 <- eBayes(fit)
	tempOutput = topTable(fit2, coef=1, n=Inf)
	nrDEG2 = na.omit(tempOutput)
	#write.csv(nrDEG2,"voom.results.csv",quote = F)
}

DEG_DESeq2 <- function(exprSet=exprSet,group_list=group_list){
	suppressMessages(library(DESeq2))
	exprSet=ceiling(exprSet)
	(colData <- data.frame(row.names=colnames(exprSet), group_list=group_list))
	dds <- DESeqDataSetFromMatrix(countData = exprSet,
								colData = colData,
								design = ~ group_list)
	dds <- DESeq(dds)
	png("qc_dispersions.png", 1000, 1000, pointsize=20)
	plotDispEsts(dds, main="Dispersion plot")
	dev.off()
	res <- results(dds)
	
	png("RAWvsNORM.png")
	rld <- rlogTransformation(dds)
	exprSet_new=assay(rld)
	par(cex = 0.7)
	n.sample=ncol(exprSet)
	if(n.sample>40) par(cex = 0.5)
	cols <- rainbow(n.sample*1.2)
	par(mfrow=c(2,2))
	boxplot(exprSet, col = cols,main="expression value",las=2)
	boxplot(exprSet_new, col = cols,main="expression value",las=2)
	hist(exprSet)
	hist(exprSet_new)
	dev.off()
	
	library(RColorBrewer)
	(mycols <- brewer.pal(8, "Dark2")[1:length(unique(group_list))])

	# Sample distance heatmap
	sampleDists <- as.matrix(dist(t(exprSet_new)))
	#install.packages("gplots",repos = "http://cran.us.r-project.org")
	library(gplots)
	png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
	heatmap.2(as.matrix(sampleDists), key=F, trace="none",
			  col=colorpanel(100, "black", "white"),
			  ColSideColors=mycols[group_list], RowSideColors=mycols[group_list],
			  margin=c(10, 10), main="Sample Distance Matrix")
	dev.off()

	png("MA.png")
	DESeq2::plotMA(res, main="DESeq2", ylim=c(-2,2))
	dev.off()
	resOrdered <- res[order(res$padj),]
	resOrdered=as.data.frame(resOrdered)
	#write.csv(resOrdered,"deseq2.results.csv",quote = F)
}

DEG_DESeq <- function(exprSet=exprSet,group_list=group_list){
	suppressMessages(library(DESeq))
	exprSet=ceiling(exprSet)
	de=newCountDataSet(exprSet,group_list)
	de=estimateSizeFactors(de)
	de=estimateDispersions(de)
	png("qc_dispersions.png", 1000, 1000, pointsize=20)
	DESeq::plotDispEsts(de, main="Dispersion plot")
	dev.off()
	
	res=nbinomTest(de,unique(group_list)[1],unique(group_list)[2])
	png("MA.png")
	DESeq::plotMA(res)
	dev.off()
	rownames(res)=res[,1]
	res=res[,-1]
	resOrdered <- res[order(res$padj),]
	#write.csv(resOrdered,"deseq.results.csv",quote = F)
}

DEG_baySeq<- function(exprSet=exprSet,group_list=group_list){
	suppressMessages(library(baySeq))
	suppressMessages(library(snow))
	cl <- makeCluster(4, "SOCK")
	groups=list(NDE = rep(1,length(group_list)),
				DE = as.numeric(factor(group_list)) 
	)
	CD <- new("countData", data = exprSet, replicates = group_list, groups = groups )
	libsizes(CD) <- getLibsizes(CD)
	png("MA.png")
	 plotMA.CD(CD, samplesA = unique(group_list)[1], samplesB = unique(group_list)[2],
           col = c(rep("red", 100), rep("black", 900)))
	dev.off()
	CD <- getPriors.NB(CD, samplesize = 1000, estimation = "QL",cl=cl)
	CD <- getLikelihoods.NB(CD, pET = 'BIC', cl = cl)
	res=topCounts(CD, group = "DE")
	res=topCounts(CD, group = "DE",number = nrow(exprSet)) 
	png("Posteriors_distribution.png")
	plotPosteriors(CD, group = "DE", col = c(rep("red", 100), rep("black", 900)))
	dev.off()
	#write.csv(resOrdered,"deseq.results.csv",quote = F)
}


#data(airway)
shinyServer(function(input, output) {
	get_example_data <- reactive({
	inFile1 <- input$gct_file
	inFile2 <- input$cls_file
    if (is.null(inFile1) || is.null(inFile2) ){
		if(input$example == 'CLL'){
			data(sCLLex)
			exprSet=exprs(sCLLex)
			samples=sampleNames(sCLLex)
			pdata=pData(sCLLex)
			group_list=as.character(pdata[,2])
			
		}else if(input$example == 'pasilla'){
			data(pasillaGenes)
			exprSet=counts(pasillaGenes)
			samples=sampleNames(pasillaGenes)
			pdata=pData(pasillaGenes)[,c("condition","type")]
			group_list=as.character(pdata[,1])
			
		}else{
		
		}

	}else{
		a=read.table(inFile1$datapath,stringsAsFactors=F,skip=2,header=T,sep="\t")
		## probeID/geneID/value1/value2~~~~~~
		if (length(unique(a[,2])) <3 ){
			exprSet=a[,-2]
		}else{
			exprSet=a[,-1]
		}
		exprSet=rmDupID(exprSet)
		rownames(exprSet)=exprSet[,1]
		exprSet=exprSet[,-1]
		sample_list=colnames(exprSet)
		b=read.table(inFile2$datapath,stringsAsFactors=F,skip=2)
		group_list=as.character(b)
		sample_group=data.frame(sample_list,group_list)
		names(sample_group)=c('sample_id','compare')
		pdata=sample_group
		rownames(pdata)=sample_list
	}
	matrix_a=exprSet[,rownames(pdata)[group_list==unique(group_list)[1]]]
	matrix_b=exprSet[,rownames(pdata)[group_list==unique(group_list)[2]]]	
	list(exprSet=exprSet,pdata=pdata,group_list=group_list,matrix_a=matrix_a,matrix_b=matrix_b)
	})

	output$experimental_information <- renderText({
		example_data=get_example_data()
		paste("there are ",ncol(example_data$exprSet)," samples","\n",
				ncol(example_data$matrix_a),"are group ",unique(example_data$group_list)[1],"\n",
				ncol(example_data$matrix_b),"are group ",unique(example_data$group_list)[2],"\n",
				"there are totally ",nrow(example_data$exprSet)," genes/probes "
		)
		
	})
	
	output$group_info  <- renderDataTable({
		example_data=get_example_data()
		example_data$pdata		
	})
	output$matrix_control  <- renderDataTable({
		example_data=get_example_data()
		tmp=example_data$matrix_a		
		tmp=cbind(ID=rownames(tmp),tmp)
		tmp
	})
	output$matrix_case  <- renderDataTable({
		example_data=get_example_data()
		tmp=example_data$matrix_b	
		tmp=cbind(ID=rownames(tmp),tmp)
		tmp		
	})
	output$DEG_results  <- renderDataTable({
		example_data=get_example_data()
		DEG_method=tolower(input$method)
		exprSet=example_data$exprSet
		group_list=example_data$group_list
		if(DEG_method=='limma'){
			DEG_results=DEG_limma(exprSet,group_list)		
		}else if (DEG_method=='voom'){
			DEG_results=DEG_voom(exprSet,group_list)
		}else if (DEG_method=='deseq'){
			DEG_results=DEG_DESeq(exprSet,group_list)
		}else if (DEG_method=='deseq2'){
			DEG_results=DEG_DESeq2(exprSet,group_list)
		}else if (DEG_method=='edger_classic'){
			DEG_results=DEG_edger_classic(exprSet,group_list)
		}else if (DEG_method=='edger_glm'){
			DEG_results=DEG_edger_glm(exprSet,group_list)
		}else if (DEG_method=='bayseq'){
			DEG_results=DEG_baySeq(exprSet,group_list)
		}else{
			stop(sprintf("we don't offer the method you want"))
		}	
		DEG_results=cbind(symbol=rownames(DEG_results),DEG_results)
		DEG_results
	})
	output$MA_plot <-renderImage({
	 list(src = "MA.png",
         contentType = 'image/png',
         width = 400,
         height = 300,
         alt = "This is alternate text")
	})
	output$MDS_plot <-renderImage({
	list(src = "MDS.png",
         contentType = 'image/png',
         width = 400,
         height = 300,
         alt = "This is alternate text")
	})
	output$BCV_plot <-renderImage({
	list(src = "BCV.png",
         contentType = 'image/png',
         width = 400,
         height = 300,
         alt = "This is alternate text")
	})
	output$volcanoplot <-renderImage({
	list(src = "volcanoplot.png",
         contentType = 'image/png',
         width = 400,
         height = 300,
         alt = "This is alternate text")
	})
	output$qc_dispersions <-renderImage({
	list(src = "qc_dispersions.png",
         contentType = 'image/png',
         width = 400,
         height = 300,
         alt = "This is alternate text")
	})
	output$RAWvsNORM <-renderImage({
	list(src = "RAWvsNORM.png",
         contentType = 'image/png',
         width = 400,
         height = 300,
         alt = "This is alternate text")
	})
		
	
})
