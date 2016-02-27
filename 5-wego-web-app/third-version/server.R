library(shiny)
annotationPKG='org.Hs.eg.db'
library(org.Hs.eg.db)
library(GOstats)
tmp=toTable(org.Hs.egPATH)
GeneID2kegg=tapply(tmp[,2],as.factor(tmp[,1]),function(x) x)
kegg2GeneID=tapply(tmp[,1],as.factor(tmp[,2]),function(x) x)
tmp=toTable(org.Hs.egGO)
GO=split(tmp,tmp[,4])
GeneID2GO_BP=tapply(GO$BP[,2],as.factor(GO$BP[,1]),function(x) x)
GO_BP2GeneID=tapply(GO$BP[,1],as.factor(GO$BP[,2]),function(x) x)
GeneID2GO_CC=tapply(GO$CC[,2],as.factor(GO$CC[,1]),function(x) x)
GO_CC2GeneID=tapply(GO$CC[,1],as.factor(GO$CC[,2]),function(x) x)
GeneID2GO_MF=tapply(GO$MF[,2],as.factor(GO$MF[,1]),function(x) x)
GO_MF2GeneID=tapply(GO$MF[,1],as.factor(GO$MF[,2]),function(x) x)
#diff_gene_list=read.table("diff_gene_list.txt",header=F,quote="")
#diff_gene_list=diff_gene_list[,1]
hyperGtest <- function(GeneID2Path,Path2GeneID,diff_gene){
	diff_gene_has_path=intersect(diff_gene,names(GeneID2Path))
	n=length(diff_gene_has_path) #306
	N=length(GeneID2Path) #5870
	options(digits = 4)
	results=c()
	for (i in names(Path2GeneID)){
	 M=length(Path2GeneID[[i]]) 
	 exp_count=n*M/N
	 k=0
	 for (j in diff_gene_has_path){
		if (i %in% GeneID2Path[[j]]) k=k+1
	 }
	 OddsRatio=k/exp_count
	 if (k==0) next
	 p=phyper(k-1,M, N-M, n, lower.tail=F)
	 print(paste(i,p,OddsRatio,exp_count,k,M,sep="    "))
	 results=rbind(results,c(i,p,OddsRatio,exp_count,k,M))
	}
	colnames(results)=c("PathwayID","Pvalue","OddsRatio","ExpCount","Count","Size")
	return(results)
}

#kegg_result=hyperGtest(GeneID2kegg,kegg2GeneID,diff_gene_list)
shinyServer(function(input, output) {
	output$FileInfo <- renderPrint({
		inFile <- input$file1    
		if (is.null(inFile))
		 return(NULL) 
		inFile
	 })
	output$KEGG <- renderTable({
	inFile <- input$file1 
	if ( !is.null(input$file1)){
		diff_gene_list=read.table(inFile$datapath,header=F,quote="")
		diff_gene_list=diff_gene_list[,1]
		Kegg_result=hyperGtest(GeneID2kegg,kegg2GeneID,diff_gene_list)
		P_value=as.numeric(as.character(Kegg_result[,2]))
		Kegg_result[P_value<input$decimal,]
	}else{
		return(NULL)
	}
	})
	output$BP <- renderTable({
	inFile <- input$file1 
	if ( !is.null(input$file1)){
		diff_gene_list=read.table(inFile$datapath,header=F,quote="")
		diff_gene_list=diff_gene_list[,1]	
		BP_result=hyperGtest(GeneID2GO_BP,GO_BP2GeneID,diff_gene_list)
		P_value=as.numeric(as.character(BP_result[,2]))
		BP_result[P_value<input$decimal,]
	}else{
		return(NULL)
	}
	})
	output$MF <- renderTable({
	inFile <- input$file1 
		if ( !is.null(input$file1)){
			diff_gene_list=read.table(inFile$datapath,header=F,quote="")
			diff_gene_list=diff_gene_list[,1]	 
			MF_result=hyperGtest(GeneID2GO_MF,GO_MF2GeneID,diff_gene_list)
			P_value=as.numeric(as.character(MF_result[,2]))
			MF_result[P_value<input$decimal,]
		}else{
			return(NULL)
		}
	})
	output$CC <- renderTable({
	inFile <- input$file1 
	if ( !is.null(input$file1)){
		diff_gene_list=read.table(inFile$datapath,header=F,quote="")
		diff_gene_list=diff_gene_list[,1]	
		CC_result=hyperGtest(GeneID2GO_CC,GO_CC2GeneID,diff_gene_list)	
		P_value=as.numeric(as.character(CC_result[,2]))
		CC_result[P_value<input$decimal,]
	}else{
		return(NULL)
	}

	})  
})