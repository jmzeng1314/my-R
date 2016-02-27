library(shiny)
annotationPKG='org.Hs.eg.db'
library(org.Hs.eg.db)
library(GOstats)
library(XML)

BP_result=NULL
CC_result=NULL
MF_result=NULL

shinyServer(function(input, output) {
	output$FileInfo <- renderPrint({
		inFile <- input$file1    
		if (is.null(inFile))
		 return(NULL) 
		inFile
	 })
	output$BP <- renderTable({
	inFile <- input$file1 
	if ( !is.null(input$file1)){
		diff_gene_list=read.table(inFile$datapath,header=F,quote="")
		diff_gene_list=diff_gene_list[,1]
		GO.hyperG.params = new("GOHyperGParams", geneIds=diff_gene_list, universeGeneIds=NULL, annotation=annotationPKG, 
								ontology='BP', pvalueCutoff=1, conditional = FALSE, testDirection = "over")
		GO.hyperG.results = hyperGTest(GO.hyperG.params);	
		BP_result=summary(GO.hyperG.results)
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
			GO.hyperG.params = new("GOHyperGParams", geneIds=diff_gene_list, universeGeneIds=NULL, annotation=annotationPKG, 
									ontology='MF', pvalueCutoff=1, conditional = FALSE, testDirection = "over")
			GO.hyperG.results = hyperGTest(GO.hyperG.params);
			MF_result=summary(GO.hyperG.results)
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
		GO.hyperG.params = new("GOHyperGParams", geneIds=diff_gene_list, universeGeneIds=NULL, annotation=annotationPKG, 
								ontology='CC', pvalueCutoff=1, conditional = FALSE, testDirection = "over")
		GO.hyperG.results = hyperGTest(GO.hyperG.params);	
		CC_result=summary(GO.hyperG.results)	
		P_value=as.numeric(as.character(CC_result[,2]))
		CC_result[P_value<input$decimal,]
	}else{
		return(NULL)
	}

	})  
})