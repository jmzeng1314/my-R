library(shiny)
annotationPKG='org.Hs.eg.db'
library(org.Hs.eg.db)
library(GOstats)
library(XML)
# Define server logic for random distribution application
shinyServer(function(input, output) {
	output$FileInfo <- renderPrint({
		inFile <- input$file1    
		if (is.null(inFile))
		 return(NULL) 
		#if ( !is.null(inFile)){
		if ( F){
			diff_gene_list=read.table(inFile$datapath,header=F,quote="")
			diff_gene_list=diff_gene_list[,1]
			GOES = c('BP','CC', 'MF');
			# first do Go pathway enrichment
			for (ontology in GOES) {
			 GO.hyperG.params = new("GOHyperGParams", geneIds=diff_gene_list, universeGeneIds=NULL, annotation=annotationPKG, 
									ontology=ontology, pvalueCutoff=1, conditional = FALSE, testDirection = "over")
			 GO.hyperG.results = hyperGTest(GO.hyperG.params);
			 outHTMLname=paste("GO_",ontology,".enrichment.html",sep="")
			 htmlReport(GO.hyperG.results, file=outHTMLname, summary.args=list("htmlLinks"=TRUE))
			 #summary(GO.hyperG.results)
			}
			#then do kegg pathway enrichment !
			hyperG.params = new("KEGGHyperGParams", geneIds=diff_gene_list, universeGeneIds=NULL, annotation=annotationPKG, 
								categoryName="KEGG", pvalueCutoff=1, testDirection = "over")
			KEGG.hyperG.results = hyperGTest(hyperG.params); #116个基因，富集到了67个通路！
			htmlReport(KEGG.hyperG.results, file="kegg.enrichment.html", summary.args=list("htmlLinks"=TRUE))		
		}
		inFile
	 })
	# Generate an HTML table view of the data
	BP_result="GO_BP.enrichment.html";
	MF_result="GO_MF.enrichment.html";
	CC_result="GO_CC.enrichment.html";
	KEGG_result="kegg.enrichment.html";
	
	output$BP <- renderTable({
	if (file.exists(BP_result)){
		tmp=readHTMLTable(BP_result)
		tmp=tmp[[1]]
		P_value=as.numeric(as.character(tmp[,2]))
		tmp[P_value<input$decimal,]
	}else{
		return(NULL)
	}

	})
	
	output$KEGG <- renderTable({
	if (file.exists(KEGG_result)){
		tmp=readHTMLTable(KEGG_result)
		tmp=tmp[[1]]
		P_value=as.numeric(as.character(tmp[,2]))
		tmp[P_value<input$decimal,]
	}else{
		return(NULL)
	}
	})
	
	output$MF <- renderTable({
	if (file.exists(MF_result)){
		tmp=readHTMLTable(MF_result)
		tmp=tmp[[1]]
		P_value=as.numeric(as.character(tmp[,2]))
		tmp[P_value<input$decimal,]
	}else{
		return(NULL)
	}
	})
	
	output$CC <- renderTable({
	if (file.exists(CC_result)){
		tmp=readHTMLTable(CC_result)
		tmp=tmp[[1]]
		P_value=as.numeric(as.character(tmp[,2]))
		tmp[P_value<input$decimal,]
	}else{
		return(NULL)
	}
	})  
})