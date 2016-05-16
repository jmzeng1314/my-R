library(shiny)

# Define UI for application that do differential expression analysis 
# By different statistical method, which are DESeq,DESeq2,edgeR,baySeq,limma.
# reflinks:http://genepattern.broadinstitute.org/gp/pages/protocols/DiffExp.html
# http://genepattern.broadinstitute.org/gp/pages/protocols/ClassDiscovery.html
# http://www.broadinstitute.org/cancer/software/genepattern/gene-expression-analysis
shinyUI(fluidPage(

  # Application title
  titlePanel("differential expression analysis!"),
  hr(),
  if(F){
	wellPanel(
	  p("Run an Differential analysis by preprocessing gene expression data and visualizing the resulting data !"),
	  p("Differential analysis, also known as marker selection, is the search for genes that are differentially expressed in distinct phenotypes."),
	  hr(),
	  wellPanel(
		p("Differential Expression Analysis---Find genes that are significantly differentially expressed between classes of samples."),
		p("Clustering---Group genes and/or samples by similar expression profiles. "),
		p("Heat Map Viewer shows you differential expression by displaying gene expression values in a heat map format.")
	  ),
	  hr()
	)
  },
  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
	  p("we provide some examples to let you familar with the differential expression analysis!"),
	  p("you can choose one of them to test the analysis method !!"),
	  selectInput("example","choose a example dataset from Bioconductor:",
					c("CLL","pasilla","airway")
	  ),
	  p("you can choose one of these analysis method !!"),
	  selectInput("method","choose a analysis method:",
					c("limma","edgeR_classic","edgeR_glm","DESeq","DESeq2","voom","baySeq")
	  ),
	  ##https://www.bioconductor.org/packages/3.3/data/experiment/
	  hr(),
	  wellPanel(
	  	p("Then you can upload the data you need analysis by yourself !"),
		p("So far,we just accept 2 files in certain formats which are GCT and CLS !"),
		hr(),
		fileInput("gct_file","choose a GCT file for the Gene expression data"),
		fileInput("cls_file","choose a CLS file for the class of each sample"),
		hr(),
		p("Gene expression data must be a ",strong("GCT")," file, and Example file: ",a("all_aml_test.gct",href="ftp://ftp.broadinstitute.org/pub/genepattern/datasets/all_aml/all_aml_test.gct")),
		p("The class of each sample must be identified in a ",strong("CLS")," file, and Example file: ",a(" all_aml_test.cls",href="ftp://ftp.broadinstitute.org/pub/genepattern/datasets/all_aml/all_aml_test.cls")),
		hr()
	  )
    ), ##end for sidebarPanel

    # Show a plot of the generated distribution
    mainPanel(
		wellPanel(
			p(),
			verbatimTextOutput("experimental_information"),
			hr()
		),
		tabsetPanel("results",
			tabPanel("group information",dataTableOutput("group_info")),
			tabPanel("group 1",dataTableOutput("matrix_control")),
			tabPanel("group 2",dataTableOutput("matrix_case")),
			tabPanel("results for the DEG",dataTableOutput("DEG_results")),
			tabPanel("quality control",imageOutput("MA_plot"),
										imageOutput("MDS_plot"),
										imageOutput("BCV_plot"),
										imageOutput("RAWvsNORM"),
										imageOutput("qc_dispersions"),
										imageOutput("volcanoplot")
										
			)
		)
    ) ## end for mainPanel 
  )### end for sidebarLayout
))
