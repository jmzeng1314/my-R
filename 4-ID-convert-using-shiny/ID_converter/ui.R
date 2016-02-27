library(shiny)

# Define UI for random distribution application 
shinyUI(fluidPage(
    
  # Application title
  titlePanel("ID converter:"),
  
  # Sidebar with controls to select the random distribution type
  # and number of observations to generate. Note the use of the
  # br() element to introduce extra vertical spacing
  sidebarLayout(
    sidebarPanel(
		h5("single gene search:"),
		textInput('id','please input a gene ID',"TP53"),	
		hr(),
		h5("choose the ID type:"),
		selectInput("IDtype",
					label = "choose the ID type",
					choices = c("HUGO Gene Symbol"  = "ID.Symbol",
								"Entrez Gene ID" 	= "ID.Entrez",
								"HUGO Gene Name"    = "ID.Name",
								"ENSEMBL Gene ID"	= "ID.ENSEMBL")),		
		#submitButton("Run"),
		hr(),		
		h5("upload a entrez ID list for batch search:"),
		fileInput('file1', 'Choose txt File',
                accept=c('text/csv', 
								'text/comma-separated-values,text/plain', 
								'.txt')),
		 tags$hr(),
		 checkboxInput('header', 'Header', TRUE)
		 
	),
    
    # Show a tabset that includes a plot, summary, and table view
    # of the generated distribution
    mainPanel(
      tabsetPanel(type = "tabs", 
          tabPanel("single gene ", 
             h3("output for single gene :"),hr(),
             wellPanel(
               h3("Gend details:"),
               hr(),
               column(5,p("Entrez Gene ID")),column(6,htmlOutput('Entrez')),br(),br(),
               column(5,p("HUGO Gene Symbol")),column(6,htmlOutput('Symbol')),br(),br(),
               column(5,p("HUGO Gene Name")),column(6,textOutput('Name')),br(),br(),
               column(5,p("hg19 Gene Locus")),column(6,textOutput('Locus')),br(),br(),
               column(5,p("ENSEMBL Gene ID ")),column(6,htmlOutput('ENSEMBL')),br(),br()
             ),
             tabsetPanel(type = "tabs", 
                         tabPanel("kegg", tableOutput("kegg")), 
                         tabPanel("go_bp", tableOutput("go_bp")),  
                         tabPanel("go_cc", tableOutput("go_cc")), 
                         tabPanel("go_mf", tableOutput("go_mf"))
             )        
           
          ), 
          tabPanel("batch search",
             h3("output for batch search:"),hr(),
             tableOutput('batch_infor')	      
          )
      )  
    )##end for mainPanel
  )
))