	library(shiny)

# Define UI for random distribution application 
shinyUI(fluidPage(
    
  # Application title
  titlePanel("Tabsets"),
  
  # Sidebar with controls to select the random distribution type
  # and number of observations to generate. Note the use of the
  # br() element to introduce extra vertical spacing
  sidebarLayout(
    sidebarPanel(
	 fileInput('file1', 'Choose text File'),
	 tags$hr(),
      br(),
	 sliderInput("decimal", "the cutoff for P-value:",
					min = 0, max = 1, value = 0.05, step= 0.01),
	 tags$hr(),
      br()
    ),
    # Show a tabset that includes three catelogy for the GO pathway !
    mainPanel(
	 verbatimTextOutput("FileInfo"),
      tabsetPanel(type = "tabs", 
		tabPanel("KEGG", tableOutput("KEGG")),
        tabPanel("BP", tableOutput("BP")), 
        tabPanel("CC", tableOutput("CC")), 
        tabPanel("MF", tableOutput("MF"))
      )
    )
  )
))

