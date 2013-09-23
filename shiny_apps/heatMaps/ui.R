library(shiny)

#main UI code
shinyUI(pageWithSidebar(
  
  #Application title
  headerPanel("Pathway Analysis : PCBC Expression Data "),

  #SIDE BAR PANEL FOR USER OPTIONS
  sidebarPanel(
    selectInput("user_selected_pathway",
                "Select a pathway:",
                choices = names(KEGG_available_pathways) #loaded from available_pathways.R
    ),
    br(),
    br(),
    
    #FILTER OPTIONS
    
    #1. filter based on sex
    checkboxGroupInput('sex','Sex',choices=sex),
    br(),
    
    #2. filter based on sample origin
    checkboxGroupInput('origin','Origin',choices=origin),
    br(),
    
    
    #. Filter 
    h5('Search on a custom gene list:'),
    tags$textarea(id="custom_gene_list",rows=5,cols=60,sample_gene_list),
    helpText("Search your own list of HUGO gene names separated by tab, comma, line or space. Note: Select custom pathway in the pathway list above"),
    
    
    #Download fitlered data
    downloadButton('downloadData','Download Filtered Exp Data')
  ),
  
  
  #Main shiny panel
  mainPanel(
    h3(textOutput("currentPathway")),
    br(),
    tableOutput("summary"),
    br(),
    plotOutput("heatMap"),
    br(),
    br()
    
    )

))
