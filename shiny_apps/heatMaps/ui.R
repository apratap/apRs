library(shiny)

#main UI code
shinyUI(pageWithSidebar(
  
  #Application title
  headerPanel("Pathway Analysis : PCBC Expression Data "),

  #SIDE BAR PANEL FOR USER OPTIONS
  sidebarPanel(
#     selectInput("user_selected_pathway",
#                 "Select a pathway:",
#                 choices = names(KEGG_available_pathways) #loaded from getDATA.R
#     ),
#     br(),
#   
    selectInput("user_selected_geneList",
                  "Significant gene lists",
                  choices = sort(names(precomputed_enrichedPathways_in_geneLists)) #loaded from getDATA.R
    ),
    
    #list of enriched pathways
    #dynamically updated based on used selected gene list
#     selectInput("enrichedPathways",
#                 "Enriched pathway/s",
#                 choices = 'select a significant gene list'
#     ),

    
    br(),
    br(),
    
    #dynamic UI
    uiOutput("enrichedPathways"),
    
    br(),
    br(),
    #FILTER OPTIONS
    h4('Custom Filters:'),
    
    #1. filter based on group level 1 differentiation state
    checkboxGroupInput('level_1_diff_state', h5('Level 1 differentiation state'),choices=level_1_diff_state),
    br(),
    
    #2. filter based on group level 3 differentiation state
    checkboxGroupInput('level_3_diff_state',h5('Level 3 differentiation state'),choices=level_3_diff_state),
    br(),
    
    #3. filter based on cell origin
    checkboxGroupInput('cell_origin',h5('Cell type of origin'),choices=cell_origin),
    br(),
    
    
    #. Filter based on user defined search list
    h5('Search on a custom gene list:'),
    tags$textarea(id="custom_gene_list",rows=5,cols=60,sample_gene_list),
    #helpText("Note: Select <Custom gene list> in the pathway list above"),
    
    #Download fitlered data
    downloadButton('downloadData','Download Filtered Exp Data')
  ),
  
  
  #Main shiny panel
  mainPanel(
    plotOutput("heatMap",height="700px"),
    br(),
    tableOutput("summary"),
    br()
    )
))
