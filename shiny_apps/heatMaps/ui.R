library(shiny)

#main UI code
shinyUI(pageWithSidebar(
  
  #Application title
  headerPanel("Pathway Analysis : PCBC Expression Data "),

  #SIDE BAR PANEL FOR USER OPTIONS
  sidebarPanel(

    selectInput("user_selected_geneList",
                  "Significant gene lists",
                  choices = sort(names(precomputed_enrichedPathways_in_geneLists)) #loaded from getDATA.R
    ),
       
    br(),
    
    #dynamically updated
    selectInput('enrichedPathways','Enriched Pathways',choices='ALL'),
    
    br(),

    #heatmap annotation labels
    checkboxGroupInput('heatmap_annotation_labels', h5('Heatmap Annotation Labels'),
                       choices  = names(heatmap_annotation_cols),
                       selected = c('Diff Name')),
    br(),
    #FILTER OPTIONS
    h4('Custom Filters:'),
    
    #1. filter based on mod_linetype
    checkboxGroupInput('mod_linetype', h5('Line Type'),choices=sort(unique(metadata$mod_linetype) ) ),
    br(),
    
    #2. filter based on diff_short_name
    checkboxGroupInput('mod_diffnameshort',h5('Differentiation Name'),choices= sort(unique(metadata$mod_diffnameshort) ) ),
    br(),
    
    #3. filter based on cell origin
    checkboxGroupInput('cell_origin',h5('Cell Origin'),choices=sort(unique(metadata$mod_origcell) ) ),
    br(),
    
    #4. filter based on induction genes
    checkboxGroupInput('induction_genes',h5('Induction Genes'),choices=sort(unique(metadata$inductiongenes) ) ),
    br(),
    
    #. Filter based on user defined search list
    h5('Search on a custom gene list:'),
    tags$textarea(id="custom_gene_list",rows=5,cols=60,sample_gene_list),
    #helpText("Note: Select <Custom gene list> in the pathway list above"),
    
    br(),
    #Download fitlered data
    downloadButton('downloadData','Download Filtered Exp Data')
  ),
  
  
  #Main shiny panel
  mainPanel(
    textOutput("msg"),
    plotOutput("heatMap",height="700px"),
    br(),
    tableOutput("summary"),
    br()
    )
))
