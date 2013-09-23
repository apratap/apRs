#load the modules
library("shiny")


#load local files
source("./getDATA.R")
source("./geneExpression_heatMap.R")


#Define the server the logic
shinyServer(function(input,output){
  
   #return pathway selected by user
   output$currentPathway <- renderText({
     paste(input$user_selected_pathway,"pathway",sep=' ')
   })
#    
   #reactive function to get the current selected pathway name
   current_pathway <- reactive({
     input$user_selected_pathway
   })
   
   #return the gene list in the pathways
   output$genes_in_currentPathway <- renderText({
     genes_in_currentPathway <- eval(parse(text=KEGG_available_pathways[[input$user_selected_pathway]])) 
   })
   
   
   #get list of genes in current pathway or user entered list
   selected_genes <- reactive({
     if(input$user_selected_pathway == "custom pathway"){
       print(user_submitted_gene_set())
     }
     else{
       eval(parse(text=KEGG_available_pathways[[input$user_selected_pathway]]))  
     }
   })
   
   #subset the gene norm counts to select only the genes found in the selected pathway or user list
   selected_geneNormCounts <- reactive({
     #filter based on user list if any
     subset(geneNormCounts, symbol %in% selected_genes())
   })
   
   
   user_submitted_gene_set <- reactive({
     temp_geneL <- unlist(strsplit(input$custom_gene_list,split=c('[\\s+,\\n+\\r+)]'),perl=T))
     #conevert everything to upper case
     temp_geneL <- toupper(temp_geneL)
   })
   
   
   
#    get the list of user submitted genes
   output$user_submitted_gene_set <- renderText({
     temp_geneL <- unlist(strsplit(input$custom_gene_list,split=c('[\\s+,\\n+\\r+)]'),perl=T))
     #conevert everything to upper case
     temp_geneL <- toupper(temp_geneL)
      print(length(temp_geneL))
      print(temp_geneL)
      temp_geneL
   })
   
   
   #return the heatMap plot
   output$heatMap <- renderPlot({
     
     #eliminate the first 3 cols to get rid of the annotation
     #getting the gene norm counts ffrom a shiny based reactive fucntion
     m <- as.matrix(selected_geneNormCounts()[4:ncol(selected_geneNormCounts())])
     
     #plot the heatmap
     get_geneExpression_heatMap(m)  #available in geneExpression_heatMap.R
   })
  
   
   
  #create summary table
  output$summary <- renderTable({     
     summary <- data.frame('Category' =  c('#genes_currentPathway', '#genes_found_in_dataset'),
                           'Value'    =  c( length(selected_genes()), nrow(selected_geneNormCounts()))
                          )
   })
   
   #prepare data for download
   output$downloadData <- downloadHandler(
     filename = function() { paste('PCBC_geneExpr_data_for_',gsub(' ','_',current_pathway()),'_pathway.csv')},
     content  = function(file){
       write.csv(selected_geneNormCounts(),file)
     })
   
})



