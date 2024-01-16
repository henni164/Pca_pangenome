#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(phytools)
library(ggplot2)
library(ggtree)
library(tidytree)

# Define UI for application that draws a histogram

ui <- fluidPage(
  
  titlePanel("Interactive Phylogenetic Tree Viewer"),
  
  sidebarLayout(
    sidebarPanel(
      radioButtons("shape","Tree shape",
                   choices = list("Circular" = 1, "Linear" = 2), selected = 1),
      radioButtons("view","Tree view",
                   choices = list("Full" = 1, "USA_Sexual" = 2, "Australian" = 3), selected = 1)),
    mainPanel(plotOutput("tree", width = "100%", height = "100%")
    )
  )
)

server <- function(input, output) {
  
  output$tree <- renderPlot({
    
    choices <- c("circular", "rectangular")
    
    pca203tree <- read.tree("C:/Users/HEN294/OneDrive - CSIRO/Documents/GitHub/Pca_pangenome/Figure1/RAxML_bipartitionsBranchLabels.all_isolates_Pca203_biallelic_500")
    midrooted <- midpoint_root(pca203tree)
    
    
    if(input$view == "1") {
      
      grouped <- groupClade(midrooted, c(492, 339, 452, 424, 421, 488, 511))
      
      ggtree(grouped, layout = choices[as.numeric(input$shape)], aes(color = group)) + scale_color_manual(values = c("black","#00B7BF","#FB0449","#FB0449","#FF8503","#FF8503","#D253FF","#00B7BF")) + geom_tiplab()

        } else if (input$view == "2") {
      
          grouped <- groupClade(midrooted, c(492, 339, 452, 424, 421, 488, 511))
          p <- ggtree(grouped, layout = choices[as.numeric(input$shape)], aes(color = group)) + scale_color_manual(values = c("black","#00B7BF","#FB0449","#FB0449","#FF8503","#FF8503","#D253FF","#00B7BF")) + geom_tiplab()
          
          p %>% collapse(492, mode = "none", clade_name = "USA") %>%
            collapse(339, mode = "none", clade_name = "AUS") %>%
            collapse(452, mode = "none", clade_name = "AUS") %>%
            collapse(424, mode = "none", clade_name = "SA") %>%
            collapse(421, mode = "none", clade_name = "SA") %>%
            collapse(488, mode = "none", clade_name = "TW") + geom_cladelab(node = 339, label = "AUS", geom = "text") +
            geom_cladelab(node = 492, label = "USA", geom = "text") +
            geom_cladelab(node = 452, label = "AUS", geom = "text") +
            geom_cladelab(node = 424, label = "SA", geom = "text") +
            geom_cladelab(node = 421, label = "SA", geom = "text") +
            geom_cladelab(node = 488, label = "TW", geom = "text")
      
    } else if (input$view == "3") {
      
      grouped <- groupClade(midrooted, c(492, 339, 452, 424, 421, 488, 511))
      p <- ggtree(grouped, layout = choices[as.numeric(input$shape)], aes(color = group)) + scale_color_manual(values = c("black","#00B7BF","#FB0449","#FB0449","#FF8503","#FF8503","#D253FF","#00B7BF")) + geom_tiplab()
      p %>% collapse(510, mode = "none", clade_name = "USA") %>% 
        collapse(492, mode = "none", clade_name = "USA") %>% 
        collapse(424, mode = "none", clade_name = "SA") %>%
        collapse(421, mode = "none", clade_name = "SA") + geom_cladelab(node = 510, label = "USA", geom = "text") +
        geom_cladelab(node = 492, label = "USA", geom = "text") +
        geom_cladelab(node = 424, label = "SA", geom = "text") +
        geom_cladelab(node = 421, label = "SA", geom = "text")
      
    }
    
  }, height = 2000, width = 2000)
  
}


# Run the application 
shinyApp(ui = ui, server = server)
