#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)


load("toyGeneNames.RData")
geneNames = toyGeneNames
rm( toyGeneNames)

#JDA
optionForClusters=c(Choose='', Exitatory_1=0, Exitatory_2=1, Exitatory_4=3,
                    Exitatory_5=2, Exitatory_6=8, Exitatory_7=4, Exitatory_8=10,
                    Inhibitory_1=6, Inhibitory_6=7, 
                    Oligodendrocyte=5, Astrocyte=9, 
                    Microglia=11, OPC=12, Endothelial=13)

#Jason
#optionForClusters = c(Choose='',unique( label_data_marker_gene@active.ident))

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  title = "Single-nuclei RNA-seq expression browser",
  titlePanel("Single-nuclei RNA-seq expression browser"),
  
  # Sidebar with a slider input for number of bins 
  fluidRow(
    column(3,
          # textInput("Gene", "Gene:", value = "APOE", width = NULL, placeholder = "Provide the gene Official Symbol" )
          selectInput('Gene', 'Gene:', geneNames, multiple=FALSE, selectize=TRUE),
          selectInput('Cluster', 'Highlight', optionForClusters, selectize=FALSE)
    ) ),
    br(),
    # Show a plot of the generated distribution
    plotOutput("clusters", width = 800, height = 600),
    plotOutput("gene", width = 800, height = 600),
    br(),
    plotOutput("cellTypeDE", width = 800, height = 300),
    br(),
    p("Add a text with the citation info")
    )
  )

