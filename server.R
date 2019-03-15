#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library("Seurat")
library(shiny)

#Load the differentially epxressed Genes (APOE and 100 other random picked genes)
load( "toyDE.RData")
differentialExpre = toyDifferentialExpre
rm( toyDifferentialExpre)
names( differentialExpre) = c('Neuron', 'Oligodendrocytes', 'Astrocytes', 'Microglia', 'OPC', 'Endothelial')

#Function to plot the differential expression
PlotDifferentialExpre = function ( targetGene) {
  ret = data.frame( )
  for ( cmp in names(differentialExpre ) ) {
    c = differentialExpre[[cmp]][targetGene,]
    if (!is.na( c$p_val)) {
      foldChange = sprintf( "%2.2f", c$avg_logFC)
      p.value = sprintf( "%2.2e", c$p_val)
      p.adj = sprintf( "%2.2e", c$p_val_adj)
      ret[cmp,"Cell"]=cmp
      ret[cmp,"FoldChange"] = foldChange
      ret[cmp,"fc"] = c$avg_logFC
      ret[cmp,"P value"] = p.value
      ret[cmp,"P adj"] = p.adj
      ret[cmp,"CellDetails"]= sprintf( "%s\nlogFC=%2.2f\np=%2.0e\nadjc=%2.0e", cmp, c$avg_logFC, c$p_val, c$p_val_adj)
    }
  }
  

  ret.gr = ggplot(ret, aes(x=CellDetails, y=fc)) +
    geom_segment( aes(x=CellDetails, xend=CellDetails, y=0, yend=fc), color="grey") +
    geom_point( color="orange", size=3) +
    theme_light(base_size=15) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    xlab("") +
    ylab("Log fold change")
  
  return( ret.gr)
}





load( 'toyCells.RData')
label_data_marker_gene = toy
rm( toy)
label_data_marker_gene <- SetAllIdent(object = label_data_marker_gene, id = "res.0.6")

#JDA
tsnePlot = TSNEPlot(object = label_data_marker_gene, pt.size=0.5, label.size = 2.8, do.label = TRUE, no.legend=TRUE)


shinyServer(function(input, output) {
  
  
  output$clusters <- renderPlot({
    if (input$Cluster == '') {
      tsnePlot
    } else {
      
      cellsHighlight=which( label_data_marker_gene@ident==input$Cluster )
      TSNEPlot(object = label_data_marker_gene, do.label = F, cells.highlight=cellsHighlight, cols.highlight = "orange", plot.order=input$Cluster) 
  }
    
    
  })
  
  output$gene <- renderPlot({
    gene= input$Gene
    FeaturePlot(object = label_data_marker_gene, features.plot = gene , cols.use = c("lightgrey", "blue"), reduction.use = "tsne", min.cutoff = 0)
  })
  
  output$cellTypeDE = renderPlot({
    PlotDifferentialExpre(input$Gene)
  })
  
})
