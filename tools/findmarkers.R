library("Seurat")  

#load the Seurat Object
load( "toyCells.RData" )
toy <- SetAllIdent(object = toy, id = "res.0.6")


head( toy@meta.data )
cell.tags  = names(toy@ident)
clusters = as.numeric(toy@ident )  #from 1 to 14 
              #N N N N N O N N N A N M P E
cell.type = c( 0,0,0,0,0,1,0,0,0,2,0,3,4,5 ) + 1
cell.names = data.frame( i=0:5, name=c('Neuron', 'Oligodendrocytes', 'Astrocytes', 'Microglia', 'OPC', 'Endothelial'), stringsAsFactors = F)
    #Neuron 0
    #Oligo  1
    #Astro  2
    #Micro  3
    #OPC    4
    #Endo   5
which( cell.type==1)-1 #OK for neurons
which( cell.type==2)-1 #OK for Oligo
which( cell.type==3)-1 #OK for Astro
which( cell.type==4)-1 #OK for Micro
which( cell.type==5)-1 #OK for OPC
which( cell.type==6)-1 #OK for Endo
length( cell.type)

clusters.grouped = sapply( clusters, function(x) cell.type[x])
head( clusters.grouped)
head( clusters)

#verify
all( clusters.grouped[which( clusters == 1)]==1)
all( clusters.grouped[which( clusters == 2)]==1)
all( clusters.grouped[which( clusters == 3)]==1)
all( clusters.grouped[which( clusters == 4)]==1)
all( clusters.grouped[which( clusters == 5)]==1)
all( clusters.grouped[which( clusters == 7)]==1)
all( clusters.grouped[which( clusters == 8)]==1)
all( clusters.grouped[which( clusters == 9)]==1)
all( clusters.grouped[which( clusters == 11)]==1)
all( clusters.grouped[which( clusters == 6)]==2)#Oligo
all( clusters.grouped[which( clusters == 10)]==3)#Astro
all( clusters.grouped[which( clusters == 12)]==4)#Migco
all( clusters.grouped[which( clusters == 13)]==5)#OPC
all( clusters.grouped[which( clusters == 14)]==6)#Endo

clusters.grouped.df2 = as.data.frame(clusters.grouped-1)
names( clusters.grouped.df2) = "celltype.grouped"
row.names( clusters.grouped.df2) = cell.tags
head( clusters.grouped.df2)
table( clusters.grouped.df2)
toy =  AddMetaData(object = toy, metadata = clusters.grouped.df2, col.name="celltype.grouped")
head( toy@meta.data)
toy = SetAllIdent(toy, id='celltype.grouped')
table( toy@ident)
differentialExpre = list()

#To get as many genes as possibly tested
  #logfc.threshold=0.01   
#For testing (this should run faster):
  logfc.threshold=0.25
for ( i in 0:5) {
  ct = cell.names[cell.names$i==i,]$name
  cat (sprintf( "Processing %i: %s\n", i, ct))
  a = FindMarkers(object = toy, ident.1 = i, logfc.threshold=logfc.threshold)
  differentialExpre[[ct]] =  a
}

#Just save the list with all of the differential expression 
#save( differentialExpre, file="diffExpre/differentialExpre.RData")

#query the differential expression
target.gene = "APOE"
lapply( differentialExpre, function(X) X[target.gene,])

