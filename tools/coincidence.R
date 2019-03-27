
#install.packages('svs')
library( svs )
#install.packages('reshape')
library( reshape)




CoincidenceSeurat = function( cluster.A, cluster.B) {
  
  #work with the cells that passed QC in both strategies
  inter<-intersect(names( cluster.A), names( cluster.B))
  
  #Get the cells in the intersection in a common order
  #First cluster.A
  cluster.A.ndx = match(inter, names(cluster.A))
  #check
  all(inter == names(cluster.A[cluster.A.ndx]))
  cluster.A.matched<-cluster.A[cluster.A.ndx]
  #Then cluster.B
  cluster.B.ndx<-match(inter, names(cluster.B))
  all(inter == names(cluster.B[cluster.B.ndx]))
  cluster.B.matched<-cluster.B[cluster.B.ndx]
  
  #check
  all(names(cluster.A.matched) == names(cluster.B.matched))
  
  
  cluster.A.matched.int<-as.integer(as.character(cluster.A.matched))
  cluster.B.matched.int<-as.integer(as.character(cluster.B.matched))
  conc<-matrix(0L, nrow = max(cluster.A.matched.int)+1, ncol = max(cluster.B.matched.int)+1)
  
  for( i in 1:length( cluster.A.matched.int)) {
    conc[cluster.A.matched.int[i]+1, cluster.B.matched.int[i]+1] =conc[cluster.A.matched.int[i]+1, cluster.B.matched.int[i]+1]+1
  }
  
  
  #Get the normalized pointwise mutual information
  npmi = pmi(conc, normalize = T) #Get the normalized pointwize mutual information
  
  #calculate jaccard
  jaccard<-matrix(0L, nrow = max( cluster.A.matched.int)+1, ncol = max( cluster.B.matched.int )+1)
  cluster.A.matched.int.clustersN<-summary( cluster.A)
  cluster.B.matched.int.clustersN<-summary( cluster.B)
  for( i in 1:nrow(conc)) {
    for( j in 1:ncol(conc)) {
      jaccard[i, j] = conc[i,j] / (cluster.A.matched.int.clustersN[i] + cluster.B.matched.int.clustersN[j] - conc[i,j])
    }
  }
  
  #format the data for ggplot
  data.toPlot<-melt(npmi)
  data.toPlot<-merge(data.toPlot, melt( jaccard), by=c("X1", "X2"))
  data.toPlot<-merge(data.toPlot, melt( conc), by=c("X1", "X2"))
  names(data.toPlot)<-c("clusterA", "clusterB", "npmi", "jaccard", "intersection")
  data.toPlot$textSize<-3*(data.toPlot$intersection>=10)
  clusterA.label<-paste(names(summary(cluster.A)), " N=", summary(cluster.A))
  clusterB.label<-paste(names(summary(cluster.B)), " N=", summary(cluster.B))
  
  p <- ggplot(data = data.toPlot, aes(x = clusterA, y = clusterB, size=jaccard, color=npmi, label=intersection)) +
    geom_point() +scale_color_gradient(low = "blue", high = "red") +
    geom_text(size=data.toPlot$textSize, color="black", fontface ="bold", angle=45) +
    scale_y_continuous(labels= clusterB.label, breaks = 1:length(clusterB.label)) +
    scale_x_continuous(labels= clusterA.label, breaks=1:length(clusterA.label))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8), axis.text.y = element_text(size = 8))
  p
}


#compare the clustering of the toy dataset using two distinct resolutions
toy.A <- SetAllIdent(object = toy, id = "res.0.4")
toy.A.ident<-toy.A@ident

toy.B <- SetAllIdent(object = toy, id = "res.0.6")
toy.B.ident = toy.B@ident

#cluster 0 (from res.0.4) is splitted into two cluster (0 and 1) in t3he res.0.6
#Note that some clusters might have distinct ID between distinct resolutions, but they are grouping the same nuclei
#For instance all nuclei of cluster 5 in res.0.4 are reasigned to cluster 7 (res.0.6)
CoincidenceSeurat( toy.A.ident, toy.B.ident)



