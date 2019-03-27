#This function calculates the entropy of the single-nuclei clusters 
#based on the distribution of donors per cluster
#it expects a Seurat object with the slotes ident and orig.ident defined
# ident: specifies the cluster assigned to each nuclei
# orig.ident: specifies the donor (sample) of the nuclei

#install.packages('entropy')
library( entropy)

EntropySeurat = function (seuratObject) {
  
  distrib = table(seuratObject@ident, seuratObject@meta.data$orig.ident)
  distrib.n = dim( distrib)[1]
  distrib.subjects = dim( distrib)[2]
  ret.entro = apply( distrib, 1, function(X) {entropy(X, unit = "log2")})
  ret.max = entropy( rep(1, distrib.subjects)/distrib.subjects, unit="log2")
  ret=list( entropy=ret.entro, max=ret.max)
  return(ret)
}


EntropySeurat( toy)
