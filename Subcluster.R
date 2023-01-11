Subcluster <- function(data, clusternums, numdims = PCA) {
  supercluster <- subset(data, idents = clusternums)
  DimPlot(supercluster, reduction = "umap")
  supercluster <- RunPCA(supercluster)
  supercluster <- FindNeighbors(supercluster, dims = 1:PCA)
  supercluster <- FindClusters(supercluster, resolution = 0.5)
  supercluster <- RunUMAP(supercluster, dims = 1:PCA)
DimPlot(supercluster, reduction = "umap")
}
