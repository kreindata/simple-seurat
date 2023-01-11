FeaturePlotAnalysis <- function(data, cluster = 0, numfeatures = 3) {
  cluster.markers <- FindMarkers(data, ident.1 = cluster, min.pct = 0.25)
  head <- row.names(head(cluster.markers, n = numfeatures))

  FeaturePlot(data, features = head)
}
