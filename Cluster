Cluster <- function(data) {
  sdata <- CreateSeuratObject(counts = data, project = "clusters", min.cells = 3, min.features = 200)

  subdata <- subset(sdata, subset = nCount_RNA > quantile(sdata$nFeature_RNA, probs = .02, names = FALSE)  & nFeature_RNA < quantile(sdata$nFeature_RNA, probs = .98, names = FALSE))
  subdata <- NormalizeData(subdata)

  subvars <- FindVariableFeatures(subdata, selection.method ='vst', nfeatures = quantile(subdata$nFeature_RNA, probs = .95))
  genes <- rownames(subvars)
  subvars <- ScaleData(subvars, features = genes)
  subvars <- RunPCA(subvars, features = VariableFeatures(object = subvars))

  plot <- ElbowPlot(subvars)

  scaledPCA <- scale(plot$data$stdev)
  scaledPCA <- subset(scaledPCA, scaledPCA[,1] > 0)
  PCA <- nrow(scaledPCA)

  clusters <- FindNeighbors(subvars, dims = 1:PCA)
  clusters <- FindClusters(clusters, resolution = 0.5)
  clusters <- RunUMAP(clusters, dims = 1:PCA)
  DimPlot(clusters, reduction = "umap")
}
