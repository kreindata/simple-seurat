Cluster <- function(data, dims, strength) {
  # test for input
  if(missing(data)) stop('Missing input')
  
  # test for dgCMatrx and convert to Seurat if applicable
  if(substr(type_sum(data), 1, 8) == "dgCMatrx") {
    data <- CreateSeuratObject(counts = data, project = "clusters", min.cells = 3, min.features = 200)
  }
  
  # test for Seurat object after any matrix conversion, stop if failed
  if(!(substr(type_sum(data), 1, 6) == "Seurat")) stop('Argument should point to data matrix or prepared Seurat object')

  # skim outliers and normalize
  subdata <- subset(sdata, subset = nCount_RNA > quantile(sdata$nFeature_RNA, probs = .02, names = FALSE)  & nFeature_RNA < quantile(sdata$nFeature_RNA, probs = .98, names = FALSE))
  subdata <- NormalizeData(subdata)

  # examine most variable features available in the data and re-normalize, then run PCA analysis
  subvars <- FindVariableFeatures(subdata, selection.method ='vst', nfeatures = quantile(subdata$nFeature_RNA, probs = .95))
  genes <- rownames(subvars)
  subvars <- ScaleData(subvars, features = genes)
  subvars <- RunPCA(subvars, features = VariableFeatures(object = subvars))

  # sort and identify components with strongest percent variance explained scores
  plot <- ElbowPlot(subvars)

  # if strength not set by user, defaults to 0, causing scaledPCA to take all above-average components by percent variance explained
  if(missing(strength)) strength <- 0
  scaledPCA <- scale(plot$data$stdev)
  scaledPCA <- subset(scaledPCA, scaledPCA[,1] > strength)
  
  # allows user to set number of dimensions, or by default takes number of dimensions that best fit
  if(missing(dims)) dims <- nrow(scaledPCA)

  # clustering analysis and visual plotting
  clusters <- FindNeighbors(subvars, dims = 1:dims)
  clusters <- FindClusters(clusters, resolution = 0.5) # edit to add names during this portion?
  clusters <- RunUMAP(clusters, dims = 1:dims)
  DimPlot(clusters, reduction = "umap", label = TRUE)
}
