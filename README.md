# simple-seurat
Simplifying Seurat data processing, clustering, and analysis. This setup assumes that your data has already been cleaned. If not, you can refer to the [Seurat Guided Clustering Tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) for advice on data cleaning. As in the Seurat tutorial, this function utilizes the Seurat, dplyr, and patchwork libraries.

The Cluster function collects and runs the Seurat functions necessary to identify an ideal number of clusters. This is done using PCA analysis, normalizing the Percent Variance Explained around zero and then taking all positive amounts. In my experience, this identifies a strong 'elbow' in the PCA plot, but may be underinclusive, depending on your goal.

The Subcluster function is meant for Seurat objects which have already been clustered, either by the cluster function above or manually. It takes as its argument the numerical cluster identifier(s). Subclustering can be run on a single cluster or on a vector of clusters. The dimensional analysis can be massaged to keep the subcluster a similar shape to the original superclusters. 

The FeaturePlotAnalysis function is a clustering quality analysis tool built on existing Seurat functions which takes as arguments the data object, which cluster to examine, and the number of features to examine. The output is a FeaturePlot of those features against the overall clustering plot, which allows researchers to quickly and easily identify the top markers associated with each cluster, as well as examine whether those clusters are appropriately associated with those markers. It can also help identify whether some clusters have been grouped around markers which are common across all clusters. This seems especially common with the initial Seurat clusters chosen by these functions (cluster 0).
