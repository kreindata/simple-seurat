# simple-seurat
Simplifying Seurat data processing, clustering, and analysis. This setup assumes that your data has already been cleaned. If not, you can refer to the [Seurat Guided Clustering Tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) for advice on data cleaning. As in the Seurat tutorial, this function utilizes the Seurat, dplyr, and patchwork libraries.

This function collects and runs the Seurat functions necessary to identify an ideal number of clusters. This is done using PCA analysis, normalizing the Percent Variance Explained around zero and then taking all positive amounts. In my experience, this identifies a strong 'elbow' in the PCA plot, but may be underinclusive, depending on your goal.

