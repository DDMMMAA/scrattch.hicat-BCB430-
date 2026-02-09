library(Seurat)
library(dplyr)
library(bigstatsr)
library(Matrix)
library(RcppParallel)
library(matrixStats)
library(ggplot2)
library(devtools)
library(roxygen2)
library(knitr)

##############
# load current r.proj as r package
load_all()

##############
# load data
# seurat object
glia_only_seurat <- readRDS("D:/BCB430/scrattch.hicat/data/Allen/glia_only_SeuratObj.rds")
# expression matrix
glia_only_counts_matrix <- LayerData(glia_only_seurat, assay = "RNA", layer = "counts")

#############
# Clustering
# hicat pipline
de.param = de_param(q1.th=0.4, 
                    q.diff.th = 0.7, 
                    de.score.th=300, 
                    min.cells=50)

select.cells = colnames(glia_only_counts_matrix)
set.seed(123)
buf = iter_clust(norm.dat = glia_only_counts_matrix, 
                 select.cells = select.cells, 
                 de.param =de.param, 
                 prefix="./data/prelim_cluster/cluster_0220705", 
                 max.cl.size=200, 
                 split.size = 50, 
                 verbose=1, 
                 sampleSize=10000)

# Sample cells for merging
sampled.cells = sample_cells(buf$cl, 100)
# Uncomment the line below on whole dataset
# sampled.cells = sample(sampled.cells, 200000, replace = F)

subset.dat <- glia_only_counts_matrix[, sampled.cells]

# Note the original paper changed min.cells=100
de.param = de_param(q1.th=0.4, 
                    q.diff.th = 0.7, 
                    de.score.th=300, 
                    min.cells=50)

merge.result = merge_cl(norm.dat=subset.dat, 
                        cl=buf$cl, 
                        rd.dat.t=subset.dat[buf$markers,], 
                        de.param=de.param, 
                        verbose=TRUE)

##################
# Insert hicat cluster result into seurat
hicat_clusters <- merge.result$cl
table(Cells(glia_only_seurat) %in% names(hicat_clusters))
glia_only_seurat <- AddMetaData(object = glia_only_seurat, 
                                  metadata = hicat_clusters, 
                                  col.name = "hicat_cluster")

# Set the identity to be hicat cluster result
Idents(glia_only_seurat) <- "hicat_cluster"

# Calculate Umap
glia_only_seurat <- FindVariableFeatures(glia_only_seurat) %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30, reduction.name = "umap_dim30")

# Visualize Umap
# 1. Define the list of variables to plot
group_vars <- c("hicat_cluster", "class_label", "subclass_label", 
                "supertype_label", "cluster_label", "roi", "sex", "age_cat")

# 2. Loop through each variable
for (var in group_vars) {
  
  # Print progress to console
  message(paste("Plotting and saving:", var))
  
  # Generate the plot
  p <- DimPlot(glia_only_seurat, group.by = var, label = TRUE, repel = TRUE, reduction = "umap_dim30") +
    ggtitle(paste(var, "on Seurat UMAP"))  + NoLegend()
  
  # 3. Save the plot
  # Construct filename and path
  filename <- paste0("data/Result/cluster/non-neuro/", "UMAP_dim30_", var, ".png")
  
  ggsave(filename = filename, plot = p, width = 10, height = 7.5, dpi = 320)
}