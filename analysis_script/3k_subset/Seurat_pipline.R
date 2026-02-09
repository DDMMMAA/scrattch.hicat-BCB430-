library(Seurat)
library(dplyr)
library(bigstatsr)
library(Matrix)
library(RcppParallel)
library(matrixStats)
library(ggplot2)

##############
# load data
# seurat object
`3k_subset_seurat` <- readRDS("D:/BCB430/scrattch.hicat-BCB430-/data/Allen/3k_subset_seurat.rds")
# expression matrix
`3k_subset_expression_matrix` <- readRDS("D:/BCB430/scrattch.hicat-BCB430-/data/Allen/3k_subset_expression_matrix.rds")

#############
# Clustering
# hicat pipline
de.param = de_param(q1.th=0.4, 
                    q.diff.th = 0.7, 
                    de.score.th=300, 
                    min.cells=50)

select.cells = colnames(`3k_subset_expression_matrix`)
set.seed(123)
buf = iter_clust(norm.dat = `3k_subset_expression_matrix`, 
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

subset.dat <- `3k_subset_expression_matrix`[, sampled.cells]

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
table(Cells(`3k_subset_seurat`) %in% names(hicat_clusters))
`3k_subset_seurat` <- AddMetaData(object = `3k_subset_seurat`, 
                                  metadata = hicat_clusters, 
                                  col.name = "hicat_cluster")

# Set the identity to be hicat cluster result
Idents(`3k_subset_seurat`) <- "hicat_cluster"

# Calculate Umap
`3k_subset_seurat` <- FindVariableFeatures(`3k_subset_seurat`) %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30)

# Visualize Umap
DimPlot(`3k_subset_seurat`, group.by = "hicat_cluster", label = TRUE, repel = TRUE) +
  ggtitle("scrattch.hicat Clusters on Seurat UMAP")

DimPlot(`3k_subset_seurat`, group.by = "class_label", label = TRUE, repel = TRUE) +
  ggtitle("class_label on Seurat UMAP")

DimPlot(`3k_subset_seurat`, group.by = "subclass_label", label = TRUE, repel = TRUE) +
  ggtitle("subclass_label on Seurat UMAP")

DimPlot(`3k_subset_seurat`, group.by = "supertype_label", label = TRUE, repel = TRUE) +
  ggtitle("supertype_labe on Seurat UMAP")

DimPlot(`3k_subset_seurat`, group.by = "cluster_label", label = TRUE, repel = TRUE) +
  ggtitle("cluster on Seurat UMAP")

DimPlot(`3k_subset_seurat`, group.by = "roi", label = TRUE, repel = TRUE) +
  ggtitle("roi on Seurat UMAP")

DimPlot(`3k_subset_seurat`, group.by = "sex", label = TRUE, repel = TRUE) +
  ggtitle("sex on Seurat UMAP")

DimPlot(`3k_subset_seurat`, group.by = "age_cat", label = TRUE, repel = TRUE) +
  ggtitle("age_cat on Seurat UMAP")
