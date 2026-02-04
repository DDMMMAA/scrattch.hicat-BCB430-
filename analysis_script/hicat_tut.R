library(tasic2016data)
# library(scrattch.hicat) comment out to use source package in this r project
# rather than downloaded package
library(dendextend)
library(dplyr)
library(matrixStats)
library(Matrix)

######################
# data preparation
######################
# Subset small set of interneuron for the sake of performance
# extract the sample name of subset cell
select.cells <- tasic_2016_anno %>%
  filter(primary_type_label != "unclassified") %>%
  filter(grepl("Igtp|Ndnf|Vip|Sncg|Smad3",primary_type_label)) %>%
  select(sample_name) %>%
  unlist()
# Use subset cell name to subset cells
ref_anno <- tasic_2016_anno %>%
  filter(sample_name %in% select.cells)

# construct df of unique primary_type id, type, color, and broad type
ref.cl.df <- ref_anno %>%
  select(primary_type_id, 
         primary_type_label, 
         primary_type_color, 
         broad_type) %>%
  unique()
# assign colname
colnames(ref.cl.df)[1:3] <- c("cluster_id", "cluster_label", "cluster_color")

# Sort by cluster_id
ref.cl.df <- arrange(ref.cl.df, cluster_id)
row.names(ref.cl.df) <- ref.cl.df$cluster_id

# Set up the ref.cl factor object
ref.cl <- setNames(factor(ref_anno$primary_type_id), ref_anno$sample_id)

######################
# Data normalization
######################
# convert raw counts to counts per million read
tasic_2016_cpm <- cpm(tasic_2016_counts[,select.cells])
# "+ 1" is a pseudo counts to avoid log2(0) and "0 + 1" is transformed into 0
norm.dat <- log2(tasic_2016_cpm + 1)

######################
# parameter setting
######################
de.param <- de_param(padj.th     = 0.1, 
                     lfc.th      = 1, 
                     low.th      = 1, 
                     q1.th       = 0.1,
                     q2.th       = NULL,
                     q.diff.th   = 0, 
                     de.score.th = 40,
                     min.cells = 10)

######################
# Dimension Filtering
######################
# Collect number of feature with expression value > 0
# Not quit sure whats is doing here? Is it for the sake of PCA?
# Why log-transform the number of gene with non-zero expression?
gene.counts <- colSums(norm.dat > 0)
rm.eigen <- matrix(log2(gene.counts), ncol = 1)
row.names(rm.eigen) <- names(gene.counts)
colnames(rm.eigen) <- "log2GeneCounts"

######################
# Clustering
######################
# One step clustering using WGCNA as dimension reduction
strict.param <- de_param(de.score.th = 500)

onestep.result <- onestep_clust(norm.dat, 
                                select.cells = select.cells, 
                                dim.method = "pca", 
                                de.param = strict.param, 
                                rm.eigen = rm.eigen)
display.result <- display_cl(onestep.result$cl, norm.dat, plot = TRUE, de.param = de.param)
