library(tasic2016data)
# library(scrattch.hicat) comment out to use source package in this r project
# rather than downloaded package
library(devtools)
library(roxygen2)
library(knitr)
library(testthat)
library(usethis)
library(dendextend)
library(dplyr)
library(matrixStats)
library(Matrix)
library(survival)
library(cluster)
library(impute)
library(mgcv)
library(foreign)
library(WGCNA)

######################
# load package in dev mode
######################
load_all()

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
de.param <- de_param(padj.th     = 0.5, 
                     lfc.th      = 1, 
                     low.th      = 1, 
                     q1.th       = 0.5,
                     q2.th       = NULL,
                     q.diff.th   = 0.7, 
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
                                dim.method = "WGCNA", 
                                de.param = strict.param, 
                                rm.eigen = rm.eigen)
display.result <- display_cl(onestep.result$cl, norm.dat, plot = TRUE, de.param = de.param)

# Iterative clustering
iter.result <- iter_clust(norm.dat, 
                          select.cells = select.cells, 
                          dim.method = "WGCNA", 
                          de.param = de.param, 
                          rm.eigen = rm.eigen,
                          result = onestep.result)
display.result <- display_cl(iter.result$cl, norm.dat, plot = TRUE, de.param = de.param)

####################
# Merging
####################
rd.dat <- t(norm.dat[iter.result$markers, select.cells])

merge.param <- de_param(de.score.th = 70) # The original value was 40.

merge.result <- merge_cl(norm.dat, 
                         cl = iter.result$cl, 
                         rd.dat = rd.dat,
                         de.param = merge.param)

display.result <- display_cl(merge.result$cl, 
                             norm.dat, 
                             plot = TRUE, 
                             de.param = merge.param)

# Compare pre- and post-merged clusters
# Set up the cl and cl.df objects for use with compare_annotate()
iter.cl <- setNames(as.factor(iter.result$cl), select.cells)
iter.cl.df <- data.frame(cluster_id = unique(iter.cl),
                         cluster_label = paste0("Pre-merge_cl_",unique(iter.cl)),
                         cluster_color = rainbow(length(unique(iter.cl))))
rownames(iter.cl.df) <- iter.cl.df$cluster_id

compare.result <- compare_annotate(merge.result$cl, iter.cl, iter.cl.df)
compare.result$g

# Generate comparison
compare.result <- compare_annotate(iter.result$cl, ref.cl, ref.cl.df)
# Output the plot
compare.result$g

