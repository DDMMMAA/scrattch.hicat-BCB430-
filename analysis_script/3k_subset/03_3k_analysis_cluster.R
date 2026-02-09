## Description/Goal:
## Perform final round of clustering
library(dplyr)
# library(scrattch.hicat)
# library(scrattch.bigcat)
library(bigstatsr)
library(Matrix)
library(RcppParallel)
library(matrixStats)


#################################################
## Load data
#################################################
# Expression matrix
load("D:/BCB430/scrattch.hicat-BCB430-/data/workspace/norm.dat_cl100_ss200Kcells.rda")

big.dat <- norm.dat

# Metadata
load("D:/BCB430/scrattch.hicat-BCB430-/data/workspace/samp.dat.filtered_prelim.20220705.rda")

anno.df.clean <- samp.dat.filtered
anno.df.clean$sample_id <- rownames(anno.df.clean)

select.cells = anno.df.clean$sample_id


#################################################
## Second clustering round
#################################################
de.param = de_param(q1.th=0.4, 
                    q.diff.th = 0.7, 
                    de.score.th=300, 
                    min.cells=50)

set.seed(123)
resulti = iter_clust(norm.dat = big.dat, 
                     select.cells = select.cells, 
                     de.param =de.param, 
                     prefix="./data/prelim_cluster/cluster_0916/cluster_intermed_files", 
                     max.cl.size=200, 
                     split.size = 50, 
                     verbose=1, 
                     sampleSize=10000)

save(resulti, file=paste0("./data/workspace/cluster_0916/resulti_cluster.rda"))


## Merge clusters
load("./data/workspace/cluster_0916/resulti_cluster.rda")
cl = resulti$cl
markers = resulti$markers

sampled.cells = sample_cells(cl, 100)
# sampled.cells = sample(sampled.cells, 200000, replace = F)
norm.dat = big.dat
select.markers = intersect(markers, rownames(big.dat))

de.param = de_param(q1.th=0.4, q.diff.th = 0.7, de.score.th=300, min.cells=100)
merge.result = merge_cl(norm.dat=norm.dat, cl=cl, rd.dat.t=norm.dat[select.markers,], de.param=de.param, verbose=TRUE)



#################################################
## Save objects for future use
#################################################
save(merge.result, file= "./data/workspace/cluster_0916/merge.result_th300_mincells100.rda") ## clustering results
save(norm.dat, file = "./data/workspace/cluster_0916/norm.dat_cluster0916_cl100.rda") ## subsampled cell-by-gene matrix