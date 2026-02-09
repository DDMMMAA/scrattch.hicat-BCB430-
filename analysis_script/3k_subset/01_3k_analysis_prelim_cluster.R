library(Seurat)
#################
# Load expression matrix
`3k_subset_expression_matrix` <- readRDS("D:/BCB430/scrattch.hicat-BCB430-/data/Allen/3k_subset_expression_matrix.rds")

################
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

merge.result = merge_cl(norm.dat=`3k_subset_expression_matrix`, 
                        cl=buf$cl, 
                        rd.dat.t=`3k_subset_expression_matrix`[buf$markers,], 
                        de.param=de.param, 
                        verbose=TRUE)

#################################################
## Save objects for future use
#################################################
save(merge.result, file="data/workspace/merge.result_th300_cluster_0220705.rda")
norm.dat <- `3k_subset_expression_matrix`
save(norm.dat, file = "data/workspace/norm.dat_cl100_ss200Kcells.rda")
