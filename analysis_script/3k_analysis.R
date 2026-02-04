library(Seurat)
#################
de.param = de_param(q1.th=0.4, q.diff.th = 0.7, de.score.th=300, min.cells=50)
select.cells = colnames(`3k_subset_expression_matrix`)
set.seed(123)
buf = iter_clust(norm.dat = `3k_subset_expression_matrix`, select.cells = select.cells, de.param =de.param, prefix="./prelim_cluster/cluster_0220705",
                                 max.cl.size=200, split.size = 50, verbose=1, sampleSize=10000)