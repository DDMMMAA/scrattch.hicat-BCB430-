# This file is for load the data
library(anndataR)
library(Seurat)
library(Matrix)

#### Data Directory constants ###########
raw_metadata_dir <- "D:/BCB430/Data/raw_metadata.csv"
metadata_dir <- "D:/BCB430/Data/metadata.csv"
cell_names_dir <- "D:/BCB430/Data/cell_names.txt"
gene_names_dir <- "D:/BCB430/Data/gene_names.txt"
cell_indices_dir <- "D:/BCB430/Data/cell_indices.npy"
hpf_expression_h5_dir <- "D:/BCB430/Data/HPF.h5ad"

#### Load data #############################
# Load Seurat Object
SeuratObj <- read_h5ad(hpf_expression_h5_dir, as = "Seurat")

# Ensure default assay is RNA
DefaultAssay(SeuratObj) <- "RNA"

# Assign layer X to both layer data and layer counts, since layer x is already normalized
LayerData(SeuratObj, assay = "RNA", layer = "data") <- LayerData(SeuratObj, assay = "RNA", layer = "X")
LayerData(SeuratObj, assay = "RNA", layer = "counts") <- LayerData(SeuratObj, assay = "RNA", layer = "X")

# Subset Non-neuronal cell
Non_neuro_seuratobj <- subset(SeuratObj, class_label == "OPC-Oligo" |
                              class_label == "Immune" | 
                              class_label == "Astro-Epen" |
                              class_label == "Vascular")
# Subset Neuronal cell
neuro_seuratobj <- subset(SeuratObj, class_label != "OPC-Oligo" &
                           class_label != "Immune" & 
                           class_label != "Astro-Epen" &
                           class_label != "Vascular")

# Subset cells with doublet score == 0
zero_doublet_seuratobj <- subset(SeuratObj, doublet_score == 0)

# Subset male cells
male_seuratobj <- subset(SeuratObj, sex == "M")

# Subset female cells
female_seuratobj <- subset(SeuratObj, sex == "F")

# Subset adult cells
adult_seuratobj <- subset(SeuratObj, age_cat == "adult")

# subset aged cell
aged_seuratobj <- subset(SeuratObj, age_cat == "aged")

# Subset cell of different roi
ENT_seuratobj <- subset(SeuratObj, roi == "ENT")

HIP_seuratobj <- subset(SeuratObj, roi == "HIP")

HIP_CA_seuratobj <- subset(SeuratObj, roi == "HIP - CA")

PAR_POST_PRE_SUB_ProS_seuratobj <- subset(SeuratObj, roi == "PAR-POST-PRE-SUB-ProS")

#### Subset 3k cell ###############
# random subset 3k cell
set.seed(123) 
random_cells <- sample(colnames(SeuratObj), size = 3000)
subset_obj <- subset(SeuratObj, cells = random_cells)

# Extract the expression matrix properly
expression_matrix <- LayerData(subset_obj, assay = "RNA", layer = "data")

#### Save obj to disk #################
# save whole seurat obj
saveRDS(SeuratObj, file = "analysis_script/mouse_aging_scRNAseq-main/my_scripts/SeuratObj.rds")

saveRDS(Non_neuro_seuratobj, file = "analysis_script/mouse_aging_scRNAseq-main/my_scripts/non_neuro_seuratobj.rds")

saveRDS(neuro_seuratobj, file = "analysis_script/mouse_aging_scRNAseq-main/my_scripts/neuro_seuratobj.rds")

saveRDS(zero_doublet_seuratobj, file = "analysis_script/mouse_aging_scRNAseq-main/my_scripts/zero_doublet_seuratobj.rds")

saveRDS(male_seuratobj, file = "analysis_script/mouse_aging_scRNAseq-main/my_scripts/male_seuratobj.rds")

saveRDS(female_seuratobj, file = "analysis_script/mouse_aging_scRNAseq-main/my_scripts/female_seuratobj.rds")

saveRDS(adult_seuratobj, file = "analysis_script/mouse_aging_scRNAseq-main/my_scripts/adult_seuratobj.rds")

saveRDS(aged_seuratobj, file = "analysis_script/mouse_aging_scRNAseq-main/my_scripts/aged_seuratobj.rds")

saveRDS(ENT_seuratobj, file = "analysis_script/mouse_aging_scRNAseq-main/my_scripts/ENT_seuratobj.rds")

saveRDS(HIP_seuratobj, file = "analysis_script/mouse_aging_scRNAseq-main/my_scripts/HIP_seuratobj.rds")

saveRDS(HIP_CA_seuratobj, file = "analysis_script/mouse_aging_scRNAseq-main/my_scripts/HIP_CA_seuratobj.rds")

saveRDS(PAR_POST_PRE_SUB_ProS_seuratobj, file = "analysis_script/mouse_aging_scRNAseq-main/my_scripts/PAR_POST_PRE_SUB_ProS_seuratobj.rds")

