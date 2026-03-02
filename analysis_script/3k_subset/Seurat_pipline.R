library(Seurat)
library(dplyr)
library(bigstatsr)
library(Matrix)
library(RcppParallel)
library(matrixStats)
library(ggplot2)
library(devtools)
library(knitr)
library(roxygen2)
library(testthat)
library(usethis)
library(patchwork)
library(impute)
library(WGCNA)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

##############
# load current r project as package
load_all()

##############
# load data
# seurat object
`3k_subset_seurat` <- readRDS("D:/scrattch.hicat-BCB430-/data/Allen/3k_subset_seurat.rds")
# expression matrix
`3k_subset_expression_matrix` <- readRDS("D:/scrattch.hicat-BCB430-/data/Allen/3k_subset_expression_matrix.rds")

############
# EDA
# 1. Define the variables you want to plot
vars_to_plot <- c("subclass_label", "roi", "sex", "age_cat")

# 2. Initialize a list to store the plots
pie_plots <- list()

# 3. Loop through each variable to create a pie chart
for (var in vars_to_plot) {
  
  # A. Summarize the data (Count and Percent)
  # We use .data[[var]] to dynamically select the column inside the loop
  plot_data <- `3k_subset_seurat`@meta.data %>%
    group_by(.data[[var]]) %>%
    summarise(count = n()) %>%
    mutate(
      prop = count / sum(count),
      percentage = round(prop * 100, 1),
      # Create a label that looks like: "150 (10.5%)"
      label_text = paste0(count, "\n(", percentage, "%)") 
    ) %>%
    arrange(desc(.data[[var]])) # Fix order for plotting
  
  # B. Create the Pie Chart
  # Pie charts in ggplot are just Bar Charts with polar coordinates
  p <- ggplot(plot_data, aes(x = "", y = prop, fill = .data[[var]])) +
    geom_bar(width = 1, stat = "identity", color = "white") +
    coord_polar("y", start = 0) +
    theme_void() + # Removes axes and grey background
    labs(
      title = paste("Distribution by", var),
      fill = var
    ) +
    # Add text labels in the middle of each slice
    geom_text(aes(label = label_text), 
              position = position_stack(vjust = 0.5), 
              size = 3.5, fontface = "bold") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right"
    )
  
  # Save to list
  pie_plots[[var]] <- p
}

# 4. Combine the 4 plots into one image
# (Requires 'patchwork' library, or use cowplot::plot_grid / gridExtra::grid.arrange)
final_plot <- wrap_plots(pie_plots, ncol = 2)

# 5. Save the result
ggsave("data/Result/cluster/3k_subset/Cell_Demographics_PieCharts.png", 
       final_plot, width = 12, height = 10, bg = "white")

# Display in RStudio
print(final_plot)

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
  p <- DimPlot(`3k_subset_seurat`, group.by = var, label = TRUE, repel = TRUE, reduction = "umap_dim30") +
    ggtitle(paste(var, "on Seurat UMAP"))  + NoLegend()
  
  # 3. Save the plot
  # Construct filename and path
  filename <- paste0("data/Result/cluster/3k_subset/", "UMAP_dim30_", var, ".png")
  
  ggsave(filename = filename, plot = p, width = 10, height = 7.5, dpi = 320)
}

#################
# DE analysis
# Identify markers of each cluster
`3k_subset_marker` <- FindAllMarkers(object = `3k_subset_seurat`)
`3k_subset_marker` <- split(`3k_subset_marker`, `3k_subset_marker`$cluster)

# filtered_markers <- lapply(`3k_subset_marker`, function(df) {
#   subset(df, 
#          pct.1 > 0.1 &           # Criteria A
#            pct.2 > 0.1 &           # Criteria B
#            p_val < 0.05 &          # Criteria C
#            (avg_log2FC > 0.25 | avg_log2FC < -0.25)) # Criteria D
# })



#`3k_subset_marker` <- lapply(`3k_subset_marker`, function(df) {
#  df$p_FC <- (1 - df$p_val) * df$avg_log2FC
#  return(df)
#})
#
#`3k_subset_marker` <- lapply(`3k_subset_marker`, function(df) {
#  df$adj_diff <- abs(df$pct.1 - df$pct.2)/max(c(df$pct.1, df$pct.2))
#  return(df)
#})
#
#`3k_subset_marker` <- lapply(`3k_subset_marker`, function(df) {
#  df$diff_pfc <- df$adj_diff * df$p_FC
#  return(df)
#})
################
# age, sex, roi associated DE AMONG whole dataset
# Goal: Find DE genes for Sex, Age, and ROI across all cells

# 1. Define the variables we want to loop through
variables_to_test <- c("sex", "age_cat", "roi")

# 2. Initialize a list to store results
# New Structure: de_results$Variable_Name (No longer nested by subclass)
global_de_results <- list()

# 3. Main Loop (Iterate through variables only)
for (var in variables_to_test) {
  
  message(paste0("\nProcessing Variable Globally: ", var, " ------------------"))
  
  # Use the full object directly
  # (No subsetting happens here)
  object_to_test <- `3k_subset_seurat`
  
  # Check 1: Does the variable exist in metadata?
  if (!var %in% colnames(object_to_test@meta.data)) {
    message(paste("    Skipping: Variable", var, "not found."))
    next
  }
  
  # Set the identity to the variable of interest
  Idents(object_to_test) <- var
  
  # Check 2: Are there at least 2 groups? 
  unique_groups <- unique(Idents(object_to_test))
  unique_groups <- unique_groups[!is.na(unique_groups)] # Remove NAs
  
  if (length(unique_groups) < 2) {
    message(paste("    Skipping: Only 1 group found for", var, "(", unique_groups, ")"))
    next
  }
  
  # C. Run DE Analysis
  tryCatch({
    markers <- FindAllMarkers(
      object = object_to_test,
      only.pos = TRUE,       # Only look for upregulated genes
      verbose = FALSE, 
      logfc.threshold = 0.1, 
      min.pct = 0.01,
      min.diff.pct = -Inf    # Keep Seurat defaults or your specific params
    )
    
    # D. Add Custom Metrics
    if (nrow(markers) > 0) {
      # Calculate p_FC
      markers$p_FC <- (1 - markers$p_val) * markers$avg_log2FC
      
      # Calculate difference in proportion
      markers$adj_diff <- abs(markers$pct.1 - markers$pct.2) / pmax(markers$pct.1, markers$pct.2)
      
      # Calculate combined score
      markers$diff_pfc <- markers$adj_diff * markers$p_FC
      
      # Save to list (Key is just the variable name now)
      global_de_results[[var]] <- markers
      message(paste("    Found", nrow(markers), "markers."))
      
    } else {
      message("    No significant markers found.")
    }
    
  }, error = function(e) {
    message(paste("    Error running DE for", var, ":", e$message))
  })
}

# Optional: Save the global results
saveRDS(global_de_results, file = "data/Result/global_DE_results.rds")

# plot number of DE
# 1. Initialize empty list for summary data
plot_data_list <- list()

# 2. Loop through the results to extract counts
# 'de_results' is the list you created in the previous step
for (var_name in names(global_de_results)) {
  
  # Get the results dataframe for this variable
  df <- global_de_results[[var_name]]
  
  if (!is.null(df) && nrow(df) > 0) {
    # Count the number of markers per group (stored in 'cluster' column)
    summary_df <- df %>%
      group_by(cluster) %>%
      summarise(Count = n()) %>%
      mutate(Variable = var_name) # Add label for the variable type
    
    plot_data_list[[var_name]] <- summary_df
  }
}

# 3. Combine into one master dataframe
plot_data <- do.call(rbind, plot_data_list)

# 4. Create the Plot
# We use 'facet_wrap' to create separate panels for Sex, Age, and ROI
p <- ggplot(plot_data, aes(x = cluster, y = Count, fill = cluster)) +
  geom_col() + 
  geom_text(aes(label = Count), vjust = -0.5, size = 3.5) + # Add numbers on top
  facet_wrap(~Variable, scales = "free_x", strip.position = "top") + # Separate panels
  theme_classic() +
  labs(
    title = "Number of DE Genes by Criteria",
    subtitle = "Non-neuronal 3k DE result",
    x = "Group",
    y = "Number of DE Genes",
    fill = "Group"
  ) +
  theme(
    legend.position = "none", # Hide legend (colors are self-explanatory)
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    strip.background = element_rect(fill = "lightgrey"), # Style the panel headers
    strip.text = element_text(face = "bold", size = 12)
  )

# 5. Save and Print
ggsave("data/Result/Global_DE_Counts.png", p, width = 10, height = 6)
print(p)

#################
# Subclass-Specific DE Analysis
# Goal: Find DE genes for Sex, Age, and ROI within each Subclass

# 1. Define the groups we want to loop through
target_subclasses <- unique(`3k_subset_seurat`$subclass_label)
variables_to_test <- c("sex", "age_cat", "roi")

# 2. Initialize a list to store results
# Structure: de_results$Subclass_Name$Variable_Name
de_results <- list()

# 3. Main Loop
for (subclass in target_subclasses) {
  
  message(paste0("\nProcessing Subclass: ", subclass, " ------------------"))
  
  # A. Subset the object to just this subclass
  # We use subset() to isolate the cells
  sub_obj <- subset(`3k_subset_seurat`, subset = subclass_label == subclass)
  
  # Initialize list for this subclass
  de_results[[subclass]] <- list()
  
  # B. Loop through the variables (Sex, Age, ROI)
  for (var in variables_to_test) {
    
    message(paste("  Testing variable:", var))
    
    # Check 1: Does the variable exist in metadata?
    if (!var %in% colnames(sub_obj@meta.data)) {
      message(paste("    Skipping: Variable", var, "not found."))
      next
    }
    
    # Set the identity to the variable of interest (e.g., set Idents to "sex")
    Idents(sub_obj) <- var
    
    # Check 2: Are there at least 2 groups? 
    # (e.g., If a subclass only exists in Males, we can't do Male vs Female)
    unique_groups <- unique(Idents(sub_obj))
    unique_groups <- unique_groups[!is.na(unique_groups)] # Remove NAs
    
    if (length(unique_groups) < 2) {
      message(paste("    Skipping: Only 1 group found for", var, "(", unique_groups, ")"))
      next
    }
    
    # C. Run DE Analysis
    # We use FindAllMarkers to handle both binary (M vs F) and multi-class (ROI A vs B vs C) cases automatically.
    # It returns markers for each group compared to the others in this subclass.
    tryCatch({
      markers <- FindAllMarkers(
        object = sub_obj,
        only.pos = TRUE,       # Only look for upregulated genes (standard for markers)
        verbose = FALSE, 
        logfc.threshold = 0.1, 
        min.pct = 0.01, 
        min.diff.pct = -Inf
      )
      
      # D. Add your Custom Metrics (p_FC, etc.)
      if (nrow(markers) > 0) {
        # Calculate p_FC
        markers$p_FC <- (1 - markers$p_val) * markers$avg_log2FC
        
        # Calculate difference in proportion (pct.1 - pct.2)
        # Note: FindAllMarkers returns pct.1 and pct.2 automatically
        markers$adj_diff <- abs(markers$pct.1 - markers$pct.2) / pmax(markers$pct.1, markers$pct.2)
        
        # Calculate combined score
        markers$diff_pfc <- markers$adj_diff * markers$p_FC
        
        # Save to list
        de_results[[subclass]][[var]] <- markers
        message(paste("    Found", nrow(markers), "markers."))
        
      } else {
        message("    No significant markers found.")
      }
      
    }, error = function(e) {
      message(paste("    Error running DE for", var, ":", e$message))
    })
  }
}

saveRDS(de_results, file = "data/Result/DE_results.rds")

#################
# extract top 15 DE gene sorted by adjusted p_value
# ------------------------------------------------------------------
# 1. Base Directory Setup
# ------------------------------------------------------------------
base_dir_top15_global <- file.path("data", "Result", "cluster", "3k_subset", "Top15_DE", "global")
base_dir_top15_subclass <- file.path("data", "Result", "cluster", "3k_subset", "Top15_DE", "subclass")

# ------------------------------------------------------------------
# 2. Extraction and Export Function
# ------------------------------------------------------------------
extract_and_save_top15 <- function(de_df, out_path) {
  
  # Validation: Ensure dataframe exists and is not empty
  if (is.null(de_df) || nrow(de_df) == 0) return(invisible())
  
  # Group by cluster to ensure we get top 15 PER condition
  # Primary sort: p_val_adj (ascending / non-decreasing)
  # Secondary sort: absolute avg_log2FC (descending) to break p-value ties
  top15_df <- de_df %>%
    group_by(cluster) %>%
    arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
    slice_head(n = 15) %>%
    ungroup()
  
  # Write to structured CSV
  write.csv(top15_df, file = out_path, row.names = FALSE)
}

# ------------------------------------------------------------------
# 3. Execute Loop: Global Top 15 DE
# ------------------------------------------------------------------
message("\nExtracting top 15 global DE genes...")

for (var in names(global_de_results)) {
  
  # Construct nested directory: .../Top15_DE/global/<Variable>/
  var_dir <- file.path(base_dir_top15_global, var)
  dir.create(var_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Define file path
  out_file <- file.path(var_dir, paste0("Global_", var, "_top15.csv"))
  
  extract_and_save_top15(global_de_results[[var]], out_file)
}

# ------------------------------------------------------------------
# 4. Execute Loop: Subclass Top 15 DE
# ------------------------------------------------------------------
message("Extracting top 15 subclass DE genes...")

for (subclass in names(de_results)) {
  
  # Enforce string sanitization for directory creation
  safe_subclass <- gsub("[ /]", "_", subclass)
  
  for (var in names(de_results[[subclass]])) {
    
    # Construct exact nested path: .../Top15_DE/subclass/<Subclass_Name>/<Criteria>/
    nested_dir <- file.path(base_dir_top15_subclass, safe_subclass, var)
    dir.create(nested_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Define file path
    out_file <- file.path(nested_dir, paste0(safe_subclass, "_", var, "_top15.csv"))
    
    extract_and_save_top15(de_results[[subclass]][[var]], out_file)
  }
}

message("Top 15 DE extraction mapped to hierarchical directories.")

#################
# 1. Extract the number of DE genes from the 'de_results' list
# Initialize an empty list to store the counts

summary_list <- list()

# Loop through the nested list
for (subclass in names(de_results)) {
  for (variable in names(de_results[[subclass]])) {
    
    # Get the marker dataframe
    markers <- de_results[[subclass]][[variable]]
    
    # Count rows (if NULL, count is 0)
    num_genes <- if (is.null(markers)) 0 else nrow(markers)
    
    # Store in the list
    summary_list[[length(summary_list) + 1]] <- data.frame(
      Subclass = subclass,
      Variable = variable,
      Count = num_genes
    )
  }
}

# 2. Combine into a clean Data Frame
plot_data <- do.call(rbind, summary_list)

# Optional: Clean up variable names for the plot (Capitalize)
plot_data$Variable <- factor(plot_data$Variable, 
                             levels = c("sex", "age_cat", "roi"),
                             labels = c("Sex", "Age", "ROI"))

# 3. Create the Plot
p <- ggplot(plot_data, aes(x = Subclass, y = Count, fill = Variable)) +
  geom_bar(stat = "identity", position = "dodge") + # 'dodge' puts bars side-by-side
  coord_flip() + # Flip coordinates to make subclass names readable
  labs(
    title = "Number of DE Genes per Subclass",
    subtitle = "Comparing Sex, Age, and ROI differences",
    x = "Subclass Label",
    y = "Number of DE Genes",
    fill = "Criteria"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10), # Adjust text size if you have many subclasses
    legend.position = "top"
  ) +
  geom_text(aes(label = Count), position = position_dodge(width = 0.9), hjust = -0.2, size = 3) # Add numbers

# 4. Display and Save
print(p)
ggsave("data/Result/cluster/3k_subset/DE_Counts_Summary.png", p, width = 10, height = 8)

###################
# GO enrichment
# 1. Define a helper function to run GO on a single DE dataframe
run_go_comparison <- function(de_df, organism_db = org.Mm.eg.db) {
  
  # Return NULL if the dataframe is empty
  if (is.null(de_df) || nrow(de_df) == 0) return(NULL)
  
  # Filter for significant genes (Adjust thresholds as needed)
  sig_df <- de_df %>%
    filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
  
  # Return NULL if there aren't enough genes to compare
  if (nrow(sig_df) < 5) {
    message("    Not enough significant genes for GO enrichment.")
    return(NULL)
  }
  
  # Run Biological Process (BP) comparison
  message("    Running GO: Biological Process (BP)...")
  go_bp <- compareCluster(
    gene ~ cluster, 
    data = sig_df, 
    fun = "enrichGO",
    OrgDb = organism_db,
    keyType = "SYMBOL", # Seurat usually outputs Gene Symbols
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  # Run Molecular Function (MF) comparison
  message("    Running GO: Molecular Function (MF)...")
  go_mf <- compareCluster(
    gene ~ cluster, 
    data = sig_df, 
    fun = "enrichGO",
    OrgDb = organism_db,
    keyType = "SYMBOL",
    ont = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  # Return a list containing both results
  return(list(BP = go_bp, MF = go_mf))
}

# ---------------------------------------------------------
# 2. Process global_de_results
# ---------------------------------------------------------
message("\nStarting Global GO Enrichment...")
global_go_results <- list()

for (var in names(global_de_results)) {
  message(paste("\nProcessing Global Variable:", var))
  global_go_results[[var]] <- run_go_comparison(global_de_results[[var]])
}

# Save global GO results
saveRDS(global_go_results, file = "data/Result/global_GO_results.rds")


# ---------------------------------------------------------
# 3. Process nested de_results (Subclass-specific)
# ---------------------------------------------------------
message("\nStarting Subclass-Specific GO Enrichment...")
subclass_go_results <- list()

for (subclass in names(de_results)) {
  message(paste("\n--- Processing Subclass:", subclass, "---"))
  
  subclass_go_results[[subclass]] <- list()
  
  for (var in names(de_results[[subclass]])) {
    message(paste("  Variable:", var))
    
    # Run the helper function and store it in the nested list
    subclass_go_results[[subclass]][[var]] <- run_go_comparison(de_results[[subclass]][[var]])
  }
}

# Save subclass GO results
saveRDS(subclass_go_results, file = "data/Result/subclass_GO_results.rds")
message("\nGO Enrichment complete!")

# GO visualization
# ------------------------------------------------------------------
# 1. Base Directory Setup
# ------------------------------------------------------------------
base_dir_global <- file.path("data", "Result", "cluster", "3k_subset", "GO_plots", "global")
base_dir_subclass <- file.path("data", "Result", "cluster", "3k_subset", "GO_plots", "subclass")

# ------------------------------------------------------------------
# 2. Refactored Plotting Function
# ------------------------------------------------------------------
generate_and_save_plots <- function(go_res_list, title_prefix, file_prefix, out_dir) {
  
  if (is.null(go_res_list)) return(invisible())
  
  for (ont in names(go_res_list)) {
    res <- go_res_list[[ont]]
    
    # Validation constraint
    if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
      
      message(paste("  Generating plots for:", title_prefix, "-", ont))
      
      # -- A. Dotplot --
      p_dot <- dotplot(res, showCategory = 5) + 
        ggtitle(paste(title_prefix, "-", ont)) +
        theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))
      
      # Save to dynamically assigned out_dir
      ggsave(filename = file.path(out_dir, paste0(file_prefix, "_", ont, "_dotplot.png")), 
             plot = p_dot, width = 10, height = 8, bg = "white")
      
      # -- B. Cnetplot --
      tryCatch({
        p_cnet <- cnetplot(res, showCategory = 4, layout = "kk", colorEdge = TRUE, node_label = "all") + 
          ggtitle(paste(title_prefix, "-", ont, "Network")) +
          theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))
        
        ggsave(filename = file.path(out_dir, paste0(file_prefix, "_", ont, "_cnetplot.png")), 
               plot = p_cnet, width = 12, height = 10, bg = "white")
        
      }, error = function(e) {
        message(paste("    [Warning] Failed cnetplot layout for", file_prefix, ont, "-", e$message))
      })
      
    } else {
      message(paste("  [Skip] No significant terms for:", title_prefix, "-", ont))
    }
  }
}

# ------------------------------------------------------------------
# 3. Execute Loop: Global Results
# ------------------------------------------------------------------
message("\nProcessing Global GO visualizations...")

for (var in names(global_go_results)) {
  
  # Construct and execute directory creation: .../global/<Variable>/
  var_dir <- file.path(base_dir_global, var)
  dir.create(var_dir, recursive = TRUE, showWarnings = FALSE)
  
  generate_and_save_plots(
    go_res_list = global_go_results[[var]],
    title_prefix = paste("Global", toupper(var)),
    file_prefix = paste0("Global_", var),
    out_dir = var_dir
  )
}

# ------------------------------------------------------------------
# 4. Execute Loop: Subclass Results (Nested Architecture)
# ------------------------------------------------------------------
message("\nProcessing Subclass-specific GO visualizations...")

for (subclass in names(subclass_go_results)) {
  
  # Critical string sanitization constraint
  safe_subclass <- gsub("[ /]", "_", subclass)
  
  for (var in names(subclass_go_results[[subclass]])) {
    
    # Construct exact nested path: .../subclass/<Subclass_Name>/<Criteria>/
    nested_dir <- file.path(base_dir_subclass, safe_subclass, var)
    
    # Enforce directory creation prior to passing to function
    dir.create(nested_dir, recursive = TRUE, showWarnings = FALSE)
    
    generate_and_save_plots(
      go_res_list = subclass_go_results[[subclass]][[var]],
      title_prefix = paste(subclass, toupper(var)),
      file_prefix = paste0(safe_subclass, "_", var),
      out_dir = nested_dir
    )
  }
}

message("\nAll visualizations mapped to hierarchical directories.")