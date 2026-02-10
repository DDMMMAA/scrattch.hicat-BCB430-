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

##############
# load current r project as package
load_all()

##############
# load data
# seurat object
`3k_subset_seurat` <- readRDS("D:/BCB430/scrattch.hicat-BCB430-/data/Allen/3k_subset_seurat.rds")
# expression matrix
`3k_subset_expression_matrix` <- readRDS("D:/BCB430/scrattch.hicat-BCB430-/data/Allen/3k_subset_expression_matrix.rds")

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
all.markers <- FindAllMarkers(object = `3k_subset_seurat`)
all.markers <- split(all.markers, all.markers$cluster)

# filtered_markers <- lapply(all.markers, function(df) {
#   subset(df, 
#          pct.1 > 0.1 &           # Criteria A
#            pct.2 > 0.1 &           # Criteria B
#            p_val < 0.05 &          # Criteria C
#            (avg_log2FC > 0.25 | avg_log2FC < -0.25)) # Criteria D
# })

all.markers <- lapply(all.markers, function(df) {
  df$p_FC <- (1 - df$p_val) * df$avg_log2FC
  return(df)
})

all.markers <- lapply(all.markers, function(df) {
  df$adj_diff <- abs(df$pct.1 - df$pct.2)/max(c(df$pct.1, df$pct.2))
  return(df)
})

all.markers <- lapply(all.markers, function(df) {
  df$diff_pfc <- df$adj_diff * df$p_FC
  return(df)
})

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
        verbose = FALSE
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