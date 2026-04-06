library(ggplot2)
library(dplyr)

gene_of_interest <- "Ddx3y"
target_var <- "sex"

# 1. Generate base plot
p_base <- VlnPlot(glia_only_seurat, features = gene_of_interest, 
                  split.by = target_var, group.by = "subclass_label",
                  pt.size = 0.1, combine = FALSE)[[1]]

plot_data <- p_base$data

# 2. Calculate coordinates natively handling sparse zeroes
y_positions <- plot_data %>%
  group_by(ident) %>%
  summarise(
    # suppressWarnings handles the -Inf warning if a group is all zeroes
    max_val = suppressWarnings(max(.data[[gene_of_interest]], na.rm = TRUE)) 
  ) %>%
  mutate(
    # Convert -Inf back to 0 for totally unexpressed groups
    max_val = ifelse(is.infinite(max_val), 0, max_val),
    # Calculate a dynamic buffer (10% of the absolute max expression) instead of a hardcoded 0.5
    y_max = max_val + (max(max_val) * 0.1)
  )

# 3. Retrieve adjusted p-values
sig_list <- list()

for (subclass in levels(plot_data$ident)) {
  
  markers <- de_results[[subclass]][[target_var]]
  
  if (!is.null(markers) && gene_of_interest %in% rownames(markers)) {
    pval <- markers[gene_of_interest, "p_val_adj"]
    star <- ifelse(pval < 0.001, "***",
                   ifelse(pval < 0.01, "**",
                          ifelse(pval < 0.05, "*", "ns")))
  } else {
    star <- "ns" 
  }
  
  sig_list[[length(sig_list) + 1]] <- data.frame(ident = subclass, sig = star)
}

sig_df <- do.call(rbind, sig_list)
anno_df <- merge(y_positions, sig_df, by = "ident")

# 4. Overlay labels and force y-axis expansion
p_final <- p_base +
  geom_text(data = anno_df, aes(x = ident, y = y_max, label = sig), 
            size = 5, fontface = "bold", color = "black", inherit.aes = FALSE) +
  # This command forces the top of the y-axis to expand by an additional 15%
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_final)