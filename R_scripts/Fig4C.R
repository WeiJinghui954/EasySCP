###########################
#### SCP Figure Code ####
###########################

#### -- Figure 4C -- ####

#  Prepare Workspace
rm(list = ls());gc()

packages <- c(
  "tidyverse",
  "ggplot2",
  "ggrepel",
  "patchwork"
  )
for (p in packages) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p)
    library(p, character.only = TRUE)
  }
}

#  Set working Directory
set.seed(42)

#  Load data
load("../Output/Variables/SCP_protein_ordered_com_values.R")

load("../Output/Variables/SCP_protein_data_sorted_by_com.R")
SCP_zonated_protein <- rownames(protein_data_sorted_by_com)

load("../Output/Variables/Rosenberger_2023_zonated_protein.R")
Rosenberger_2023_zonated_protein <- Rosenberger_2023_zonated_protein

load("../Output/Variables/Ben_Moshe_2019_zonated_protein.R")
Ben_Moshe_2019_zonated_protein <- Ben_Moshe_2019_zonated_protein

#  Calculate the intersection
only_in_SCP <- setdiff(
  SCP_zonated_protein,
  union(Ben_Moshe_2019_zonated_protein, Rosenberger_2023_zonated_protein)
)
SCP_Ben_Moshe <- intersect(SCP_zonated_protein, Ben_Moshe_2019_zonated_protein)
SCP_Rosenberger <- intersect(SCP_zonated_protein, Rosenberger_2023_zonated_protein)

ordered_com_values$protein <- factor(ordered_com_values$protein, levels = ordered_com_values$protein)
ordered_com_values$index   <- seq_len(nrow(ordered_com_values))

#  Define groups
ordered_com_values$group <- "Other"
ordered_com_values$group[rownames(ordered_com_values) %in% only_in_SCP] <- "only_in_SCP"
ordered_com_values$group[rownames(ordered_com_values) %in% SCP_Ben_Moshe]   <- "SCP_Ben_Moshe"
ordered_com_values$group[rownames(ordered_com_values) %in% SCP_Rosenberger]     <- "SCP_Rosenberger"

#  Save plot

highlight_genes <- c(
  "Glul", "Gulo","Cyp2e1", "Rhbg", "Slc13a3", "Oat", 
  "Ftcd","Cyp2f2", "Tkfc","Rida", 
  "Slc25a19","Plcb1","Psmf1", "Ctsz","Pdhb",
  "Nt5e","Cyp1a2", "Aldh3a2","Preb",            
  "Cdh1","Gls2","Tomm40l","Gldc"
) %>% unique()
ordered_com_values$label <- ifelse(rownames(ordered_com_values) %in% highlight_genes,
                                   rownames(ordered_com_values), NA)

group_colors <- c(
  "SCP_Ben_Moshe"    = "#1f77b4", 
  "only_in_SCP" = "#ff7f0e", 
  "SCP_Rosenberger"      = "#2ca02c"
)

scatter_plot <- ggplot(ordered_com_values, aes(x = index, y = COM)) +
  geom_point(aes(color = group), size = 1) +
  geom_point(
    data = subset(ordered_com_values, !is.na(label)),
    aes(x = index, y = COM, color = group),
    size = 1.5, stroke = 1.2
  ) +
  geom_text_repel(
    aes(label = label, colour = group),
    size = 5, max.overlaps = Inf, box.padding = 0.3, fontface = "bold"
  ) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.title.x  = element_blank(),
    panel.grid    = element_blank(),
    axis.title.y  = element_text(size = 18),
    axis.text.y   = element_text(color = "black", size = 14),
    axis.line     = element_line(color = "black", size = 0.6),
    axis.ticks    = element_line(color = "black"),
    panel.border  = element_rect(color = "black", fill = NA, size = 1),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title    = element_text(hjust = 0.5, size = 16),
    legend.title  = element_blank(),
    legend.text   = element_text(size = 14),
    plot.margin   = margin(0, 1, 0, 1)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  ylab("COM")

focus_groups <- c("SCP_Ben_Moshe", "SCP_Rosenberger", "only_in_SCP")
group_y_map <- setNames(seq(0, -1, length.out = length(focus_groups) + 1), focus_groups)

bar_df <- ordered_com_values %>% filter(group %in% focus_groups) %>%
  mutate(
    ymin = sapply(group, function(x) group_y_map[x]),
    ymax = ymin + 1 / length(focus_groups)
  )

color_bar_plot <- ggplot(ordered_com_values, aes(x = index, y = COM)) +
  geom_rect(
    data = bar_df,
    aes(xmin = index - 0.4, xmax = index + 0.4, ymin = ymin, ymax = ymax, fill = group),
    inherit.aes = FALSE, show.legend = FALSE
  ) +
  scale_fill_manual(values = group_colors) +
  theme_minimal() +
  theme(
    axis.text   = element_blank(),
    axis.ticks  = element_blank(),
    axis.title  = element_blank(),
    panel.grid  = element_blank(),
    panel.border= element_rect(color = "black", fill = NA, size = 1),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(0, 1, 0, 1)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)))
combined_plot <- scatter_plot / color_bar_plot + plot_layout(heights = c(3, 1))

ggsave(filename = "../Output/Figure/Fig4C.png", plot = combined_plot, width = 8, height = 6,
       dpi = 300,limitsize = FALSE)
