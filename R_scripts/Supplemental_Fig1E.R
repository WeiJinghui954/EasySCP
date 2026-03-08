###########################
#### SCP Figure Code ####
###########################

#### -- Supplemental Figure 1E -- ####

#  Prepare Workspace
rm(list = ls());gc()

packages <- c(
  "Seurat",
  "ggplot2",
  "dplyr"
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
load("../Output/Variables/HEK293T_50cells_data.R")

# CV
data_CV <- data_1_numeric %>%
  rowwise() %>%
  mutate(
    mean_value = mean(c_across(everything()), na.rm = TRUE),
    sd_value = sd(c_across(everything()), na.rm = TRUE),
    CV = ifelse(mean_value != 0, (sd_value / mean_value) * 100, NA)
  ) %>%
  ungroup() %>%
  select(CV)
density_data <- density(data_CV$CV, na.rm = TRUE)
median_CV <- median(data_CV$CV, na.rm = TRUE)

#  Save plot
CV_plot <- ggplot(data_CV, aes(x = CV)) +
  geom_density(fill = "#037ef3", color = "black", alpha = 0.3) +
  geom_vline(xintercept = median_CV, color = "red", linetype = "dashed", linewidth = 0.6) +
  annotate(
    "text", 
    x = median_CV + 20, 
    y = max(density_data$y) * 0.9,  
    label = sprintf("%.2f", median_CV),  
    color = "red", 
    vjust = -1,
    size = 6
  ) +
  theme_minimal() +
  labs(
    x = "CV (%)",
    y = "Density"
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + 
  coord_cartesian(xlim = c(0, max(data_CV$CV)), 
                  ylim = c(0, max(density_data$y) + 0.02)) + 
  theme(
    axis.text.x = element_text(size = 14, angle = 0, hjust = 1),  
    axis.text.y = element_text(size = 14),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16),  
    plot.title = element_text(size = 16, hjust = 0.5),  
    panel.grid = element_blank(),  
    axis.line = element_line(color = "black")  
  )
ggsave(file = paste0("../Output/Figure/Supplemental Fig1E.png"), plot = CV_plot, width = 6, height = 6, bg = "white")






