###########################
#### SCP Figure Code ####
###########################

#### -- Figure 5C -- ####


#  Prepare Workspace
rm(list = ls());gc()

packages <- c(
  "tidyverse",
  "Seurat",
  "GSVA",
  "ggplot2"
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
load("../Output/Variables/SCP_filtered_data.R")
load("../Output/Variables/Conserved_zonation_markers.R")
load("../Output/Variables/SCP_HSS_gene_sets.R")

#  SCP HSS score
expr_matrix_SCP_normalized <- filtered_data[ , !grepl("20241220", colnames(filtered_data))]
expr_matrix_SCP_normalized <- expr_matrix_SCP_normalized[rownames(expr_matrix_SCP_normalized) %in% conserved_zonation_markers, ]

expr_matrix <- as.matrix(expr_matrix_SCP_normalized)
params <- gsvaParam(
  expr_matrix,                 
  HSS_gene_sets,         
  kcdf = "Gaussian"                 
)
scores <- gsva(params)
HSS_score <- scores["CV", ] - scores["PV", ]
sorted_scores <- sort(HSS_score, decreasing = F)
sorted_scores_matrix <- cbind(Cell = names(sorted_scores),
                              Score = as.numeric(sorted_scores))

#  Save plot
df <- as.data.frame(sorted_scores_matrix)
df$Score <- as.numeric(df$Score)
inner_gap <- 0.1   
df$norm_score <- (df$Score - min(df$Score)) / (max(df$Score) - min(df$Score))
df$r <- inner_gap + (1 - df$norm_score) * (1 - inner_gap) 
set.seed(123)
df$theta <- runif(nrow(df), 0, 2*pi)
df$x <- df$r * cos(df$theta)
df$y <- df$r * sin(df$theta)

inside_hex <- function(x,y,R=1){
  (abs(x) <= R) & (abs(y) <= sqrt(3)*(R - abs(x)/2))
}
df <- df[inside_hex(df$x, df$y), ]
hex_border <- function(R=1){
  angles <- seq(0, 2*pi, length.out=7)
  data.frame(x=R*cos(angles), y=R*sin(angles), group=R)
}

n_layers <- 8
hex_list <- lapply(seq(inner_gap, 1, length.out=n_layers), hex_border)
hex_df <- do.call(rbind, hex_list)
p <- ggplot() +
  geom_polygon(data=hex_df, aes(x=x, y=y, group=group), 
               fill=NA, color="grey60", size=0.6, linetype="dashed") +
  geom_point(data=df, aes(x=x, y=y, color=Score), size=4, alpha=0.9) +
  scale_color_gradientn(colors = c("#FDE725", "#B5DE2B", "#28AE80", "#20908C", "#440154")) +
  coord_equal() +
  theme_void() +
  labs(color="HSS") +
  theme(plot.title = element_text(hjust=0.5, size=16, face="bold"))

ggsave(filename = "../Output/Figure/Fig5C.png", plot = p, width = 7, height = 6,bg = "white",limitsize = FALSE)







