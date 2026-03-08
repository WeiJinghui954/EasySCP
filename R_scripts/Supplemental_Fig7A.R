###########################
#### SCP Figure Code ####
###########################

#### -- Supplemental_Fig7A -- ####


#  Prepare Workspace
rm(list = ls());gc()

packages <- c(
  "readxl",       
  "VennDiagram", 
  "dplyr",
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
load("../Output/Variables/SCP_protein_data_sorted_by_com.R")
SCP_zonated_protein <- rownames(protein_data_sorted_by_com)

load("../Data/Guilliams_2022/ST_Visium_Zonated.RData")
Guilliams_2022_zonated_gene <- ST_Visium_union_zonated_gene

load("../Data/Xu_2024/ST_Stereoseq_Zonated.RData")
Xu_2024_zonated_gene <- ST_Stereoseq_union_zonated_gene

load("../Output/Variables/Rosenberger_2023_zonated_protein.R")
Rosenberger_2023_zonated_protein <-na.omit(Rosenberger_2023_zonated_protein) 

load("../Output/Variables/Ben_Moshe_2019_zonated_protein.R")
Ben_Moshe_2019_zonated_protein <- Ben_Moshe_2019_zonated_protein

#  Save variables
conserved_zonation_markers <- Reduce(intersect, list(
  SCP_zonated_protein, Rosenberger_2023_zonated_protein,
  Guilliams_2022_zonated_gene,Xu_2024_zonated_gene,Ben_Moshe_2019_zonated_protein
))
save(conserved_zonation_markers,file = "../Output/Variables/Conserved_zonation_markers.R")


#  Save plot
flog.threshold(ERROR)
venn.plot <- venn.diagram(
  x = list(SCP_zonated_protein, Rosenberger_2023_zonated_protein,
           Guilliams_2022_zonated_gene,Xu_2024_zonated_gene,Ben_Moshe_2019_zonated_protein), 
  category.names = c("SCP","Rosenberger_2023", 
                     "Guilliams_2022","Xu_2024","Ben_Moshe_2019"),
  filename = NULL,   # 设置为 NULL 可以在 R 中显示图形
  output = TRUE,
  col = "black",     # Venn图边框颜色
  fill = c("#fc636b","#ffb900", "#6a67ce","#1aafd0","#3be8b0"), # 填充颜色
  alpha = 0.6,       # 填充透明度
  cat.col = c("#fc636b","#ffb900", "#6a67ce","#1aafd0","#3be8b0"),  # 分类标签颜色
  cat.cex = 1.5,     # 分类标签字体大小
  cex = 2,         # 数字字体大小
  height = 480,      # 高度
  width = 480,       # 宽度
  resolution = 300,  # 分辨率
  layout = "circle",  # 设置布局为圆形
  fontfamily = "Arial",      # 数字字体 Arial
  scaled = TRUE,      # 确保两个圆的大小相等
)
file_name <- paste0("../Output/Figure/Supplemental Fig7A.png")
ggsave(file_name, venn.plot, width = 10, height = 10)






