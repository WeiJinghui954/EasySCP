###########################
#### SCP Figure Code ####
###########################

#### -- Data process -- ####

#  Prepare Workspace
rm(list = ls());gc()

#  Rosenberger_2023 data process
## Read relevant data
load("../Data/Rosenberger_2023/d.R")
load("../Data/Rosenberger_2023/meta_distances.R")
load("../Data/Rosenberger_2023/meta_pg.R")
load("../Data/Rosenberger_2023/SA_incl_all.R")

img_meta <- read_csv("../Data/Rosenberger_2023/meta_img-proteome.csv") %>%
  mutate(cell_ID = str_replace(cell_ID, "DimethNter", "target"))

## Define number of classes
classes = 20

## Subset to 90% complete proteins
SA_incl_heps <- d %>%
  filter(cell_ID %in% meta_distances$cell_ID) %>%
  distinct(cell_ID) %>%
  pull(cell_ID)

data.frame(cell_ID = meta_distances$cell_ID, ratio = meta_distances$ratio) %>%
  mutate(range = cut_interval(ratio, n = classes))  -> meta_distances_bins

meta_distances_bins %>%
  filter(cell_ID %in% SA_incl_heps) %>%
  distinct(range) %>%
  arrange(range) %>%
  mutate(bin = c(1:classes)) %>%
  right_join(meta_distances_bins) %>%
  filter(cell_ID %in% SA_incl_heps) %>%
  column_to_rownames("cell_ID") %>%
  mutate(bin = abs(bin - (classes + 1))) -> meta_distances_bins

d %>%
  filter(cell_ID %in% rownames(meta_distances_bins)) %>%
  drop_na(int_core) %>%
  group_by(Protein) %>%
  summarise(completeness = n()/length(SA_incl_heps)) %>%
  filter(completeness >= 0.7) %>%
  pull(Protein) -> proteome_90_heps

d %>%
  filter(cell_ID %in% SA_incl_heps) %>%
  filter(Protein %in% proteome_90_heps) %>%
  filter(cell_ID %in% rownames(meta_distances_bins)) %>%
  dplyr::select(cell_ID, int_core, Protein) %>%
  spread(cell_ID, int_core) %>%
  arrange(Protein) %>%
  column_to_rownames("Protein") -> d_wide_90

## Limma statistics
design <- model.matrix(~meta_distances_bins[colnames(d_wide_90),]$bin)

fit <- lmFit(d_wide_90, design)
fit <- eBayes(fit)
limma_8bins_90complete <- topTable(fit, number = Inf, confint = TRUE, coef = 2, adjust.method = "fdr") %>%
  rownames_to_column("Protein") %>%
  left_join(meta_pg) %>%
  arrange(logFC) %>%
  mutate(FC_rank = c(1:nrow(.))) %>%
  mutate(significant = adj.P.Val < 0.05)

d %>%
  dplyr::select(cell_ID, int_core, Protein) %>%
  left_join(meta_distances_bins %>% rownames_to_column("cell_ID") %>% drop_na(bin)) %>%
  group_by(Protein, bin) %>%
  summarise(int = log2(median(2^int_core))) %>%
  drop_na(bin) %>%
  spread(bin, int) %>%
  column_to_rownames("Protein") -> d_wide_90

#  Heatmappping data
d_heatmap <- d_wide_90[limma_8bins_90complete %>% arrange(logFC) %>% pull(Protein),]
heatmap_proteins <- rownames(d_heatmap)
heatmap_symbols <- d$Symbol[match(heatmap_proteins, d$Protein)]
Rosenberger_2023_zonated_protein <- heatmap_symbols

#  Save variables
save(Rosenberger_2023_zonated_protein,file = "../Output/Variables/Rosenberger_2023_zonated_protein.R")


#  Ben_Moshe_2019 data process
file_path <- "../Data/Ben_Moshe_2019/42255_2019_109_MOESM6_ESM.xlsx"
protein_data <- readxl::read_excel(file_path)
colnames(protein_data) <- as.character(protein_data[1, ])
protein_data <- protein_data[-1, ]
protein_data <- protein_data %>%
  mutate(`Protein FDR` = as.numeric(`Protein FDR`))

Ben_Moshe_2019_zonated_protein <- protein_data %>%
  filter(`Protein FDR` < 0.05) %>%
  pull(Gene)

#  Save variables
save(Ben_Moshe_2019_zonated_protein,file = "../Output/Variables/Ben_Moshe_2019_zonated_protein.R")








