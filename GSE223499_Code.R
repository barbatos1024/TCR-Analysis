
rm(list = ls())

getwd()
setwd('/Users/naturexv/Datasets/TCR/GSE223503_NSCLC')

list.files("./GSE223499_RAW/")

library(Seurat)
library(stringr)
library(data.table)
library(dplyr)
library(vegan)
library(ggpubr)
library(ggplot2)

TCR_merge <- read.table("GSE223500_Merged_TCR_data.csv.gz", header = T, sep = ",", row.names = 1)
metadata <- read.table("GSE223499_sn_integrated_data.csv.gz", header = T, sep = ",", row.names = 1)

data_dir <- "./GSE223499_RAW"

files <- list.files(data_dir, pattern = "h5$", full.names = TRUE)
seu.list <- lapply(files, function(f){
  sample_id <- basename(f) %>%
    str_remove("_sn_raw_feature_bc_matrix.h5")
  print(sample_id)
  counts <- Read10X_h5(f)        # 读取10X
  colnames(counts) <- paste0(sample_id, "_", colnames(counts))       # 修改barcode
  keep.cells <- intersect(colnames(counts), rownames(metadata))        # 当前样本metadata中的细胞
  print(paste0(
    sample_id,
    ": ",
    length(keep.cells),
    " cells kept"
  ))
  counts <- counts[, keep.cells]         # 只保留metadata存在的细胞
  seu <- CreateSeuratObject(counts = counts, project = sample_id)     # 创建Seurat
  seu$sample <- sample_id               # 添加metadata
  seu <- AddMetaData(
    seu,
    metadata = metadata[keep.cells, ]    # 添加你的metadata信息
  )
  return(seu)
})

scrna <- merge(seu.list[[1]], y = seu.list[-1], add.cell.ids = NULL)
table(metadata$orig.ident)
table(scrna@meta.data$sample)    # 作者少上传两个样本？ KRAS_13 N561
dim(scrna)        
dim(metadata)

### 添加metadata
common.cells <- intersect(colnames(scrna), rownames(metadata))
length(common.cells)

scrna <- subset(scrna, cells = common.cells)
metadata <- metadata[common.cells, ]
scrna <- AddMetaData(scrna, metadata = metadata)

### 添加TCR
dim(TCR_merge)
common.cells <- intersect(colnames(scrna), rownames(TCR_merge))

length(common.cells)
TCR_meta <- TCR_merge[common.cells, ]
scrna <- AddMetaData(scrna, metadata = TCR_meta)
colnames(scrna@meta.data)

table(scrna$Frequency)
summary(scrna$Frequency)
table(scrna$cloneType)

### IGSF3 levels vs TCR expansion per sample
igsf3_df <- FetchData(scrna, vars=c("IGSF3", "orig.ident"))
igsf3_sample <- igsf3_df %>%
  group_by(orig.ident) %>%
  summarise(
    IGSF3_mean = mean(IGSF3),
    IGSF3_pos_ratio = mean(IGSF3 > 0)
  )

colnames(TCR_meta)
table(TCR_meta$cloneType)
table(TCR_meta$orig.ident)
TCR_meta$expanded <- TCR_meta$cloneType %in% c("Medium (0.001 < X <= 0.01)", "Large (0.01 < X <= 0.1)", "Hyperexpanded (0.1 < X <= 1)")
expansion_sample <- TCR_meta %>%
  group_by(orig.ident) %>%
  summarise(
    total_TCR_cells = n(),
    expanded_cells = sum(expanded),
    expanded_fraction = expanded_cells / total_TCR_cells
  )
cor_df <- merge(igsf3_sample, expansion_sample, by = "orig.ident")

# IGSF3表达水平 vs TCR expansion
cor.test(cor_df$IGSF3_mean, cor_df$expanded_fraction, method = "spearman")

p1 = ggplot(cor_df, aes(x = IGSF3_mean, y = expanded_fraction)) +
  geom_point(size = 2.5, shape = 1) +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., sep = "~`  `~")),
           label.x = 0.3, label.y = 0.4, size = 4) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title.x = element_text(family = "Arial", size = 13),
        axis.title.y = element_text(family = "Arial", size = 13),
        axis.text.x = element_text(family = "Arial", size = 11),
        axis.text.y = element_text(family = "Arial", size = 11)) +
  labs(x = "IGSF3 mean expression levels", y = "TCR expanded fraction")         

cor.test(cor_df$IGSF3_pos_ratio, cor_df$expanded_fraction, method = "spearman")

p2 = ggplot(cor_df, aes(x = IGSF3_pos_ratio, y = expanded_fraction)) +
  geom_point(size = 2.5, shape = 1, ) +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., sep = "~`  `~")),
           label.x = 0.15, label.y = 0.4, size = 4) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title.x = element_text(family = "Arial", size = 13),
        axis.title.y = element_text(family = "Arial", size = 13),
        axis.text.x = element_text(family = "Arial", size = 11),
        axis.text.y = element_text(family = "Arial", size = 11)) +
  labs(x = "IGSF3 positive cells ratio", y = "TCR expanded fraction")         # IGSF3表达水平 vs TCR expansion

p1 | p2
# ggsave(filename = "TCR_correlation.png", width = 6, height = 3, dpi = 300)

### TCR diversity
clone_id <- "CTstrict"
tcr_count <- TCR_meta %>%
  filter(!is.na(CTstrict)) %>%
  group_by(orig.ident, CTstrict) %>%
  summarise(
    clone_size = n(),
    .groups="drop"
  )
tcr_freq <- tcr_count %>%
  group_by(orig.ident) %>%
  mutate(
    proportion = clone_size / sum(clone_size)
  )
shannon <- tcr_freq %>%
  group_by(orig.ident) %>%
  summarise(
    Shannon = diversity(clone_size, index="shannon")
  )
simpson <- tcr_freq %>%
  group_by(orig.ident) %>%
  summarise(
    Simpson = diversity(clone_size,index="simpson")
  )
TCR_diversity <- left_join(shannon, simpson, by="orig.ident")

cor_df_tcr <- merge(igsf3_sample, TCR_diversity, by="orig.ident")
p3 = ggplot(cor_df_tcr, aes(x = IGSF3_mean, y = Shannon)) +
  geom_point(size = 2.5, shape = 1) +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., sep = "~`  `~")),
           label.x = 0.1, label.y = -1, size = 4) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title.x = element_text(family = "Arial", size = 13),
        axis.title.y = element_text(family = "Arial", size = 13),
        axis.text.x = element_text(family = "Arial", size = 11),
        axis.text.y = element_text(family = "Arial", size = 11)) +
  labs(x = "IGSF3 mean expression levels", y = "TCR Diversity")        
p3

p4 = ggplot(cor_df_tcr, aes(x = IGSF3_pos_ratio, y = Shannon)) +
  geom_point(size = 2.5, shape = 1) +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., sep = "~`  `~")),
           label.x = 0.1, label.y = -1, size = 4) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title.x = element_text(family = "Arial", size = 13),
        axis.title.y = element_text(family = "Arial", size = 13),
        axis.text.x = element_text(family = "Arial", size = 11),
        axis.text.y = element_text(family = "Arial", size = 11)) +
  labs(x = "IGSF3 positive cells ratio", y = "TCR Diversity")        
p4
p3 | p4

p1 | p3
# ggsave(filename = "TCR_expansion&diversity.png", width = 6.5, height = 3.25, dpi = 300)

