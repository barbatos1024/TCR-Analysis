
rm(list = ls())

getwd()
setwd('/Users/naturexv/Datasets/TCR/GSE286827_HNSCC')

list.files("./GSE286827_RAW/")

library(Seurat)
library(Matrix)
library(dplyr)
library(stringr)
library(ggpubr)

metadata_Bcells <- readRDS("GSE286827_metadata_Bcells.rds")
metadata_CD4T_TCR <- readRDS("GSE286827_metadata_CD4T_TCR.rds")
metadata_CD8T_TCR <- readRDS("GSE286827_metadata_CD8T_TCR.rds")
metadata_immunecells <- readRDS("GSE286827_metadata_immunecells_broad.rds")
metadata_Malignant <- readRDS("GSE286827_metadata_Malignant.rds")
metadata_TAMs <- readRDS("GSE286827_metadata_TAMs.rds")
metadata_TsCD8T <- readRDS("GSE286827_metadata_TumorSpecific_CD8T.rds")

metadata_immunecells[1:5,]
metadata_Malignant[1:5,]

metadata_Malignant <- metadata_Malignant %>%
  dplyr::rename(manual_celltype_broad = pan_cancer_subcluster)                  # 修改列名

metadata <- dplyr::bind_rows(metadata_immunecells, metadata_Malignant)      # 上下拼接
metadata[1:3,]


data_dir <- "./GSE286827_RAW/"

samples <- list.dirs(data_dir, full.names = FALSE, recursive = FALSE)

seu.list <- lapply(samples, function(sample_id) {
  message("Processing: ", sample_id)
  counts <- Read10X(data.dir = file.path(data_dir, sample_id))         # 读取10X

  if (is.list(counts)) {
    counts <- counts[["Gene Expression"]]
  }

  old_bc <- colnames(counts)             # 原始barcode
  old_bc_clean <- sub("-1$", "", old_bc)         # 去掉10X barcode末尾的 -1
  new_bc <- paste0(sample_id, "_", old_bc_clean)        # 加sample ID
  colnames(counts) <- new_bc
  keep_cells <- intersect(new_bc, rownames(metadata))      # 和metadata匹配
  
  message(
    "10X cells: ", length(new_bc),
    " | metadata cells: ",
    sum(grepl(
      paste0("^", sample_id, "_"),
      rownames(metadata)
    )),
    " | matched: ", length(keep_cells)
  )
  
  # 如果没有匹配，直接停止
  if (length(keep_cells) == 0) {
    stop(
      "No cells matched for ",
      sample_id,
      "\n请检查barcode格式。"
    )
  }
  
  counts <- counts[, keep_cells, drop = FALSE]        # 保留metadata存在的细胞
  meta <- metadata[keep_cells, , drop = FALSE]          # 提取metadata
  meta <- meta[colnames(counts), , drop = FALSE]         # 保证顺序完全一致
  seu <- CreateSeuratObject(counts = counts, meta.data = meta, project = sample_id, min.cells = 0, min.features = 0)     # 创建Seurat
  
  return(seu)
})

names(seu.list) <- samples
scrna <- merge(x = seu.list[[1]], y = seu.list[-1])
rm(seu.list)
mem.maxVSize(vsize = 1e9)
scrna@assays[["RNA"]] <- JoinLayers(scrna@assays[["RNA"]])

sum(rownames(metadata_CD8T_TCR) %in% colnames(scrna))
sum(rownames(metadata_CD4T_TCR) %in% colnames(scrna))

metadata_CD8T_TCR <- metadata_CD8T_TCR[rownames(metadata_CD8T_TCR) %in% colnames(scrna), , drop = FALSE]
metadata_CD4T_TCR <- metadata_CD4T_TCR[rownames(metadata_CD4T_TCR) %in% colnames(scrna), , drop = FALSE]

table(metadata_CD4T_TCR$batch)
table(metadata_CD8T_TCR$batch)

colnames(metadata_CD4T_TCR)
colnames(metadata_CD8T_TCR)

colnames(metadata_CD4T_TCR)[10] <- "subcluster"
colnames(metadata_CD8T_TCR)[10] <- "subcluster"

metadata_CD8T_TCR$expansion_status <- NULL

colnames(metadata_CD4T_TCR)
colnames(metadata_CD8T_TCR)

metadata_CD4T_TCR$cell_type <- "CD4"
metadata_CD8T_TCR$cell_type <- "CD8"

metadata_T_TCR <- rbind(metadata_CD4T_TCR, metadata_CD8T_TCR)
colnames(metadata_T_TCR)

table(metadata_T_TCR$batch)

clonotype_expansion <- metadata_T_TCR %>%
  filter(!is.na(CTstrict), CTstrict != "") %>%
  group_by(batch, CTstrict) %>%
  summarise(
    clone_size = n(),
    .groups = "drop"
  ) %>%
  group_by(batch) %>%
  mutate(
    clone_freq = clone_size / sum(clone_size)
  ) %>%
  ungroup()



clonotype_expansion <- clonotype_expansion %>%
  mutate(
    expansion_status = if_else(
      clone_size >= 2,
      "expansion",
      "non-expansion"
    )
  )


expansion_ratio <- clonotype_expansion %>%
  group_by(batch) %>%
  summarise(
    expansion_cells = sum(
      clone_size[clone_size >= 2]
    ),
    
    total_cells = sum(clone_size),
    
    expansion_ratio = expansion_cells / total_cells,
    
    .groups = "drop"
  )



### CD274 expression levels
CD274_sample <- FetchData(
  scrna,
  vars = c("CD274", "batch")
) %>%
  group_by(batch) %>%
  summarise(
    CD274_mean = mean(CD274, na.rm = TRUE),
    CD274_pos_ratio = mean(CD274 > 0, na.rm = TRUE),
    .groups = "drop"
  )
names(table(CD274_sample$batch)) %in% names(table(scrna@meta.data$batch))


CD274_sample <- CD274_sample %>%
  left_join(
    expansion_ratio,
    by = "batch"
  )

CD274_sample <- CD274_sample %>%
  mutate(
    Time = case_when(
      str_ends(batch, "_pre") ~ "pre",
      str_ends(batch, "_post") ~ "post"
    )
  )

ggplot(CD274_sample, aes(x = CD274_mean, y = expansion_ratio)) +
  geom_point(size = 2.5, shape = 1) +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., sep = "~`  `~")),
           label.x = 0.05, label.y = 1, size = 4) + 
  theme_bw() + xlim(0, 0.2) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title.x = element_text(family = "Arial", size = 13),
        axis.title.y = element_text(family = "Arial", size = 13),
        axis.text.x = element_text(family = "Arial", size = 11),
        axis.text.y = element_text(family = "Arial", size = 11)) +
  labs(x = "CD274 mean levels", y = "TCR expansion fraction")     

