
rm(list = ls())

getwd()
setwd("/Users/naturexv/Datasets/TCR/GSE156728_TCR/")

list.files("./GSE156728_CD4_CD8/")         # Chol为SmartSeq, 其余为10X

library(Seurat)
library(dplyr)
library(ggpubr)
library(Matrix)

metadata <- read.table("GSE156728_metadata.txt.gz", header = T, sep = "\t", row.names = 1, check.names = F)

list.files("./GSE156728_CD4_CD8/")

dir <- "./GSE156728_CD4_CD8/"

files <- list.files(dir, pattern = "_CD[48]\\.counts\\.txt\\.gz$", full.names = TRUE)   # 获取除 CHOL所有 counts 文件
files <- files[!grepl("CHOL", files)]

count_list <- list()

for (f in files) {
  cat("Reading:", basename(f), "\n")

  counts <- read.table(f, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)     # 读取 counts
  cell_ids <- colnames(counts)             # 当前文件的 cell ID
  keep_cells <- intersect(cell_ids, rownames(metadata))      # 只保留 metadata 中存在的 cell
  
  cat("Total cells:", length(cell_ids), "| Matched:", length(keep_cells), "\n")
  
  counts <- counts[, keep_cells, drop = FALSE]       # 保留匹配细胞
  count_list[[basename(f)]] <- counts        # 保存
}

mem.maxVSize(vsize = 1e9)
seu_list <- lapply(names(count_list), function(x) {
  counts <- count_list[[x]]
  common_cells <- intersect(colnames(counts), rownames(metadata))        # 只保留 count 和 metadata 都存在的细胞
  
  cat(basename(x), " | count cells =", ncol(counts), " | matched =", length(common_cells), "\n")

  counts <- counts[, common_cells, drop = FALSE]        # 保留共同细胞
  meta <- metadata[common_cells, , drop = FALSE]          # 提取对应 metadata，并保证顺序与 counts 完全一致
  cancer <- sub("_CD[48]\\.counts\\.txt\\.gz$", "", basename(x))        # 从文件名提取癌种和 CD4/CD8
  celltype <- sub("^.*_(CD4|CD8)\\.counts\\.txt\\.gz$", "\\1", basename(x))

  seu <- CreateSeuratObject(counts = counts, meta.data = meta, project = cancer, min.cells = 0, min.features = 0)
  seu$Tcell_type <- celltype   # 添加 CD4/CD8 信息
return(seu)
})

names(seu_list) <- names(count_list)
scrna <- merge(x = seu_list[[1]], y = seu_list[-1])
rm(counts, count_list)
scrna@assays[["RNA"]] <- JoinLayers(scrna@assays[["RNA"]])


library(scRepertoire)
TCR <- read.table("GSE156728_10X_VDJ.merge.txt.gz", header = T, sep = "\t", check.names = F)
TCR_list <- split(TCR, TCR$library.id)
combined.TCR <- combineTCR(
  input.data = TCR_list,
  removeNA = FALSE,
  removeMulti = FALSE,
  filterMulti = FALSE,
  filterNonproductive = FALSE
)
names(combined.TCR) <- names(TCR_list)
table(scrna@meta.data$libraryID)
for (i in seq_along(combined.TCR)) {
  combined.TCR[[i]]$libraryID <- names(combined.TCR)[i]
}       # 在每个combine.TCR里都加入libraryID
scrna <- combineExpression(combined.TCR, 
                           scrna, 
                           cloneCall = "CTstrict", 
                           group.by = "libraryID", 
                           proportion = TRUE)
table(scrna@meta.data$cloneSize)
table(scrna@meta.data$libraryID)

### Hyperexpanded, Large, Medium as expansion
expansion_summary <- scrna@meta.data %>%
  group_by(libraryID) %>%
  summarise(
    Hyperexpanded = sum(cloneSize == "Hyperexpanded (0.1 < X <= 1)", na.rm = TRUE) / n(),
    Large = sum(cloneSize == "Large (0.01 < X <= 0.1)", na.rm = TRUE) / n(),
    Medium = sum(cloneSize == "Medium (0.001 < X <= 0.01)", na.rm = TRUE) / n(),
    .groups = "drop"
  ) %>%
  mutate(
    Hyperexpanded = Hyperexpanded * 100,
    Large = Large * 100,
    Medium = Medium * 100
  )
head(expansion_summary)

expansion_summary <- expansion_summary %>%
  mutate(expansion = Hyperexpanded + Large + Medium)

bulkrna <- read.table("GSE156728_bulk_RNA.counts.txt.gz", header = T, row.names = 2)
bulkrna <- bulkrna[,-1]

colnames(bulkrna) <- gsub("\\.normal$", "-N", colnames(bulkrna))
colnames(bulkrna) <- gsub("\\.tumor$", "-T", colnames(bulkrna))
colnames(bulkrna) <- gsub("\\.", "-", colnames(bulkrna))

colnames(bulkrna)
table(expansion_summary$libraryID)

igsf3_df <- data.frame(
  libraryID = colnames(bulkrna),
  IGSF3 = as.numeric(bulkrna["IGSF3", ])
)             # 提取 IGSF3 表达

expansion_summary <- expansion_summary %>%
  left_join(igsf3_df, by = "libraryID")         # 合并到 expansion_summary

expansion_summary <- expansion_summary %>%
  mutate(
    IGSF3_log = log1p(IGSF3)
  )

cor.test(
  expansion_summary$IGSF3_log,
  expansion_summary$expansion,
  method = "spearman"
)

p1 = ggplot(expansion_summary, aes(x = IGSF3_log, y = expansion)) +
  geom_point(size = 2.5, shape = 1) +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., sep = "~`  `~")),
           label.x = 5, label.y = 80, size = 4) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title.x = element_text(family = "Arial", size = 13),
        axis.title.y = element_text(family = "Arial", size = 13),
        axis.text.x = element_text(family = "Arial", size = 11),
        axis.text.y = element_text(family = "Arial", size = 11)) +
  labs(x = "IGSF3 mean expression levels", y = "TCR expansion fraction")        
p1

### clonalFrequency for expansion calculation 
colnames(scrna@meta.data)
table(scrna@meta.data$clonalFrequency)
scrna@meta.data$expansion <- ifelse(
  scrna@meta.data$clonalFrequency > 2, "Expanded", "Non-expanded")

table(scrna@meta.data$expansion)
table(scrna@meta.data$libraryID)

expansion_summary <- scrna@meta.data %>%
  group_by(libraryID) %>%
  summarise(
    total_cells = n(),
    expanded_cells = sum(expansion == "Expanded", na.rm = TRUE),
    non_expanded_cells = sum(expansion == "Non-expanded", na.rm = TRUE),
    expansion_percent = expanded_cells / total_cells * 100,
    .groups = "drop"
  )

expansion_summary <- expansion_summary %>%
  left_join(
    igsf3_df,
    by = "libraryID"
  )
expansion_summary <- expansion_summary %>%
  mutate(
    IGSF3_log = log1p(IGSF3))
expansion_summary <- expansion_summary %>%
  mutate(
    IGSF3_log = log1p(IGSF3)
  )

p2 = ggplot(expansion_summary, aes(x = IGSF3_log, y = expansion_percent)) +
  geom_point(size = 2.5, shape = 1) +
  geom_smooth(method = "lm", se = TRUE) +
  stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., sep = "~`  `~")),
           label.x = 5, label.y = 70, size = 4) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.title.x = element_text(family = "Arial", size = 13),
        axis.title.y = element_text(family = "Arial", size = 13),
        axis.text.x = element_text(family = "Arial", size = 11),
        axis.text.y = element_text(family = "Arial", size = 11)) +
  labs(x = "IGSF3 mean expression levels", y = "TCR expansion fraction")        
p2

### TCR diversity
clonalDiversity(
  scrna,
  cloneCall = "CTstrict",
  group.by = "libraryID",
  metric = "shannon")

tcr_diversity <- lapply(names(combined.TCR), function(id) {
  
  dat <- combined.TCR[[id]]
  
  clone_counts <- table(dat$CTstrict)            # 每个 clonotype 的细胞数
  
  clone_counts <- clone_counts[
    !names(clone_counts) %in% c("None", NA)]         # 去掉 NA / None
  
  p <- as.numeric(clone_counts) / sum(clone_counts)       # 转成频率
  shannon <- -sum(p * log(p))       # Shannon
  simpson <- 1 - sum(p^2)           # Simpson
  richness <- length(p)           # clonotype richness
  clonality <- 1 - shannon / log(richness)       # clonality
  
  data.frame(
    libraryID = id,
    TCR_cells = sum(clone_counts),
    clonotype_number = richness,
    Shannon = shannon,
    Simpson = simpson,
    Clonality = clonality
  )
  
}) %>%
  bind_rows()

cor_df <- tcr_diversity %>%
  left_join(igsf3_df, by = "libraryID")

cor_df <- cor_df %>%
  mutate(
    IGSF3_log = log1p(IGSF3))
cor.test(
  cor_df$IGSF3_log,
  cor_df$Shannon,
  method = "spearman"
)
