
rm(list = ls())

getwd()
setwd("/Users/naturexv/Datasets/PD-1/PD-1_LAG-3/GSE288199_HNSCC")

list.files("./GSE288199/")
list.files("./GSE288199/HN01_pre/")

library(Seurat)

root_dir <- "./GSE288199"
samples <- list.dirs(root_dir, recursive = FALSE, full.names = FALSE)
seu.list <- vector("list", length(samples))
names(seu.list) <- samples

for (sample in samples) {
  cat("Reading:", sample, "\n")
  counts <- Read10X(file.path(root_dir, sample))
  # RNA
  rna <- counts[["Gene Expression"]]
  # ADT
  adt <- counts[["Antibody Capture"]]
  # 创建Seurat对象
  seu <- CreateSeuratObject(
    counts = rna,
    project = sample
  )
  # 添加ADT assay
  seu[["ADT"]] <- CreateAssayObject(counts = adt)
  # metadata
  seu$sample <- sample
  seu$patient <- sub("_(pre|post)$", "", sample)
  seu$timepoint <- sub("^HN[0-9]+_", "", sample)
  seu.list[[sample]] <- seu
}

scrna <- merge(
  x = seu.list[[1]],
  y = seu.list[-1],
  add.cell.ids = names(seu.list)
)   

rm(seu.list, seu, rna, adt, counts)

mem.maxVSize(vsize = 1e9)
scrna@assays[["RNA"]] <- JoinLayers(scrna@assays[["RNA"]])


