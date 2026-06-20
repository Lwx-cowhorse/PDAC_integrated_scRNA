############################################################
## 🧬 05_scRNA_QC_doublet.R
## 🎯 QC + DoubletFinder + Harmony 整合（Linux 64GB 完整版）
############################################################

library(Seurat)
library(harmony)
library(DoubletFinder)
library(dplyr)

# ===== scDblFinder（已改为可选加载，避免崩溃）=====
if (requireNamespace("scDblFinder", quietly = TRUE)) {
  library(scDblFinder)
}

IN_DIR  <- "result/step2/try2.0"
OUT_DIR <- "result/step2/try2.0"
EXPECTED_DOUBLET_RATE <- 0.075

cat("===== 05_scRNA_QC_doublet =====\n\n")

# =========================
# 1. 读取 4 个数据集
# =========================
obj1 <- readRDS(file.path(IN_DIR, "CRA001160.rds"))
obj2 <- readRDS(file.path(IN_DIR, "GSE155698.rds"))
obj3 <- readRDS(file.path(IN_DIR, "GSE154778.rds"))
obj4 <- readRDS(file.path(IN_DIR, "GSE197177.rds"))

obj1$orig.ident <- "CRA001160"
obj2$orig.ident <- "GSE155698"
obj3$orig.ident <- "GSE154778"
obj4$orig.ident <- "GSE197177"

# =========================
# 2. QC 过滤
# =========================
qc_filter <- function(obj, name) {
  nc_before <- ncol(obj)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- subset(obj, subset = nFeature_RNA > 200 &
                  nFeature_RNA < 6000 &
                  percent.mt < 15)
  cat(sprintf("  %s: %d → %d\n", name, nc_before, ncol(obj)))
  return(obj)
}

obj1 <- qc_filter(obj1, "CRA001160")
obj2 <- qc_filter(obj2, "GSE155698")
obj3 <- qc_filter(obj3, "GSE154778")
obj4 <- qc_filter(obj4, "GSE197177")

# =========================
# 3. DoubletFinder
# =========================
run_df <- function(obj, name, dbl_rate = 0.075) {

  ncells <- ncol(obj)
  cat(sprintf("\n%s: %d cells\n", name, ncells))

  obj_sweep <- obj
  obj_sweep <- NormalizeData(obj_sweep, verbose = FALSE)
  obj_sweep <- FindVariableFeatures(obj_sweep, nfeatures = 2000, verbose = FALSE)
  obj_sweep <- ScaleData(obj_sweep, verbose = FALSE)
  obj_sweep <- RunPCA(obj_sweep, npcs = 30, verbose = FALSE)
  obj_sweep <- FindNeighbors(obj_sweep, dims = 1:10, verbose = FALSE)
  obj_sweep <- FindClusters(obj_sweep, resolution = 0.5, verbose = FALSE)

  sweep_res   <- paramSweep(obj_sweep, PCs = 1:10, sct = FALSE)
  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
  bcmvn       <- find.pK(sweep_stats)
  pK_opt <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:10, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)

  homotypic_prop <- modelHomotypic(obj@active.ident)
  nExp <- round(ncells * dbl_rate)
  nExp_adj <- round(nExp * (1 - homotypic_prop))

  obj <- doubletFinder(obj,
                       PCs = 1:10,
                       pN = 0.25,
                       pK = pK_opt,
                       nExp = nExp_adj,
                       reuse.pANN = NULL,
                       sct = FALSE)

  df_col <- grep("^DF\\.classifications", colnames(obj@meta.data), value = TRUE)[1]
  obj$DoubletFinder_label <- obj@meta.data[[df_col]]

  # ===== scDblFinder fallback（已改为安全版本）=====
  if (requireNamespace("scDblFinder", quietly = TRUE)) {
    sce <- as.SingleCellExperiment(obj)
    sce <- scDblFinder(sce)
    obj$scDblFinder_label <- sce$scDblFinder.class
  } else {
    obj$scDblFinder_label <- "Not_Installed"
  }

  return(obj)
}

obj1 <- run_df(obj1, "CRA001160", EXPECTED_DOUBLET_RATE)
obj2 <- run_df(obj2, "GSE155698", EXPECTED_DOUBLET_RATE)
obj3 <- run_df(obj3, "GSE154778", EXPECTED_DOUBLET_RATE)
obj4 <- run_df(obj4, "GSE197177", EXPECTED_DOUBLET_RATE)

# =========================
# 4. 去除 Doublet
# =========================
remove_doublets <- function(obj, name) {
  nc_before <- ncol(obj)
  obj <- subset(obj, subset = DoubletFinder_label == "Singlet")
  cat(sprintf("%s: %d → %d\n", name, nc_before, ncol(obj)))
  return(obj)
}

obj1 <- remove_doublets(obj1, "CRA001160")
obj2 <- remove_doublets(obj2, "GSE155698")
obj3 <- remove_doublets(obj3, "GSE154778")
obj4 <- remove_doublets(obj4, "GSE197177")

# =========================
# 5. Harmony整合
# =========================
objs <- list(obj1, obj2, obj3, obj4)

objs <- lapply(objs, function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, nfeatures = 3000, verbose = FALSE)
  return(x)
})

combined <- merge(objs[[1]], y = objs[2:4])

combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunHarmony(combined, group.by.vars = "orig.ident", verbose = FALSE)

combined <- RunUMAP(combined, reduction = "harmony", dims = 1:20, verbose = FALSE)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:20, verbose = FALSE)
combined <- FindClusters(combined, resolution = 0.5, verbose = FALSE)

saveRDS(combined, file.path(OUT_DIR, "PDAC_integrated_scRNA.rds"))

cat("\nDONE\n")
