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

  cat("\n============================\n")
  cat("Dataset:", name, "\n")
  cat("============================\n")

  ncells <- ncol(obj)

  # =========================
  # 1. preprocessing
  # =========================
  obj_s <- obj
  obj_s <- NormalizeData(obj_s, verbose = FALSE)
  obj_s <- FindVariableFeatures(obj_s, nfeatures = 2000, verbose = FALSE)
  obj_s <- ScaleData(obj_s, verbose = FALSE)
  obj_s <- RunPCA(obj_s, npcs = 30, verbose = FALSE)

  obj_s <- FindNeighbors(obj_s, dims = 1:10, verbose = FALSE)
  obj_s <- FindClusters(obj_s, resolution = 0.5, verbose = FALSE)

  # =========================
  # 2. CLUSTER-AWARE GATE (核心)
  # =========================
  cluster_tab <- table(obj_s$seurat_clusters)
  cat("Cluster distribution:\n")
  print(cluster_tab)

  # ❗ Gate 1: too few clusters
  if (length(cluster_tab) < 2) {
    cat("⚠️ Too few clusters → skip DoubletFinder, fallback singlet\n")
    obj$DoubletFinder_label <- "Singlet"
    return(obj)
  }

  # ❗ Gate 2: cluster too small
  if (min(cluster_tab) < 10) {
    cat("⚠️ Small clusters detected → lowering resolution\n")
    obj_s <- FindClusters(obj_s, resolution = 0.2, verbose = FALSE)
    cluster_tab <- table(obj_s$seurat_clusters)
  }

  # ❗ Gate 3: PCA sanity check
  pca_var <- obj_s[["pca"]]@stdev[1:10]
  if (mean(pca_var) < 0.5) {
    cat("⚠️ Weak PCA structure → skip DoubletFinder\n")
    obj$DoubletFinder_label <- "Singlet"
    return(obj)
  }

  # =========================
  # 3. SAFE paramSweep
  # =========================
  sweep_res <- tryCatch({
    paramSweep(obj_s, PCs = 1:10, sct = FALSE)
  }, error = function(e) {
    cat("⚠️ paramSweep failed → skip dataset\n")
    return(NULL)
  })

  if (is.null(sweep_res)) {
    obj$DoubletFinder_label <- "Singlet"
    return(obj)
  }

  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
  bcmvn <- find.pK(sweep_stats)

  pK_opt <- as.numeric(as.character(
    bcmvn$pK[which.max(bcmvn$BCmetric)]
  ))

  cat("Selected pK:", pK_opt, "\n")

  # =========================
  # 4. final DF run (SAFE)
  # =========================
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:10, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)

  homotypic_prop <- modelHomotypic(obj@active.ident)
  nExp <- round(ncells * dbl_rate)
  nExp_adj <- round(nExp * (1 - homotypic_prop))

  obj <- doubletFinder(
    obj,
    PCs = 1:10,
    pN = 0.25,
    pK = pK_opt,
    nExp = nExp_adj,
    reuse.pANN = NULL,
    sct = FALSE
  )

  df_col <- grep("^DF\\.classifications", colnames(obj@meta.data), value = TRUE)[1]

  if (length(df_col) == 0) {
    cat("❌ DF output missing → fallback singlet\n")
    obj$DoubletFinder_label <- "Singlet"
  } else {
    obj$DoubletFinder_label <- obj@meta.data[[df_col]]
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
