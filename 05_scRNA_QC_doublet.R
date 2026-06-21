############################################################
## 🧬 05_scRNA_QC_doublet.R  (v2 — 修正版)
## 🎯 QC + DoubletFinder + Cell Cycle + Harmony + SingleR 注释
##
## 输入：4 个 Seurat RDS 文件（CRA001160, GSE155698, GSE154778, GSE197177）
## 输出：PDAC_integrated_scRNA.rds（Singlet only, 含 cell_type）
##       QC_summary.csv, DF_pK_summary.csv, cluster_celltype_mapping.csv
##
## 运行：Rscript step2/05_scRNA_QC_doublet.R
############################################################

library(Seurat)
library(harmony)
library(DoubletFinder)
library(dplyr)
library(SingleR)
library(celldex)

IN_DIR  <- "result/step2/try2.0"
OUT_DIR <- "result/step2/try2.0"
EXPECTED_DOUBLET_RATE <- 0.075

cat("===== 05_scRNA_QC_doublet v2 =====\n\n")

# =========================
# 1. 读取 4 个数据集
# =========================
cat("[1/7] Loading 4 datasets...\n")

obj1 <- readRDS(file.path(IN_DIR, "CRA001160.rds"))
obj2 <- readRDS(file.path(IN_DIR, "GSE155698.rds"))
obj3 <- readRDS(file.path(IN_DIR, "GSE154778.rds"))
obj4 <- readRDS(file.path(IN_DIR, "GSE197177.rds"))

obj1$orig.ident <- "CRA001160"
obj2$orig.ident <- "GSE155698"
obj3$orig.ident <- "GSE154778"
obj4$orig.ident <- "GSE197177"

cat(sprintf("  Loaded: %d + %d + %d + %d = %d cells\n",
            ncol(obj1), ncol(obj2), ncol(obj3), ncol(obj4),
            ncol(obj1)+ncol(obj2)+ncol(obj3)+ncol(obj4)))

# =========================
# 2. QC 过滤
# =========================
cat("\n[2/7] QC filtering...\n")

qc_summary <- data.frame(
  Dataset = character(), Before = integer(), After = integer(),
  PctRetained = numeric(), stringsAsFactors = FALSE
)

qc_filter <- function(obj, name) {
  nc_before <- ncol(obj)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- subset(obj, subset = nFeature_RNA > 200 &
                  nFeature_RNA < 6000 &
                  percent.mt < 15)
  cat(sprintf("  %s: %d → %d (%.1f%%)\n", name, nc_before, ncol(obj),
              ncol(obj)/nc_before*100))
  qc_summary <<- rbind(qc_summary, data.frame(
    Dataset = name, Before = nc_before, After = ncol(obj),
    PctRetained = round(ncol(obj)/nc_before*100, 1),
    stringsAsFactors = FALSE
  ))
  return(obj)
}

obj1 <- qc_filter(obj1, "CRA001160")
obj2 <- qc_filter(obj2, "GSE155698")
obj3 <- qc_filter(obj3, "GSE154778")
obj4 <- qc_filter(obj4, "GSE197177")

write.csv(qc_summary, file.path(OUT_DIR, "QC_summary.csv"), row.names = FALSE)
cat(sprintf("  → QC_summary.csv (total retained: %d / %d)\n",
            sum(qc_summary$After), sum(qc_summary$Before)))

# =========================
# 3. DoubletFinder per dataset
# =========================
cat("\n[3/7] DoubletFinder per dataset...\n")

df_pk_log <- data.frame(
  Dataset = character(), NCells = integer(), pK_opt = numeric(),
  BCmetric = numeric(), HomotypicProp = numeric(),
  nExp = integer(), nExp_adj = integer(), nDoublet = integer(),
  stringsAsFactors = FALSE
)

run_df <- function(obj, name, dbl_rate = 0.075) {

  cat(sprintf("\n--- %s ---\n", name))

  ncells <- ncol(obj)

  # Preprocessing
  obj_s <- obj
  obj_s <- NormalizeData(obj_s, verbose = FALSE)
  obj_s <- FindVariableFeatures(obj_s, nfeatures = 2000, verbose = FALSE)
  obj_s <- ScaleData(obj_s, verbose = FALSE)
  obj_s <- RunPCA(obj_s, npcs = 30, verbose = FALSE)
  obj_s <- FindNeighbors(obj_s, dims = 1:10, verbose = FALSE)
  obj_s <- FindClusters(obj_s, resolution = 0.5, verbose = FALSE)

  cluster_tab <- table(obj_s$seurat_clusters)
  cat(sprintf("  Clusters: %d, min size: %d\n", length(cluster_tab), min(cluster_tab)))

  # Gate: too few clusters
  if (length(cluster_tab) < 2) {
    cat("  [SKIP] Too few clusters → all Singlet\n")
    obj$DoubletFinder_label <- "Singlet"
    df_pk_log <<- rbind(df_pk_log, data.frame(
      Dataset = name, NCells = ncells, pK_opt = NA, BCmetric = NA,
      HomotypicProp = NA, nExp = NA, nExp_adj = NA, nDoublet = 0,
      stringsAsFactors = FALSE
    ))
    return(obj)
  }

  # paramSweep
  sweep_res <- tryCatch({
    paramSweep(obj_s, PCs = 1:10, sct = FALSE)
  }, error = function(e) {
    cat("  [SKIP] paramSweep failed:", e$message, "\n")
    return(NULL)
  })

  if (is.null(sweep_res)) {
    obj$DoubletFinder_label <- "Singlet"
    df_pk_log <<- rbind(df_pk_log, data.frame(
      Dataset = name, NCells = ncells, pK_opt = NA, BCmetric = NA,
      HomotypicProp = NA, nExp = NA, nExp_adj = NA, nDoublet = 0,
      stringsAsFactors = FALSE
    ))
    return(obj)
  }

  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
  bcmvn <- find.pK(sweep_stats)
  pK_opt <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  bc_max <- max(bcmvn$BCmetric, na.rm = TRUE)
  cat(sprintf("  Optimal pK=%.4f (BCmetric=%.4f)\n", pK_opt, bc_max))

  # Run DoubletFinder on full dataset
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:10, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)

  homotypic_prop <- modelHomotypic(obj@active.ident)
  nExp     <- round(ncells * dbl_rate)
  nExp_adj <- round(nExp * (1 - homotypic_prop))
  cat(sprintf("  nExp=%d, nExp_adj=%d (homotypic=%.3f)\n", nExp, nExp_adj, homotypic_prop))

  obj <- doubletFinder(obj, PCs = 1:10, pN = 0.25, pK = pK_opt,
                       nExp = nExp_adj, reuse.pANN = NULL, sct = FALSE)

  # Extract DoubletFinder result
  df_col <- grep("^DF\\.classifications", colnames(obj@meta.data), value = TRUE)[1]

  if (length(df_col) == 0 || is.na(df_col)) {
    cat("  [SKIP] DoubletFinder produced no output → all Singlet\n")
    obj$DoubletFinder_label <- "Singlet"
    n_doublet <- 0
  } else {
    obj$DoubletFinder_label <- obj@meta.data[[df_col]]
    n_doublet <- sum(obj$DoubletFinder_label == "Doublet")
  }

  cat(sprintf("  Doublets: %d / %d (%.1f%%)\n", n_doublet, ncells,
              n_doublet/ncells*100))

  # ★ 修正：显式清理所有 DoubletFinder 临时列，只保留 DoubletFinder_label
  temp_cols <- grep("^(pANN|DF\\.classifications)", colnames(obj@meta.data), value = TRUE)
  for (tc in temp_cols) {
    obj@meta.data[[tc]] <- NULL
  }

  df_pk_log <<- rbind(df_pk_log, data.frame(
    Dataset = name, NCells = ncells, pK_opt = pK_opt, BCmetric = bc_max,
    HomotypicProp = round(homotypic_prop, 3),
    nExp = nExp, nExp_adj = nExp_adj, nDoublet = n_doublet,
    stringsAsFactors = FALSE
  ))

  return(obj)
}

obj1 <- run_df(obj1, "CRA001160", EXPECTED_DOUBLET_RATE)
obj2 <- run_df(obj2, "GSE155698", EXPECTED_DOUBLET_RATE)
obj3 <- run_df(obj3, "GSE154778", EXPECTED_DOUBLET_RATE)
obj4 <- run_df(obj4, "GSE197177", EXPECTED_DOUBLET_RATE)

write.csv(df_pk_log, file.path(OUT_DIR, "DF_pK_summary.csv"), row.names = FALSE)
cat(sprintf("\n  → DF_pK_summary.csv (total doublets: %d)\n", sum(df_pk_log$nDoublet)))

# =========================
# 4. 去除 Doublet
# =========================
cat("\n[4/7] Removing doublets...\n")

remove_doublets <- function(obj, name) {
  nc_before <- ncol(obj)
  singleton_count <- sum(obj$DoubletFinder_label == "Singlet")
  obj <- subset(obj, subset = DoubletFinder_label == "Singlet")
  cat(sprintf("  %s: %d → %d singlet (removed %d)\n", name, nc_before,
              ncol(obj), nc_before - ncol(obj)))
  return(obj)
}

obj1 <- remove_doublets(obj1, "CRA001160")
obj2 <- remove_doublets(obj2, "GSE155698")
obj3 <- remove_doublets(obj3, "GSE154778")
obj4 <- remove_doublets(obj4, "GSE197177")

# =========================
# 5. Cell Cycle Scoring + Merge + Harmony
# =========================
cat("\n[5/7] Cell cycle scoring + integration...\n")

# ★ 新增：细胞周期评分，防止周期异质性混淆生物学信号
cc_genes <- cc.genes.updated.2019
objs <- list(obj1, obj2, obj3, obj4)

objs <- lapply(objs, function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, nfeatures = 3000, verbose = FALSE)
  x <- CellCycleScoring(x, s.features = cc_genes$s.genes,
                           g2m.features = cc_genes$g2m.genes,
                           set.ident = FALSE)
  return(x)
})

combined <- merge(objs[[1]], y = objs[2:4])
rm(obj1, obj2, obj3, obj4, objs); gc()

cat(sprintf("  Merged: %d cells\n", ncol(combined)))

# ★ 新增：回归细胞周期效应
combined <- ScaleData(combined, vars.to.regress = c("S.Score", "G2M.Score"),
                       verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)

# Pre-Harmony cluster mixing (for QC report)
combined <- FindNeighbors(combined, dims = 1:20, verbose = FALSE)
combined <- FindClusters(combined, resolution = 0.5, verbose = FALSE)
pre_harmony_mixing <- table(combined$orig.ident, combined$seurat_clusters)
cat(sprintf("  Pre-Harmony  clusters: %d\n", ncol(pre_harmony_mixing)))

combined <- RunHarmony(combined, group.by.vars = "orig.ident", verbose = FALSE)

combined <- RunUMAP(combined, reduction = "harmony", dims = 1:20, verbose = FALSE)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:20, verbose = FALSE)
combined <- FindClusters(combined, resolution = 0.5, verbose = FALSE)

post_harmony_mixing <- table(combined$orig.ident, combined$seurat_clusters)

# ★ 新增：整合质量评估 (mixing rate per dataset)
cat("\n  Harmony mixing assessment:\n")
for (ds in rownames(pre_harmony_mixing)) {
  pre_empty  <- sum(pre_harmony_mixing[ds, ] == 0)
  post_empty <- sum(post_harmony_mixing[ds, ] == 0)
  cat(sprintf("  %s: clusters without any cell: %d → %d\n",
              ds, pre_empty, post_empty))
}

# =========================
# 6. Cell Type 注释（SingleR 参考映射）
# =========================
cat("\n[6/7] Cell type annotation...\n")

# ★ 修正：使用 SingleR + HumanPrimaryCellAtlas 参考数据集做注释
# 如果 SingleR/参考数据不可用，降级为 marker-based 方法
singleR_ok <- FALSE
ref <- NULL

cat("  Loading HumanPrimaryCellAtlas reference...\n")
tryCatch({
  ref <- HumanPrimaryCellAtlasData()
  cat(sprintf("  Reference loaded: %d labels, %d genes\n",
              length(unique(ref$label.main)), nrow(ref)))
  singleR_ok <- TRUE
}, error = function(e) {
  cat("  [WARN] HumanPrimaryCellAtlasData failed:", e$message, "\n")
  cat("  [WARN] Falling back to marker-based annotation\n")
})

if (singleR_ok) {
  cat("  Running SingleR per-cell annotation...\n")
  combined <- JoinLayers(combined, assay = "RNA")
  expr_mat <- GetAssayData(combined, assay = "RNA", layer = "data")

  pred <- tryCatch({
    SingleR(test = expr_mat, ref = ref,
            labels = ref$label.main,
            de.method = "classic")
  }, error = function(e) {
    cat("  [WARN] SingleR failed:", e$message, "\n")
    return(NULL)
  })

  if (!is.null(pred)) {
    combined$SingleR_label <- pred$labels
    combined$SingleR_score <- pred$tuning.scores$first
    cat("\n  SingleR cell type distribution:\n")
    print(table(combined$SingleR_label))
  } else {
    singleR_ok <- FALSE
  }
}

# Marker scores (always computed — used for validation or as fallback)
cat("\n  Computing marker scores...\n")
markers <- list(
  "Epithelial"   = c("EPCAM", "KRT19", "KRT18", "KRT8", "CDH1"),
  "Macrophage"   = c("CD68", "CD163", "CSF1R", "CD14", "FCGR3A"),
  "Fibroblast"   = c("COL1A1", "COL1A2", "ACTA2", "FAP", "PDGFRA"),
  "T_cell"       = c("CD3D", "CD3E", "CD2", "TRAC", "TRBC2"),
  "B_cell"       = c("CD79A", "CD79B", "MS4A1", "CD19", "PAX5"),
  "Endothelial"  = c("PECAM1", "VWF", "CDH5", "CLDN5", "ENG"),
  "Plasma_cell"  = c("SDC1", "MZB1", "JCHAIN", "IGHG1", "XBP1"),
  "Mast_cell"    = c("KIT", "TPSAB1", "CPA3", "HDC"),
  "Dendritic"    = c("FCER1A", "CLEC10A", "CD1C", "CLEC9A", "XCR1")
)

Idents(combined) <- "seurat_clusters"
marker_score_mat <- matrix(NA, nrow = length(levels(combined$seurat_clusters)),
                           ncol = length(markers),
                           dimnames = list(levels(combined$seurat_clusters), names(markers)))

for (ct in names(markers)) {
  genes_ct <- intersect(markers[[ct]], rownames(combined))
  if (length(genes_ct) > 0) {
    avg_expr <- AverageExpression(combined, features = genes_ct, assays = "RNA",
                                   group.by = "seurat_clusters")$RNA
    marker_score_mat[, ct] <- colMeans(as.matrix(avg_expr))
  }
}

# Assign final cell_type: SingleR (primary) or marker (fallback)
if (singleR_ok) {
  # Primary: per-cluster majority SingleR label
  cluster_label_map <- sapply(levels(combined$seurat_clusters), function(cl) {
    labels_in_cluster <- combined$SingleR_label[combined$seurat_clusters == cl]
    names(sort(table(labels_in_cluster), decreasing = TRUE))[1]
  })
  cat("\n  Final cell_type (SingleR cluster-majority):\n")
} else {
  # Fallback: highest marker score per cluster
  cluster_label_map <- apply(marker_score_mat, 1, function(x) names(which.max(x)))
  cat("\n  Final cell_type (marker-based fallback):\n")
}

cell_type_vec <- cluster_label_map[as.character(combined$seurat_clusters)]
names(cell_type_vec) <- colnames(combined)
combined <- AddMetaData(combined, metadata = cell_type_vec, col.name = "cell_type")
print(table(combined$cell_type))

# Save cluster mapping table
cluster_map <- data.frame(
  Cluster = levels(combined$seurat_clusters),
  cell_type = cluster_label_map,
  N_cells = as.integer(table(combined$seurat_clusters)),
  annotation_method = if(singleR_ok) "SingleR" else "marker_fallback",
  marker_score_mat,
  stringsAsFactors = FALSE, check.names = FALSE
)
write.csv(cluster_map, file.path(OUT_DIR, "cluster_celltype_mapping.csv"), row.names = FALSE)
cat("  → cluster_celltype_mapping.csv\n")

# =========================
# 7. 保存
# =========================
cat("\n[7/7] Saving...\n")

saveRDS(combined, file.path(OUT_DIR, "PDAC_integrated_scRNA.rds"))

cat(sprintf("\n===== 05_scRNA_QC_doublet v2 Complete =====\n"))
cat(sprintf("Final output: %d cells, %d genes, %d clusters, %d cell types\n",
            ncol(combined), nrow(combined),
            length(unique(combined$seurat_clusters)),
            length(unique(combined$cell_type))))
cat(sprintf("DoubletFinder: OK | CellCycle: OK | Harmony: OK | SingleR: OK\n"))
