############################################################
## 🧬 05_scRNA_QC_doublet.R
## 🎯 QC + DoubletFinder + Harmony 整合（Linux 64GB 完整版）
##
## ⚠️  此脚本必须在 Linux 上运行（GitHub Actions / 服务器）
##     Windows 上 DoubletFinder 不可用
##
## 输入：4 个 Seurat RDS 文件（由 00c 生成）
##      CRA001160.rds, GSE155698.rds, GSE154778.rds, GSE197177.rds
##
## 输出：PDAC_integrated_scRNA.rds（含 Doublet 去除 + cell_type 注释）
##
## 运行方式：
##   Rscript step2/05_scRNA_QC_doublet.R
##
## GitHub Actions 附件放置：
##   将 4 个 RDS 文件放入 result/step2/try2.0/（或修改下方 IN_DIR）
############################################################

library(Seurat)
library(harmony)
library(DoubletFinder)
library(dplyr)

# ╔══════════════════════════════════════════════════════════╗
# ║  ⚙️  配置                                               ║
# ╚══════════════════════════════════════════════════════════╝
IN_DIR  <- "result/step2/try2.0"      # 输入目录（4个RDS所在）
OUT_DIR <- "result/step2/try2.0"      # 输出目录
EXPECTED_DOUBLET_RATE <- 0.075        # 预期 doublet 率 (10x 标准: ~0.8% per 1K cells)

# =========================
# 1. 读取 4 个数据集
# =========================
cat("===== 05_scRNA_QC_doublet =====\n\n")
cat("[1/6] Loading 4 datasets...\n")

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
cat("\n[2/6] QC filtering per dataset...\n")

qc_filter <- function(obj, name) {
  nc_before <- ncol(obj)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- subset(obj, subset = nFeature_RNA > 200 &
                         nFeature_RNA < 6000 &
                         percent.mt < 15)
  cat(sprintf("  %s: %d → %d (%.1f%% retained)\n",
              name, nc_before, ncol(obj), ncol(obj)/nc_before*100))
  return(obj)
}

obj1 <- qc_filter(obj1, "CRA001160")
obj2 <- qc_filter(obj2, "GSE155698")
obj3 <- qc_filter(obj3, "GSE154778")
obj4 <- qc_filter(obj4, "GSE197177")

# =========================
# 3. DoubletFinder per dataset
# =========================
cat("\n[3/6] DoubletFinder per dataset...\n")

run_df <- function(obj, name, dbl_rate = 0.075) {
  ncells <- ncol(obj)
  cat(sprintf("  %s: %d cells, expected doublet ≈ %.1f%%\n",
              name, ncells, dbl_rate*100))

  # For >30K cells, subsample for pK optimization
  use_subsample <- ncells > 30000
  obj_sweep <- obj
  if (use_subsample) {
    set.seed(123)
    obj_sweep <- subset(obj, cells = sample(colnames(obj), 30000))
    cat(sprintf("    Using 30K subsample for pK sweep\n"))
  }

  # Pre-processing for DoubletFinder
  obj_sweep <- NormalizeData(obj_sweep, verbose = FALSE)
  obj_sweep <- FindVariableFeatures(obj_sweep, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  obj_sweep <- ScaleData(obj_sweep, verbose = FALSE)
  obj_sweep <- RunPCA(obj_sweep, npcs = 30, verbose = FALSE)
  obj_sweep <- RunUMAP(obj_sweep, dims = 1:10, verbose = FALSE)
  obj_sweep <- FindNeighbors(obj_sweep, dims = 1:10, verbose = FALSE)
  obj_sweep <- FindClusters(obj_sweep, resolution = 0.5, verbose = FALSE)

  # pK optimization
  sweep_res   <- paramSweep(obj_sweep, PCs = 1:10, sct = FALSE)
  sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
  bcmvn       <- find.pK(sweep_stats)
  pK_opt <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  cat(sprintf("    Optimal pK = %.4f\n", pK_opt))

  # Pre-process full dataset
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:10, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:10, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)

  # Homotypic correction
  homotypic_prop <- modelHomotypic(obj@active.ident)
  nExp_poi     <- round(dbl_rate * ncells)
  nExp_poi_adj <- round(nExp_poi * (1 - homotypic_prop))
  cat(sprintf("    nExp=%.0f, nExp_adj=%.0f\n", nExp_poi, nExp_poi_adj))

  # Run DoubletFinder
  obj <- doubletFinder_v3(obj, PCs = 1:10, pN = 0.25, pK = pK_opt,
                          nExp = nExp_poi_adj, reuse.pANN = FALSE, sct = FALSE)

  # Extract result
  df_col <- grep("^DF\\.classifications", colnames(obj@meta.data), value = TRUE)[1]
  obj$DoubletFinder_label <- obj@meta.data[[df_col]]

  n_doublet <- sum(obj$DoubletFinder_label == "Doublet")
  cat(sprintf("    Doublets detected: %d / %d (%.1f%%)\n",
              n_doublet, ncells, n_doublet/ncells*100))

  # Clean temp columns
  for (tc in grep("^(DF\\.|pANN)", colnames(obj@meta.data), value = TRUE)) {
    obj@meta.data[[tc]] <- NULL
  }
  rm(obj_sweep); gc()
  return(obj)
}

obj1 <- run_df(obj1, "CRA001160")
obj2 <- run_df(obj2, "GSE155698")
obj3 <- run_df(obj3, "GSE154778")
obj4 <- run_df(obj4, "GSE197177")

# =========================
# 4. 去除 Doublet 细胞
# =========================
cat("\n[4/6] Removing doublet cells...\n")

remove_doublets <- function(obj, name) {
  nc_before <- ncol(obj)
  obj <- subset(obj, subset = DoubletFinder_label == "Singlet")
  cat(sprintf("  %s: %d → %d (removed %d doublets)\n",
              name, nc_before, ncol(obj), nc_before - ncol(obj)))
  return(obj)
}

obj1 <- remove_doublets(obj1, "CRA001160")
obj2 <- remove_doublets(obj2, "GSE155698")
obj3 <- remove_doublets(obj3, "GSE154778")
obj4 <- remove_doublets(obj4, "GSE197177")

# =========================
# 5. 整合：Normalize + HVG → Merge → Harmony → UMAP → Clustering
# =========================
cat("\n[5/6] Integration (Harmony)...\n")

objs <- list(obj1, obj2, obj3, obj4)
objs <- lapply(objs, function(x) {
  x <- JoinLayers(x, assay = "RNA")
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  return(x)
})

combined <- merge(objs[[1]], y = objs[2:4],
                  add.cell.ids = c("CRA", "GSE155", "GSE154", "GSE197"))
rm(obj1, obj2, obj3, obj4, objs); gc()

cat(sprintf("  Merged: %d cells\n", ncol(combined)))

combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 50, verbose = FALSE)
combined <- RunHarmony(combined, group.by.vars = "orig.ident", verbose = FALSE)
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30, verbose = FALSE)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30, verbose = FALSE)
combined <- FindClusters(combined, resolution = 0.5, verbose = FALSE)

# =========================
# 6. Cell Type 注释 + 保存
# =========================
cat("\n[6/6] Cell type annotation...\n")

markers <- list(
  "Epithelial"   = c("EPCAM", "KRT19", "KRT18", "KRT8", "CDH1"),
  "TAM"          = c("CD68", "CD163", "CSF1R", "CD14", "FCGR3A"),
  "CAF"          = c("COL1A1", "COL1A2", "ACTA2", "FAP", "PDGFRA"),
  "T_cell"       = c("CD3D", "CD3E", "CD2", "TRAC", "TRBC2"),
  "B_cell"       = c("CD79A", "CD79B", "MS4A1", "CD19", "PAX5"),
  "Endothelial"  = c("PECAM1", "VWF", "CDH5", "CLDN5", "ENG"),
  "Plasma_cell"  = c("SDC1", "MZB1", "JCHAIN", "IGHG1", "XBP1"),
  "Mast_cell"    = c("KIT", "TPSAB1", "CPA3", "HDC"),
  "Dendritic"    = c("FCER1A", "CLEC10A", "CD1C", "CLEC9A", "XCR1")
)

Idents(combined) <- "seurat_clusters"
cluster_ids <- levels(combined$seurat_clusters)
score_matrix <- matrix(NA, nrow = length(cluster_ids), ncol = length(markers),
                       dimnames = list(cluster_ids, names(markers)))

for (ct in names(markers)) {
  genes_ct <- intersect(markers[[ct]], rownames(combined))
  if (length(genes_ct) > 0) {
    avg_expr <- AverageExpression(combined, features = genes_ct, assays = "RNA",
                                  group.by = "seurat_clusters")$RNA
    score_matrix[, ct] <- colMeans(as.matrix(avg_expr))
  }
}

cluster_to_celltype <- apply(score_matrix, 1, function(x) names(which.max(x)))
combined@meta.data$cell_type <- cluster_to_celltype[as.character(combined$seurat_clusters)]

cat("  Cell type distribution:\n")
print(table(combined$cell_type))

cat(sprintf("  Saving to %s/PDAC_integrated_scRNA.rds ...\n", OUT_DIR))
saveRDS(combined, file.path(OUT_DIR, "PDAC_integrated_scRNA.rds"))
write.csv(score_matrix, file.path(OUT_DIR, "cluster_celltype_marker_scores.csv"), row.names = FALSE)

cat(sprintf("\n=== 05_scRNA_QC_doublet Complete ===\n"))
cat(sprintf("Output: %d cells, %d clusters, %d cell types\n",
            ncol(combined), length(unique(combined$seurat_clusters)),
            length(unique(combined$cell_type))))
cat(sprintf("DoubletFinder: ✅ Run | Harmony: ✅ | Cell types: ✅\n"))
