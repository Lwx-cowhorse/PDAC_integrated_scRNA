############################################################
## 🧬 05_scRNA_QC_doublet_STABLE.R
## 🎯 Publication-grade QC + Doublet detection pipeline
## 🔥 Robust version (DoubletFinder + scDblFinder fallback)
############################################################

library(Seurat)
library(harmony)
library(DoubletFinder)
library(scDblFinder)
library(SingleCellExperiment)
library(dplyr)

IN_DIR  <- "result/step2/try2.0"
OUT_DIR <- "result/step2/try2.0"
EXPECTED_DOUBLET_RATE <- 0.075

cat("===== STABLE scRNA QC PIPELINE =====\n")

# =========================
# 1. Load datasets
# =========================
objs <- list(
  CRA001160 = readRDS(file.path(IN_DIR, "CRA001160.rds")),
  GSE155698 = readRDS(file.path(IN_DIR, "GSE155698.rds")),
  GSE154778 = readRDS(file.path(IN_DIR, "GSE154778.rds")),
  GSE197177 = readRDS(file.path(IN_DIR, "GSE197177.rds"))
)

# =========================
# 2. QC function (stable)
# =========================
qc_filter <- function(obj) {
  obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- subset(obj,
                subset = nFeature_RNA > 200 &
                         nFeature_RNA < 6000 &
                         percent.mt < 15)
  return(obj)
}

objs <- lapply(objs, qc_filter)

# =========================
# 3. SAFE clustering wrapper
# =========================
safe_cluster <- function(obj, res = 0.4) {

  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)

  obj <- FindNeighbors(obj, dims = 1:15, verbose = FALSE)
  obj <- FindClusters(obj, resolution = res, verbose = FALSE)

  # ---- CLUSTER SANITY CHECK ----
  tab <- table(obj$seurat_clusters)

  if (length(tab) < 2 || min(tab) < 5) {
    cat("⚠️ Cluster unstable → lowering resolution\n")
    obj <- FindClusters(obj, resolution = 0.2, verbose = FALSE)
  }

  return(obj)
}

# =========================
# 4. Doublet detection (ROBUST)
# =========================
run_doublet <- function(obj, name, rate = 0.075) {

  cat("\n---", name, "---\n")

  ncells <- ncol(obj)
  nExp <- round(ncells * rate)

  # ========== STEP 1: try DoubletFinder ==========
  result <- tryCatch({

    obj_s <- safe_cluster(obj, res = 0.4)

    obj_s <- RunUMAP(obj_s, dims = 1:15, verbose = FALSE)

    # SAFE paramSweep
    sweep <- tryCatch({
      paramSweep(obj_s, PCs = 1:15, sct = FALSE)
    }, error = function(e) {
      cat("⚠️ paramSweep failed → fallback scDblFinder\n")
      return(NULL)
    })

    if (is.null(sweep)) stop("DF failed")

    stats <- summarizeSweep(sweep, GT = FALSE)
    bcmvn <- find.pK(stats)

    pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
    cat("pK:", pK, "\n")

    obj_s <- doubletFinder(
      obj_s,
      PCs = 1:15,
      pN = 0.25,
      pK = pK,
      nExp = nExp,
      reuse.pANN = NULL,
      sct = FALSE
    )

    df_col <- grep("DF.classifications", colnames(obj_s@meta.data), value = TRUE)[1]
    obj_s$doublet <- obj_s@meta.data[[df_col]]

    return(obj_s)

  }, error = function(e) {

    # ========== STEP 2: fallback ==========
    cat("❌ DoubletFinder failed → using scDblFinder\n")

    sce <- as.SingleCellExperiment(obj)
    sce <- scDblFinder(sce)

    obj$doublet <- sce$scDblFinder.class

    return(obj)
  })

  return(result)
}

# =========================
# 5. Run per dataset
# =========================
objs <- mapply(run_doublet,
               objs,
               names(objs),
               MoreArgs = list(rate = EXPECTED_DOUBLET_RATE),
               SIMPLIFY = FALSE)

# =========================
# 6. Remove doublets
# =========================
objs <- lapply(objs, function(x) {
  subset(x, subset = doublet == "Singlet")
})

# =========================
# 7. Integration (Harmony stable)
# =========================
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
combined <- FindClusters(combined, resolution = 0.4, verbose = FALSE)

# =========================
# 8. Save
# =========================
saveRDS(combined, file.path(OUT_DIR, "PDAC_integrated_scRNA_STABLE.rds"))

cat("\n✅ PIPELINE COMPLETE (NO CRASH VERSION)\n")
