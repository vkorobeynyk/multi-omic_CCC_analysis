# =============================================================================
# helper_functions.R
#
# Shared functions for 03_RCTD_SpotDeconvolution.qmd.
# Source this at the top of the .qmd with: source("helper_functions.R")
#
# Sections:
#   1. Style constants (colors)
#   2. RCTD reference building / running
#   3. Composition-aware colocalization (CLR transform, permutation z-score,
#      partial correlation across all cell types in a region)
#   4. Marker validation helpers
#   5. Ligand-receptor vs colocalization: core statistics
#   6. Ligand-receptor vs colocalization: one-call-per-pair pipeline
#   7. Plotting helpers (generic, reused across cell-type pairs)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(forcats)
  library(ggrepel)
  library(ppcor)
})

# -----------------------------------------------------------------------------
# 1. Style constants (Okabe-Ito colorblind-safe palette, used everywhere)
# -----------------------------------------------------------------------------

AGE_COLORS <- c(young = "#0072B2", middle = "#009E73", old = "#D55E00")
SIG_COLOR  <- "#D55E00"
SIG_DIRECTION_COLORS <- c(Positive = "#D55E00", Negative = "#0072B2")

# -----------------------------------------------------------------------------
# 2. RCTD reference building / running
# -----------------------------------------------------------------------------

#' Build a spacexr::Reference restricted to one region's candidate cell types
build_region_reference <- function(celltypes_keep, ref_counts, cell_types) {
  keep_cells <- names(cell_types)[cell_types %in% celltypes_keep]

  counts_sub <- ref_counts[, keep_cells]
  types_sub  <- droplevels(cell_types[keep_cells])

  Reference(counts_sub, types_sub, n_max_cells = 5000)
}

#' Run RCTD (full/doublet mode) on one sample x region subset of a Seurat object
run_rctd_region <- function(sp_obj_sub, reference, max_cores = 4) {
  sp_counts <- GetAssayData(sp_obj_sub, assay = "Spatial", layer = "counts")

  coords <- GetTissueCoordinates(sp_obj_sub)
  coords <- coords[, c("x", "y")]
  colnames(coords) <- c("x", "y")
  rownames(coords) <- colnames(sp_counts)

  nUMI_sp <- colSums(sp_counts)
  names(nUMI_sp) <- colnames(sp_counts)

  puck <- SpatialRNA(coords, sp_counts, nUMI_sp)

  myRCTD <- create.RCTD(puck, reference, max_cores = max_cores)
  run.RCTD(myRCTD, doublet_mode = "full")
}

#' Export one sample's merged (across regions) coords/expression/weights for LIANA+ Moran's I
#'
#' Pulls every RCTD run belonging to one physical slide, takes the union of
#' spots + real coordinates across regions (preserves true region-border
#' adjacency), and writes matrix/features/barcodes/coords/weights files.
export_rctd_slice <- function(sample_name, rctd_results, merged_obj, export_dir) {
  runs <- rctd_results[sapply(rctd_results, function(r) r$sample == sample_name)]

  coords_list <- lapply(runs, function(r) r$rctd_obj@spatialRNA@coords)
  names(coords_list) <- NULL   # prevents do.call(rbind, ...) from prefixing rownames with list names

  coords_all <- do.call(rbind, coords_list)
  coords_all <- coords_all[!duplicated(rownames(coords_all)), ]
  barcodes_all <- rownames(coords_all)

  norm_expr <- GetAssayData(merged_obj, assay = "SCT", layer = "data")[, barcodes_all]

  # weights from each region-run, merged by barcode — DG-only columns (GC, Neurogenic_lineage)
  # will be NA for ML spots and vice versa; that's expected, since those cell types were
  # never candidates in the other region
  weights_list <- lapply(runs, function(r) {
    w <- as.data.frame(as.matrix(normalize_weights(r$rctd_obj@results$weights)))
    w$barcode <- rownames(w)
    w
  })
  weights_sample <- Reduce(function(x, y) merge(x, y, by = "barcode", all = TRUE), weights_list)

  region_df <- do.call(rbind, lapply(runs, function(r) {
    data.frame(barcode = rownames(r$rctd_obj@spatialRNA@coords), region = r$region)
  }))
  weights_sample <- merge(weights_sample, region_df, by = "barcode", all.x = TRUE)

  slice_dir <- file.path(export_dir, sample_name)
  dir.create(slice_dir, showWarnings = FALSE)

  Matrix::writeMM(as(norm_expr, "sparseMatrix"), file.path(slice_dir, "matrix.mtx"))
  writeLines(rownames(norm_expr), file.path(slice_dir, "features.tsv"))
  writeLines(colnames(norm_expr), file.path(slice_dir, "barcodes.tsv"))

  write.csv(data.frame(barcode = rownames(coords_all), x = coords_all$x, y = coords_all$y),
            file.path(slice_dir, "coords.csv"), row.names = FALSE)
  write.csv(weights_sample, file.path(slice_dir, "weights.csv"), row.names = FALSE)

  invisible(slice_dir)
}

# -----------------------------------------------------------------------------
# 3. Composition-aware colocalization across ALL cell types in a region
# -----------------------------------------------------------------------------

#' Centered log-ratio transform: subtract each row's mean log = divide by
#' that row's geometric mean. Makes standard correlation valid for RCTD
#' weights, which sum to 1 per spot and are otherwise not independent.
clr_transform <- function(mat, pseudo = 1e-4) {
  mat <- as.matrix(mat) + pseudo
  log_mat <- log(mat)
  sweep(log_mat, 1, rowMeans(log_mat), "-")
}

#' Empirical null for CLR + correlation: shuffle each cell type's weights
#' across spots independently (breaking real pairing, preserving that
#' type's own distribution), recompute CLR + correlation, repeat n_perm times.
#' Seeds internally (rather than relying on a set.seed() call elsewhere in the
#' script) so this is reproducible regardless of what ran before it or what
#' order/subset of chunks got (re-)executed.
permute_null_clr_cor <- function(mat, n_perm = 500, pseudo = 1e-4, method = "pearson", seed = 1) {
  set.seed(seed)
  ct <- colnames(mat); K <- length(ct); n <- nrow(mat)
  null_cor <- array(NA, dim = c(K, K, n_perm), dimnames = list(ct, ct, NULL))

  for (p in seq_len(n_perm)) {
    perm_mat <- mat
    for (j in seq_len(K)) {
      perm_mat[, j] <- mat[sample(n), j]
    }
    null_cor[, , p] <- cor(clr_transform(perm_mat, pseudo = pseudo), method = method)
  }
  null_cor
}

#' Full colocalization analysis for one region: observed CLR correlation,
#' permutation null, z-score, and empirical two-sided p-value (BH-adjusted).
analyze_coloc_for_region <- function(weights_df, region_name, region_celltypes_list,
                                      n_perm = 500, pseudo = 1e-4, seed = 1) {
  ct_cols <- region_celltypes_list[[region_name]]
  sub     <- weights_df %>% filter(region == region_name)
  raw_mat <- as.matrix(sub[, ct_cols])

  stopifnot(!any(is.na(raw_mat)))   # only this region's own candidates should be present

  clr_mat <- clr_transform(raw_mat, pseudo = pseudo)
  obs_cor <- cor(clr_mat, method = "pearson")

  null_dist <- permute_null_clr_cor(raw_mat, n_perm = n_perm, pseudo = pseudo, seed = seed)
  null_mean <- apply(null_dist, c(1, 2), mean)
  null_sd   <- apply(null_dist, c(1, 2), sd)
  z_mat     <- (obs_cor - null_mean) / null_sd

  K <- length(ct_cols)
  p_vals <- matrix(NA, K, K, dimnames = list(ct_cols, ct_cols))
  for (i in seq_len(K)) for (j in seq_len(K)) {
    if (i >= j) next
    nv <- null_dist[i, j, ]; ob <- obs_cor[i, j]
    p_vals[i, j] <- p_vals[j, i] <- mean(abs(nv - mean(nv)) >= abs(ob - mean(nv)))
  }
  p_adj <- matrix(p.adjust(p_vals, method = "BH"), nrow = K, dimnames = dimnames(p_vals))

  list(obs_cor = obs_cor, z_mat = z_mat, p_adj = p_adj, n_spots = nrow(sub))
}

#' Long-format z-score table from analyze_coloc_for_region(), diagonal dropped
make_z_df <- function(coloc_res, region_name) {
  as.data.frame(as.table(coloc_res$z_mat)) %>%
    rename(ct1 = Var1, ct2 = Var2, z = Freq) %>%
    mutate(region = region_name) %>%
    filter(ct1 != ct2)
}

#' Partial correlation (CLR-transformed weights) between every pair of cell
#' types in a region, controlling for all other candidate cell types.
compute_pcor_region <- function(weights_df, region_name, region_celltypes_list) {
  ct_cols <- region_celltypes_list[[region_name]]

  sub     <- weights_df %>% dplyr::filter(region == region_name)
  raw_mat <- as.matrix(sub[, ct_cols])
  clr_mat <- clr_transform(raw_mat)

  pc <- pcor(as.data.frame(clr_mat), method = "pearson")
  colnames(pc$estimate) <- rownames(pc$estimate) <- colnames(clr_mat)
  colnames(pc$p.value)  <- rownames(pc$p.value)  <- colnames(clr_mat)
  pc
}

#' Long-format partial-correlation table with BH-adjusted p and significance stars
make_pcor_df <- function(pc, region_name) {
  est_df <- as.data.frame(as.table(pc$estimate)) %>% rename(ct1 = Var1, ct2 = Var2, estimate = Freq)
  p_df   <- as.data.frame(as.table(pc$p.value))  %>% rename(ct1 = Var1, ct2 = Var2, p = Freq)

  est_df %>%
    left_join(p_df, by = c("ct1", "ct2")) %>%
    filter(ct1 != ct2) %>%
    mutate(p_adj = p.adjust(p, method = "BH"),
           sig = case_when(p_adj < 0.001 ~ "***",
                            p_adj < 0.01  ~ "**",
                            p_adj < 0.05  ~ "*",
                            TRUE ~ ""),
           region = region_name)
}

# -----------------------------------------------------------------------------
# 4. Marker validation helpers
# -----------------------------------------------------------------------------

#' AddMetaData a Seurat subset with its RCTD weights for one sample
get_sample_with_weights <- function(sample_id, merged_obj, weights_df, celltype_cols,
                                     regions = c("DG", "ML")) {
  sp_sub <- subset(merged_obj, subset = sample == sample_id & region %in% regions)

  w <- weights_df %>% filter(sample == sample_id)
  rownames(w) <- w$barcode
  w <- w[colnames(sp_sub), celltype_cols, drop = FALSE]

  AddMetaData(sp_sub, metadata = w)
}

#' Thin wrapper around cor.test(method = "spearman") returning a one-row tibble
spearman_cor <- function(x, y) {
  test <- suppressWarnings(cor.test(x, y, method = "spearman"))
  tibble(rho = unname(test$estimate), p_value = test$p.value)
}

# -----------------------------------------------------------------------------
# 5. Ligand-receptor vs colocalization: core statistics
# -----------------------------------------------------------------------------

#' Colocalization strength for a pair of cell types: product of their RCTD
#' weights, per spot. NA automatically where either cell type wasn't a
#' candidate for that spot's region (e.g. Cajal-Retzius is ML-only), which
#' is exactly the set of spots that should be excluded from the pair's
#' downstream analysis.
compute_coloc_strength <- function(weights_df, ct1, ct2) {
  stopifnot(all(c(ct1, ct2, "barcode") %in% colnames(weights_df)))
  setNames(weights_df[[ct1]] * weights_df[[ct2]], weights_df$barcode)
}

#' Same-spot ligand-receptor co-expression, partial Spearman correlation
#' controlling for a covariate (nUMI, to remove the sequencing-depth confound).
safe_pcor <- function(ligand, receptor, expr_mat, covariate) {
  if (!(ligand %in% rownames(expr_mat)) || !(receptor %in% rownames(expr_mat)))
    return(tibble(rho = NA_real_, p = NA_real_))
  pc <- pcor.test(expr_mat[ligand, ], expr_mat[receptor, ], covariate, method = "spearman")
  tibble(rho = pc$estimate, p = pc$p.value)
}

#' LR product (ligand x receptor expression) vs colocalization strength,
#' partial Spearman correlation controlling for a covariate (nUMI).
safe_pcor_coloc <- function(ligand, receptor, expr_mat, coloc, covariate) {
  if (!(ligand %in% rownames(expr_mat)) || !(receptor %in% rownames(expr_mat)))
    return(tibble(rho_vs_coloc = NA_real_, p = NA_real_))
  lr_product <- expr_mat[ligand, ] * expr_mat[receptor, ]
  pc <- pcor.test(lr_product, coloc, covariate, method = "spearman")
  tibble(rho_vs_coloc = pc$estimate, p = pc$p.value)
}

#' Rank-based ANCOVA testing whether the (coloc -> LR product) relationship
#' differs by age. Reference level is "young", so coloc_est is the young-group
#' slope and *_int_est are deviations from it (add them to get the other
#' groups' own slopes — see slope_young/middle/old in run_lr_coloc_pipeline()).
safe_lm_interaction <- function(ligand, receptor, expr_mat, coloc, nUMI_vec, age_vec) {
  na_out <- tibble(coloc_est = NA_real_, coloc_p = NA_real_,
                    mid_int_est = NA_real_, mid_int_p = NA_real_,
                    old_int_est = NA_real_, old_int_p = NA_real_)
  if (!(ligand %in% rownames(expr_mat)) || !(receptor %in% rownames(expr_mat))) return(na_out)

  lr_product <- expr_mat[ligand, ] * expr_mat[receptor, ]
  df <- tibble(
    lr_product = rank(lr_product),
    coloc = rank(coloc),
    nUMI = rank(nUMI_vec),
    age = factor(age_vec, levels = c("young", "middle", "old"))
  )
  fit <- lm(lr_product ~ coloc * age + nUMI, data = df)
  s <- summary(fit)$coefficients

  get_coef <- function(name) {
    if (name %in% rownames(s)) c(s[name, "Estimate"], s[name, "Pr(>|t|)"]) else c(NA_real_, NA_real_)
  }
  coloc_c <- get_coef("coloc")
  mid_c   <- get_coef("coloc:agemiddle")
  old_c   <- get_coef("coloc:ageold")

  tibble(coloc_est = coloc_c[1], coloc_p = coloc_c[2],
         mid_int_est = mid_c[1], mid_int_p = mid_c[2],
         old_int_est = old_c[1], old_int_p = old_c[2])
}

#' Fraction of given spots where BOTH ligand and receptor are detected (> 0)
both_detect_rate <- function(ligand, receptor, barcodes, expr_mat) {
  if (!(ligand %in% rownames(expr_mat)) || !(receptor %in% rownames(expr_mat))) return(NA_real_)
  mean(expr_mat[ligand, barcodes] > 0 & expr_mat[receptor, barcodes] > 0)
}

#' Detection enrichment for one age group: co-detection rate in the top
#' quartile of coloc_strength minus the bottom quartile, WITHIN that age
#' group only (the quantile split must not pool ages — pooling would let
#' one age's spots dominate the hi/lo cutoffs for pairs selected precisely
#' because their relationship to coloc_strength differs by age).
compute_detection_by_age <- function(age_label, pair_df, expr_mat, coloc_strength, age_vec) {
  age_barcodes <- names(age_vec)[age_vec == age_label]
  age_barcodes <- intersect(age_barcodes, colnames(expr_mat))
  coloc_age <- coloc_strength[age_barcodes]

  hi_thresh <- quantile(coloc_age, 0.75, na.rm = TRUE)
  lo_thresh <- quantile(coloc_age, 0.25, na.rm = TRUE)
  hi_barcodes <- names(coloc_age)[coloc_age >= hi_thresh]
  lo_barcodes <- names(coloc_age)[coloc_age <= lo_thresh]

  pair_df %>%
    mutate(
      pair = paste0(source_genesymbol, "-", target_genesymbol),
      both_hi = map2_dbl(source_genesymbol, target_genesymbol, both_detect_rate,
                          barcodes = hi_barcodes, expr_mat = expr_mat),
      both_lo = map2_dbl(source_genesymbol, target_genesymbol, both_detect_rate,
                          barcodes = lo_barcodes, expr_mat = expr_mat),
      enrichment = both_hi - both_lo,
      age = age_label
    )
}

#' Reversal / fading / strengthening classification from a young value to an old value
classify_pattern <- function(young, old) {
  case_when(
    sign(young) != sign(old) & young != 0 & old != 0 ~ "reversal",
    abs(old) > abs(young) ~ "strengthening (same direction)",
    abs(old) < abs(young) ~ "fading (same direction)",
    TRUE ~ "unchanged"
  )
}

# -----------------------------------------------------------------------------
# 6. Ligand-receptor vs colocalization: one-call-per-pair pipeline
# -----------------------------------------------------------------------------

#' Run the full LR x colocalization pipeline for one pair of cell types.
#'
#' This is the function that makes the analysis reusable across cell-type
#' pairs (Astrocyte/Vascular, Vascular/Microglia, OPC/Cajal-Retzius, ...):
#' call it once per pair, then feed the returned pieces to the plotting
#' helpers below.
#'
#' @param weights_df   RCTD weights table (columns: barcode, age, region, one column per cell type)
#' @param expr_mat     normalized expression matrix (genes x spots), e.g. SCT "data", full spot set
#' @param lr_resource  tibble with source_genesymbol / target_genesymbol columns (mouse symbols)
#' @param ct1, ct2     names of the two RCTD weight columns to multiply into a colocalization score
#' @param nUMI_vec     named vector (names = spot barcodes), sequencing depth per spot
#' @param age_vec      named vector (names = spot barcodes), age group per spot
#'
#' @return a list: pair_label, coloc_strength, expr_sub, nUMI, age, same_spot,
#'   lr_vs_coloc, same_spot_vs_coloc (joined quadrant table), age_model,
#'   candidate_pairs, detection_by_age, combined_shift, n_spots_used
run_lr_coloc_pipeline <- function(weights_df, expr_mat, lr_resource, ct1, ct2, nUMI_vec, age_vec) {

  stopifnot(all(c(ct1, ct2) %in% colnames(weights_df)))

  coloc_strength <- compute_coloc_strength(weights_df, ct1, ct2)

  # restrict to spots where BOTH cell types were valid RCTD candidates —
  # drops NA automatically (e.g. Cajal-Retzius is only a candidate in ML,
  # so DG spots are excluded here rather than crashing pcor.test downstream)
  barcodes_used <- intersect(names(coloc_strength)[!is.na(coloc_strength)], colnames(expr_mat))
  coloc_strength <- coloc_strength[barcodes_used]
  expr_sub <- expr_mat[, barcodes_used]
  nUMI_sub <- nUMI_vec[barcodes_used]
  age_sub  <- age_vec[barcodes_used]

  same_spot <- lr_resource %>%
    mutate(test = map2(source_genesymbol, target_genesymbol, safe_pcor,
                        expr_mat = expr_sub, covariate = nUMI_sub)) %>%
    unnest(test) %>% filter(!is.na(rho)) %>%
    mutate(p_adj = p.adjust(p, method = "BH")) %>%
    arrange(desc(rho))

  lr_vs_coloc <- lr_resource %>%
    mutate(test = map2(source_genesymbol, target_genesymbol, safe_pcor_coloc,
                        expr_mat = expr_sub, coloc = coloc_strength, covariate = nUMI_sub)) %>%
    unnest(test) %>% filter(!is.na(rho_vs_coloc)) %>%
    mutate(p_adj = p.adjust(p, method = "BH")) %>%
    arrange(desc(rho_vs_coloc))

  same_spot_vs_coloc <- same_spot %>%
    dplyr::select(source_genesymbol, target_genesymbol, rho_same_spot = rho, p_adj_same_spot = p_adj) %>%
    inner_join(
      lr_vs_coloc %>% dplyr::select(source_genesymbol, target_genesymbol, rho_vs_coloc, p_adj_vs_coloc = p_adj),
      by = c("source_genesymbol", "target_genesymbol")
    ) %>%
    mutate(pair = paste0(source_genesymbol, "-", target_genesymbol),
           both_sig = p_adj_same_spot < 0.05 & p_adj_vs_coloc < 0.05)

  age_model <- lr_resource %>%
    mutate(test = map2(source_genesymbol, target_genesymbol, safe_lm_interaction,
                        expr_mat = expr_sub, coloc = coloc_strength,
                        nUMI_vec = nUMI_sub, age_vec = age_sub)) %>%
    unnest(test) %>%
    filter(!is.na(coloc_p)) %>%
    mutate(coloc_p_adj = p.adjust(coloc_p, method = "BH"),
           mid_int_p_adj = p.adjust(mid_int_p, method = "BH"),
           old_int_p_adj = p.adjust(old_int_p, method = "BH"),
           slope_young  = coloc_est,
           slope_middle = coloc_est + mid_int_est,
           slope_old    = coloc_est + old_int_est,
           pair = paste0(source_genesymbol, "-", target_genesymbol),
           sig_old = old_int_p_adj < 0.05,
           sig_mid = mid_int_p_adj < 0.05)

  candidate_pairs <- age_model %>%
    filter(sig_old | sig_mid) %>%
    distinct(source_genesymbol, target_genesymbol, pair)

  detection_by_age <- map_dfr(c("young", "middle", "old"), compute_detection_by_age,
                               pair_df = candidate_pairs, expr_mat = expr_sub,
                               coloc_strength = coloc_strength, age_vec = age_sub)

  enrichment_wide <- detection_by_age %>%
    dplyr::select(pair, age, enrichment) %>%
    pivot_wider(names_from = age, values_from = enrichment, names_prefix = "enrich_")

  combined_shift <- age_model %>%
    dplyr::select(pair, slope_young, slope_middle, slope_old, old_int_p_adj, mid_int_p_adj) %>%
    inner_join(enrichment_wide, by = "pair") %>%
    mutate(
      slope_shift = slope_old - slope_young,
      enrichment_shift = enrich_old - enrich_young,
      same_direction = sign(slope_shift) == sign(enrichment_shift),
      combined_magnitude = abs(slope_shift) * abs(enrichment_shift)
    )

  list(
    ct1 = ct1, ct2 = ct2,
    pair_label = paste0(ct1, "-", ct2),
    coloc_strength = coloc_strength,
    expr_sub = expr_sub,
    nUMI = nUMI_sub,
    age = age_sub,
    same_spot = same_spot,
    lr_vs_coloc = lr_vs_coloc,
    same_spot_vs_coloc = same_spot_vs_coloc,
    age_model = age_model,
    candidate_pairs = candidate_pairs,
    detection_by_age = detection_by_age,
    combined_shift = combined_shift,
    n_spots_used = length(barcodes_used)
  )
}

# -----------------------------------------------------------------------------
# 7. Plotting helpers (generic — take a pair_label for titles instead of
#    hardcoding cell-type names, so the same call works for any pair)
# -----------------------------------------------------------------------------

#' Top N positive and top N negative LR pairs vs coloc strength, as a bar chart
plot_top_lr_bar <- function(lr_vs_coloc, pair_label, n = 10) {
  df <- lr_vs_coloc %>%
    filter(p_adj < 0.05) %>%
    mutate(direction = ifelse(rho_vs_coloc > 0, "Positive", "Negative"),
           pair = paste0(source_genesymbol, " -> ", target_genesymbol)) %>%
    group_by(direction) %>%
    slice_max(order_by = abs(rho_vs_coloc), n = n) %>%
    ungroup()

  ggplot(df, aes(x = rho_vs_coloc, y = fct_reorder(pair, rho_vs_coloc), fill = direction)) +
    geom_col(width = 0.65) +
    scale_fill_manual(values = SIG_DIRECTION_COLORS, name = NULL) +
    labs(x = "Partial Spearman rho (LR product vs coloc strength, adj. for nUMI)",
         y = NULL,
         title = paste0("Top ", n, " positive/negative LR pairs vs coloc strength: ", pair_label)) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.position = "top")
}

#' Volcano plot of every tested LR pair vs coloc strength
plot_lr_volcano <- function(lr_vs_coloc, pair_label, n_label = 20) {
  df <- lr_vs_coloc %>%
    mutate(neglog10p = -log10(p_adj),
           sig = p_adj < 0.05,
           pair = paste0(source_genesymbol, "-", target_genesymbol))

  top_labels <- df %>%
    filter(sig) %>%
    mutate(direction = ifelse(rho_vs_coloc > 0, "Positive", "Negative")) %>%
    group_by(direction) %>%
    slice_max(order_by = abs(rho_vs_coloc), n = n_label) %>%
    ungroup()

  df <- df %>% mutate(label = ifelse(pair %in% top_labels$pair, pair, NA))

  ggplot(df, aes(x = rho_vs_coloc, y = neglog10p, color = sig)) +
    geom_point(alpha = 0.6, size = 1.8) +
    ggrepel::geom_text_repel(aes(label = label), size = 3, max.overlaps = 30, color = "black",
                              segment.size = 0.3, segment.color = "grey50") +
    scale_color_manual(values = c(`TRUE` = SIG_COLOR, `FALSE` = "grey70"), guide = "none") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    labs(x = "Partial Spearman rho (LR product vs coloc strength, adj. for nUMI)",
         y = expression(-log[10](p[adj])),
         title = paste0("All tested LR pairs vs ", pair_label, " colocalization strength"),
         subtitle = paste0("Labeled: top ", n_label, " positive and top ", n_label, " negative by |rho| (p_adj < 0.05)")) +
    theme_minimal(base_size = 12)
}

#' Scatter of one LR pair's product vs coloc strength, with its adjusted rho in the subtitle
plot_lr_vs_coloc_scatter <- function(ligand, receptor, expr_mat, coloc_strength, rho_table, pair_label) {
  lr_product <- expr_mat[ligand, ] * expr_mat[receptor, ]
  df <- tibble(lr_product = lr_product, coloc_strength = coloc_strength[names(lr_product)])
  rho_adj <- rho_table %>%
    filter(source_genesymbol == ligand, target_genesymbol == receptor) %>%
    pull(rho_vs_coloc)

  ggplot(df, aes(x = coloc_strength, y = lr_product)) +
    geom_point(alpha = 0.4, color = "#0072B2") +
    geom_smooth(method = "loess", color = "#D55E00", se = TRUE, linewidth = 1) +
    labs(x = paste0(pair_label, " RCTD weight product (per spot)"),
         y = paste0(ligand, " x ", receptor, " expression product"),
         title = paste0(ligand, "-", receptor, " vs ", pair_label, " coloc strength"),
         subtitle = paste0("Partial Spearman rho (adj. for nUMI) = ", round(rho_adj, 3))) +
    theme_minimal(base_size = 12)
}

#' Generic quadrant scatter with dashed zero-lines, optional y=x diagonal,
#' quadrant text labels, and repel-labeled significant points. Used both for
#' same-spot-rho vs coloc-tracking-rho, and for slope_young vs slope_old.
plot_quadrant_scatter <- function(df, x, y, sig_col, quad_labels, xlab, ylab, title,
                                   subtitle = NULL, diagonal = FALSE, label_col = "pair") {
  p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]], color = .data[[sig_col]])) +
    geom_point(alpha = 0.6, size = 1.8)

  if (diagonal) {
    p <- p + geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey40")
  }

  p +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_text(data = quad_labels, aes(x = x, y = y, label = label), inherit.aes = FALSE,
              size = 3, color = "grey40", fontface = "italic") +
    ggrepel::geom_text_repel(data = filter(df, .data[[sig_col]]), aes(label = .data[[label_col]]),
                              size = 3, color = "black", max.overlaps = 15) +
    scale_color_manual(values = c(`TRUE` = SIG_COLOR, `FALSE` = "grey70"), guide = "none") +
    labs(x = xlab, y = ylab, title = title, subtitle = subtitle) +
    theme_minimal(base_size = 12)
}

#' Build the long-format age-trajectory table for the N pairs with the
#' largest young-to-old shift in a single metric (e.g. slope alone).
build_age_trajectory_df <- function(age_model, sig_col = "sig_old", n = 12) {
  top_pairs <- age_model %>%
    dplyr::filter(.data[[sig_col]]) %>%
    mutate(shift_magnitude = abs(slope_old - slope_young)) %>%
    slice_max(order_by = shift_magnitude, n = n) %>%
    pull(pair)

  age_model %>%
    dplyr::filter(pair %in% top_pairs) %>%
    dplyr::select(pair, slope_young, slope_middle, slope_old) %>%
    pivot_longer(-pair, names_to = "age", values_to = "value") %>%
    mutate(age = factor(recode(age, slope_young = "young", slope_middle = "middle", slope_old = "old"),
                         levels = c("young", "middle", "old")),
           metric = "slope",
           pair = fct_reorder(pair, value, .fun = function(v) abs(diff(range(v))), .desc = TRUE))
}

#' Build the long-format trajectory table (slope + detection enrichment) for
#' the top-N pairs by combined magnitude, both metrics moving in the same direction.
build_combined_shift_trajectory_df <- function(combined_shift, n = 30, p_thresh = 0.05) {
  top_combined <- combined_shift %>%
    filter(same_direction, old_int_p_adj < p_thresh) %>%
    arrange(desc(combined_magnitude)) %>%
    slice_head(n = n)

  slope_traj <- top_combined %>%
    dplyr::select(pair, young = slope_young, middle = slope_middle, old = slope_old) %>%
    pivot_longer(-pair, names_to = "age", values_to = "value") %>%
    mutate(metric = "Coloc -> LR product slope")

  enrich_traj <- top_combined %>%
    dplyr::select(pair, young = enrich_young, old = enrich_old) %>%
    pivot_longer(-pair, names_to = "age", values_to = "value") %>%
    mutate(metric = "Detection enrichment")

  bind_rows(slope_traj, enrich_traj) %>%
    mutate(age = factor(age, levels = c("young", "middle", "old")),
           pair = fct_reorder(pair, value, .fun = function(v) max(v) - min(v), .desc = TRUE))
}

#' Generic multi-panel (one facet per pair) trajectory-across-age line plot.
#' Works for a single metric (e.g. slope alone) or several (color = metric).
plot_metric_trajectory <- function(long_df, color_map, title, ncol = 5, y_lab = NULL,
                                    show_legend = length(color_map) > 1) {
  ggplot(long_df, aes(x = age, y = value, group = interaction(pair, metric), color = metric)) +
    geom_line(alpha = 0.8, linewidth = 0.8) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    facet_wrap(~pair, scales = "free_y", ncol = ncol) +
    scale_color_manual(values = color_map) +
    labs(x = NULL, y = y_lab, title = title) +
    theme_minimal(base_size = 11) +
    theme(legend.position = if (show_legend) "top" else "none",
          panel.grid.minor = element_blank(),
          strip.text = element_text(size = 10, face = "bold"))
}

#' Quadrant text labels for the same-spot-rho vs coloc-tracking-rho scatter.
#' `scale` places the labels a fraction of the way to each corner — adjust
#' per pair if the rho range differs a lot from the Astrocyte-Vascular case.
make_same_spot_quad_labels <- function(scale = 0.6) {
  tibble(
    x = c(scale, -scale, scale, -scale),
    y = c(scale, scale, -scale, -scale),
    label = c("Co-expressed &\nenriched in coloc niche",
              "Not co-expressed,\nbut tracks coloc (check)",
              "Co-expressed,\nbut elsewhere in tissue",
              "Neither\n(background)")
  )
}

#' Quadrant text labels for the slope_young vs slope_old scatter.
make_age_shift_quad_labels <- function(scale = 0.2) {
  tibble(
    x = c(scale, -scale, scale, -scale),
    y = c(scale, scale, -scale, -scale),
    label = c("Positive in both\n(stable, coupled)",
              "Negative in young,\npositive in old (emerges)",
              "Positive in young,\nnegative in old (lost)",
              "Negative in both\n(stable, excluded)")
  )
}

#' Violin + boxplot of a numeric column across age, faceted (free y scale)
plot_violin_by_age <- function(df, y_col, facet_col, title, y_lab, age_colors = AGE_COLORS) {
  ggplot(df, aes(x = age, y = .data[[y_col]], fill = age)) +
    geom_violin(trim = FALSE, alpha = 0.85, linewidth = 0.3) +
    geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA) +
    facet_wrap(vars(.data[[facet_col]]), scales = "free_y") +
    scale_fill_manual(values = age_colors) +
    labs(x = NULL, y = y_lab, title = title) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none", panel.grid.minor = element_blank())
}

# -----------------------------------------------------------------------------
# 8. Cross-modality comparison: spatial LR-coloc pipeline vs LIANA+ scRNA CCC
# -----------------------------------------------------------------------------

#' Compare one spatial LR-coloc pipeline result (from run_lr_coloc_pipeline())
#' to its matching age-split LIANA+ scRNA CCC result (a CSV with columns
#' source, target, ligand_complex, receptor_complex, magnitude_rank,
#' specificity_rank, age), and check 3-way sign agreement on the young-to-old
#' shift across spatial slope, spatial detection enrichment, and scRNA evidence.
compare_spatial_to_liana <- function(spatial_result, liana_csv_path) {
  liana_by_age <- read_csv(liana_csv_path, show_col_types = FALSE) %>%
    mutate(pair = paste0(ligand_complex, "-", receptor_complex))

  # a pair can appear twice per age (both interaction directions, e.g.
  # Vasculature->Microglia and Microglia->Vasculature) — keep the direction
  # with the strongest evidence, but record which direction that was
  # (source = the cell type expressing the ligand, target = the receptor)
  # so downstream tables can report which way the signal is predicted to flow
  liana_by_age_dedup <- liana_by_age %>%
    group_by(pair, age) %>%
    slice_min(magnitude_rank, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(direction = paste0(source, " -> ", target))

  slope_long <- spatial_result$age_model %>%
    dplyr::select(pair, slope_young, slope_middle, slope_old) %>%
    pivot_longer(-pair, names_to = "age", values_to = "slope") %>%
    mutate(age = recode(age, slope_young = "young", slope_middle = "middle", slope_old = "old"))

  combined_evidence <- slope_long %>%
    inner_join(spatial_result$detection_by_age %>% dplyr::select(pair, age, enrichment), by = c("pair", "age")) %>%
    inner_join(liana_by_age_dedup %>% dplyr::select(pair, age, ligand_complex, receptor_complex, magnitude_rank, specificity_rank, direction), by = c("pair", "age")) %>%
    mutate(
      scrna_evidence = -log10(magnitude_rank),   # flip direction so higher = stronger evidence, matching slope/enrichment's convention
      age = factor(age, levels = c("young", "middle", "old"))
    )

  # keep only pairs with all 3 ages present in every joined table
  pairs_complete <- combined_evidence %>%
    group_by(pair) %>%
    filter(n_distinct(age) == 3) %>%
    ungroup()

  dropped_pairs <- setdiff(unique(combined_evidence$pair), unique(pairs_complete$pair))

  agreement_check <- pairs_complete %>%
    group_by(pair) %>%
    summarize(
      slope_shift = slope[age == "old"] - slope[age == "young"],
      enrichment_shift = enrichment[age == "old"] - enrichment[age == "young"],
      scrna_shift = scrna_evidence[age == "old"] - scrna_evidence[age == "young"],
      spatial_agrees = sign(slope_shift) == sign(enrichment_shift),
      scrna_agrees = sign(slope_shift) == sign(scrna_shift),
      .groups = "drop"
    ) %>%
    filter(spatial_agrees, scrna_agrees) %>%
    arrange(desc(abs(slope_shift) + abs(scrna_shift)))

  # which direction (source -> target) had the stronger scRNA evidence, per age —
  # lets you state e.g. "Ptn-Sdc4 is Astrocyte -> Vascular in young and old,
  # not the reverse" directly from the LIANA+ output, rather than assuming a
  # direction from gene identity alone
  direction_wide <- combined_evidence %>%
    filter(pair %in% agreement_check$pair) %>%
    dplyr::select(pair, age, direction) %>%
    pivot_wider(id_cols = pair, names_from = age, values_from = direction, names_glue = "direction_{age}")

  # ligand_complex/receptor_complex are constant for a given "pair" string
  # (pair is literally paste0(ligand_complex, "-", receptor_complex)) -- grab
  # one row per pair rather than splitting the "pair" string back apart later,
  # since some mouse gene symbols (e.g. the H2-* MHC genes) contain a literal
  # "-" and would make that split ambiguous.
  complex_lookup <- combined_evidence %>%
    distinct(pair, ligand_complex, receptor_complex)

  pair_detail_wide <- combined_evidence %>%
    filter(pair %in% agreement_check$pair) %>%
    dplyr::select(pair, age, slope, enrichment, scrna_evidence) %>%
    pivot_wider(id_cols = pair, names_from = age, values_from = c(slope, enrichment, scrna_evidence),
                names_glue = "{.value}_{age}") %>%
    left_join(agreement_check %>% dplyr::select(pair, slope_shift, enrichment_shift, scrna_shift), by = "pair") %>%
    left_join(direction_wide, by = "pair") %>%
    left_join(complex_lookup, by = "pair") %>%
    arrange(desc(abs(slope_shift) + abs(scrna_shift))) %>%
    mutate(
      slope_pattern = classify_pattern(slope_young, slope_old),
      enrichment_pattern = classify_pattern(enrichment_young, enrichment_old),
      # TRUE if the dominant direction is the same across all 3 ages; FALSE flags
      # a pair where the predicted sender/receiver actually flips with age
      direction_consistent = (direction_young == direction_middle) & (direction_middle == direction_old)
    )

  list(
    combined_evidence = combined_evidence,
    pairs_dropped_incomplete = dropped_pairs,
    agreement_check = agreement_check,
    pair_detail_wide = pair_detail_wide
  )
}

#' Trajectory plot (slope / enrichment / scRNA evidence) for the pairs where
#' all 3 metrics agree on the direction of their young-to-old shift.
plot_cross_modality_trajectory <- function(compare_result, title, ncol = 5) {
  traj_three <- compare_result$combined_evidence %>%
    filter(pair %in% compare_result$agreement_check$pair) %>%
    dplyr::select(pair, age, slope, enrichment, scrna_evidence) %>%
    pivot_longer(c(slope, enrichment, scrna_evidence), names_to = "metric", values_to = "value") %>%
    mutate(metric = recode(metric,
                            slope = "Spatial: coloc -> LR slope",
                            enrichment = "Spatial: detection enrichment",
                            scrna_evidence = "scRNA: -log10(magnitude_rank)"))

  plot_metric_trajectory(traj_three,
                          color_map = c("Spatial: coloc -> LR slope" = "#0072B2",
                                        "Spatial: detection enrichment" = "#D55E00",
                                        "scRNA: -log10(magnitude_rank)" = "#009E73"),
                          title = title, ncol = ncol)
}

# -----------------------------------------------------------------------------
# 9. Spatial visualization of high/low colocalization spots
# -----------------------------------------------------------------------------

#' Classify spots into top/bottom colocalization quartile, computed
#' SEPARATELY within each age group — the same convention used by
#' compute_detection_by_age(), so "High"/"Low" here are exactly the
#' hi_barcodes/lo_barcodes that feed the detection-enrichment analysis.
classify_coloc_quartile <- function(coloc_strength, age_vec) {
  barcodes <- names(coloc_strength)
  age_sub  <- age_vec[barcodes]

  group <- rep("Mid", length(barcodes))
  names(group) <- barcodes

  for (age_label in unique(age_sub)) {
    age_barcodes <- barcodes[age_sub == age_label]
    coloc_age <- coloc_strength[age_barcodes]

    hi_thresh <- quantile(coloc_age, 0.75, na.rm = TRUE)
    lo_thresh <- quantile(coloc_age, 0.25, na.rm = TRUE)

    group[age_barcodes[coloc_age >= hi_thresh]] <- "High"
    group[age_barcodes[coloc_age <= lo_thresh]] <- "Low"
  }

  factor(group, levels = c("Low", "Mid", "High"))
}

#' Spatial plot of one sample's spots, colored by top/bottom colocalization
#' quartile for a given cell-type-pair pipeline result. Automatically
#' restricted to whichever spots that pipeline actually used (e.g. ML-only
#' for OPC-Cajal-Retzius), since coloc_strength is NA everywhere else.
plot_coloc_spatial <- function(sample_id, spatial_result, merged_obj, regions = c("DG", "ML")) {
  sp_sub <- subset(merged_obj, subset = sample == sample_id & region %in% regions)

  coloc_group <- classify_coloc_quartile(spatial_result$coloc_strength, spatial_result$age)
  keep_barcodes <- intersect(names(coloc_group), colnames(sp_sub))

  sp_sub <- subset(sp_sub, cells = keep_barcodes)
  sp_sub <- AddMetaData(sp_sub, metadata = coloc_group[keep_barcodes], col.name = "coloc_group")

  SpatialDimPlot(sp_sub, group.by = "coloc_group") +
    scale_fill_manual(values = c(Low = "#0072B2", Mid = "grey85", High = "#D55E00"), name = NULL) +
    ggtitle(paste0(sample_id, ": ", spatial_result$pair_label))
}

#' Cross-tabulate two cell-type-pair pipelines' High/Mid/Low colocalization
#' classification, spot by spot, over whichever spots both pipelines used.
#' Reports the full contingency table plus Jaccard overlap for the two most
#' interpretable cells (High-High: spots colocalized for both pairs at once;
#' Low-Low: excluded from both) and a chi-squared test of independence as a
#' quick check of whether the two classifications associate more than chance.
compute_coloc_overlap <- function(result1, result2, label1 = result1$pair_label, label2 = result2$pair_label) {
  group1 <- classify_coloc_quartile(result1$coloc_strength, result1$age)
  group2 <- classify_coloc_quartile(result2$coloc_strength, result2$age)

  shared_barcodes <- intersect(names(group1), names(group2))

  overlap_df <- tibble(
    barcode = shared_barcodes,
    group1  = group1[shared_barcodes],
    group2  = group2[shared_barcodes],
    age     = result1$age[shared_barcodes]
  )

  contingency <- table(group1 = overlap_df$group1, group2 = overlap_df$group2)

  jaccard <- function(a, b) {
    inter <- sum(a & b); uni <- sum(a | b)
    if (uni == 0) return(NA_real_)
    inter / uni
  }

  chisq <- suppressWarnings(chisq.test(contingency))

  list(
    label1 = label1, label2 = label2,
    overlap_df = overlap_df,
    n_shared_spots = length(shared_barcodes),
    contingency = contingency,
    proportions = prop.table(contingency),
    jaccard_high_high = jaccard(overlap_df$group1 == "High", overlap_df$group2 == "High"),
    jaccard_low_low   = jaccard(overlap_df$group1 == "Low",  overlap_df$group2 == "Low"),
    chisq_p = chisq$p.value
  )
}

#' Heatmap of a compute_coloc_overlap() contingency table, with the Jaccard
#' overlaps and shared-spot count in the subtitle.
plot_coloc_overlap <- function(overlap_result) {
  df <- as.data.frame(overlap_result$contingency)

  ggplot(df, aes(x = group2, y = group1, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq), size = 4) +
    scale_fill_gradient(low = "white", high = "#0072B2", name = "spots") +
    labs(x = overlap_result$label2, y = overlap_result$label1,
         title = paste0("Spot overlap: ", overlap_result$label1, " vs ", overlap_result$label2),
         subtitle = paste0("n = ", overlap_result$n_shared_spots, " shared spots  |  ",
                            "Jaccard(High-High) = ", round(overlap_result$jaccard_high_high, 2),
                            ", Jaccard(Low-Low) = ", round(overlap_result$jaccard_low_low, 2))) +
    theme_minimal(base_size = 12)
}
