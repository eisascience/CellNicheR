
#' Create a proportional scale bar for spatial CosMx plots
#'
#' Generates a list of `annotation_custom()` layers that can be added to a ggplot
#' to display a scale bar with appropriate length and units (µm or mm) based on
#' pixel coordinates. The scale bar automatically adjusts to roughly one-quarter
#' of the x-axis span.
#'
#' @param x_vals Numeric vector of x coordinates (in pixels).
#' @param y_vals Numeric vector of y coordinates (in pixels).
#' @param microns_per_pixel Numeric. Conversion factor from pixels to microns.
#'   Default is `0.12028`, corresponding to the commercial CosMx platform.
#'
#' @return A named list with three ggplot annotation layers:
#'   \itemize{
#'     \item `bg` — a semi-transparent white rectangle background.
#'     \item `rect` — a solid black rectangle representing the scale bar.
#'     \item `label` — a text label indicating the physical length (µm or mm).
#'   }
#'
#' @details
#' This helper is designed for consistent inclusion of scale bars across
#' CosMx spatial plots. It estimates a visually balanced bar length based on
#' the x-axis range and applies an appropriate unit label automatically.
#'
#' @examples
#' scale_bar <- make_scale_bar_r(
#'   x_vals = cell_meta$CenterX_global_px,
#'   y_vals = cell_meta$CenterY_global_px
#' )
#' ggplot() + scale_bar$bg + scale_bar$rect + scale_bar$label
#'
#' @export
make_scale_bar_r <- function(x_vals, y_vals, microns_per_pixel = 0.12028) {
  
  # Adds a scale bar to a ggplot.
  
  #Parameters:
  #x_vals: vector of x coordinates in pixels
  #y_vals: vector of y coordinates in pixels
  #microns_per_pixel: conversion factor from pixels to microns. Default is conversion factor for commercial CosMx instrument
  
  # Example usage:
  # scale_bar = make_scale-bar_r(x_vals = cell_meta$CenterX_global_px, y_vals = cell_meta$CenterY_global_px)
  # ggplot() + scale_bar$bg + scale_bar$rect + scale_bar$label
  
  
  # Calculate x-axis range
  x_range <- range(x_vals, na.rm = TRUE)
  x_length <- diff(x_range)
  x_length_um <- x_length * microns_per_pixel
  
  # Target scale length ~1/4 of the x-axis
  target <- x_length_um / 4
  
  # Compute order of magnitude
  order <- 10^floor(log10(target))
  mantissa <- target / order
  
  # Round mantissa to nearest 1, 2, or 5
  nice_mantissa <- if (mantissa < 1.5) {
    1
  } else if (mantissa < 3.5) {
    2
  } else if (mantissa < 7.5) {
    5
  } else {
    10
  }
  
  # Final scale length in pixels
  scale_length_um <- nice_mantissa * order
  scale_length_px <- scale_length_um / microns_per_pixel
  
  # Format label
  scale_label <- if (scale_length_um >= 1000) {
    paste0(scale_length_um / 1000, " mm")
  } else {
    paste0(scale_length_um, " µm")
  }
  
  # Set coordinates for the scale bar 
  x_start <- x_range[2] - scale_length_px * 1.1
  x_end <- x_range[2] - scale_length_px * 0.1
  y_pos <- min(y_vals, na.rm = TRUE) + scale_length_px * 0.1
  
  # Generate scale bar background, scale bar and annotation to return
  list(
    bg = annotation_custom(
      grob = rectGrob(gp = gpar(fill = "white", alpha = 0.8, col = NA)),
      xmin = x_start- scale_length_px*0.05, xmax = x_end+ scale_length_px*0.05,
      ymin = y_pos- scale_length_px*0.05, ymax = y_pos + scale_length_px*0.3
    ),
    rect = annotation_custom(
      grob = rectGrob(gp = gpar(fill = "black")),
      xmin = x_start, xmax = x_end,
      ymin = y_pos, ymax = y_pos + scale_length_px * 0.05
    ),
    label = annotation_custom(
      grob = textGrob(scale_label, gp = gpar(col = "black"), just = "center", vjust = 0),
      xmin = (x_start + x_end)/2, xmax = (x_start + x_end)/2,
      ymin = y_pos + scale_length_px * 0.1, ymax = y_pos + scale_length_px * 0.1
    )
  )
}



#' Plot a cell neighborhood around a cell of interest (COI)
#'
#' Builds a local spatial neighborhood view around a selected cell of interest
#' and its nearest neighbors. Polygon boundaries and metadata are joined to
#' visualize viral status, cell type, or expression summaries.
#'
#' @param cell_meta Data frame of per-cell metadata containing at least
#'   `cell`, `CenterX_global_px`, and `CenterY_global_px`.
#' @param cell_polygons Data frame of per-cell polygon coordinates containing
#'   `cell`, `x_global_px`, and `y_global_px`.
#' @param expr_matrix Optional numeric expression matrix (cells × features),
#'   required if `SIVsum` is not present in `cell_meta`.
#' @param coi Character. Cell ID for the cell of interest. If `NULL`, the cell
#'   with the highest `SIVsum` is selected automatically.
#' @param k Integer. Number of nearest neighbors to include (default: `50`).
#' @param restrict_to_sample Logical. If `TRUE`, restricts the neighborhood to
#'   cells from the same `SampleID` (default: `TRUE`).
#' @param cluster_var Character. Column in `cell_meta` giving cluster or cell
#'   type labels.
#' @param fill_by Character. Variable used to color polygons: one of
#'   `"VirusPos"`, `"cluster"`, or `"SIVsum"`.
#' @param alpha_coi,alpha_bg Numeric transparency levels for the COI and
#'   background cells.
#' @param outline_coi,outline_bg Colors for polygon outlines.
#' @param linewidth Numeric. Polygon border width.
#' @param scale_bar_fun Function generating a scale bar (default:
#'   `make_scale_bar_r`). Use `NULL` to omit.
#' @param seed Integer. Random seed for reproducibility of any sampling steps.
#'
#' @return A `ggplot2` object showing the cell neighborhood polygons.
#'   The returned object includes attributes:
#'   \itemize{
#'     \item `"coi"` — the cell of interest ID.
#'     \item `"neighbors"` — vector of neighboring cell IDs.
#'   }
#'
#' @details
#' The function uses `FNN::get.knnx()` for efficient nearest-neighbor search in
#' global pixel coordinates. It automatically overlays a scale bar and adapts
#' fill colors based on the selected variable.
#'
#' @examples
#' p <- neighborhood_plot(cell_meta, cell_polygons, k = 50, fill_by = "VirusPos")
#' print(p)
#'
#' @export
neighborhood_plot <- function(
    cell_meta,
    cell_polygons,
    expr_matrix = NULL,                 # only needed if SIVsum not already in cell_meta
    coi = NULL,                         # e.g., "c_1_100_555"; if NULL choose max SIVsum
    k = 50,                             # COI + (k-1) neighbors
    restrict_to_sample = TRUE,          # only neighbors from same SampleID
    cluster_var = "RNA_Pre.Napari.v3.0.071825_Cell.Typing.InSituType.2_1_clusters",
    fill_by = c("VirusPos","cluster","SIVsum"),  # what to color polygons by
    alpha_coi = 1, alpha_bg = 0.6,
    outline_coi = "black", outline_bg = "grey40",
    linewidth = 0.25,
    scale_bar_fun = make_scale_bar_r,   # your helper; set to NULL to skip
    seed = 2025                         # for any internal sampling if added later
) {
  stopifnot(requireNamespace("dplyr", quietly = TRUE),
            requireNamespace("FNN",   quietly = TRUE),
            requireNamespace("ggplot2", quietly = TRUE))
  
  fill_by <- match.arg(fill_by)
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
  })
  
  # --- Ensure key columns exist ---
  req_meta <- c("cell","CenterX_global_px","CenterY_global_px")
  missing_meta <- setdiff(req_meta, names(cell_meta))
  if (length(missing_meta)) {
    stop("cell_meta is missing required columns: ", paste(missing_meta, collapse = ", "))
  }
  req_poly <- c("cell","x_global_px","y_global_px")
  missing_poly <- setdiff(req_poly, names(cell_polygons))
  if (length(missing_poly)) {
    stop("cell_polygons is missing required columns: ", paste(missing_poly, collapse = ", "))
  }
  
  # --- Ensure SIVsum and VirusPos exist in cell_meta (derive if needed) ---
  if (!("SIVsum" %in% names(cell_meta))) {
    if (is.null(expr_matrix)) {
      stop("SIVsum not found in cell_meta and expr_matrix is NULL. ",
           "Provide expr_matrix with a 'SIVsum' column or add SIVsum to cell_meta.")
    }
    if (!("SIVsum" %in% colnames(expr_matrix))) {
      stop("expr_matrix does not contain a 'SIVsum' column.")
    }
    # align by cell IDs
    if (!all(cell_meta$cell %in% rownames(expr_matrix))) {
      warning("Some cell_meta$cell are not present in rownames(expr_matrix); they will get NA SIVsum.")
    }
    cell_meta$SIVsum <- expr_matrix[cell_meta$cell, "SIVsum"]
  }
  if (!("VirusPos" %in% names(cell_meta))) {
    cell_meta$VirusPos <- ifelse(is.na(cell_meta$SIVsum), NA_character_,
                                 ifelse(cell_meta$SIVsum > 0, "SIVpos", "SIVneg"))
  }
  
  # --- Choose COI ---
  if (is.null(coi)) {
    coi <- cell_meta %>% arrange(dplyr::desc(SIVsum)) %>% slice(1) %>% pull(cell)
  } else {
    if (!(coi %in% cell_meta$cell)) stop("Provided 'coi' not found in cell_meta$cell.")
  }
  
  # --- Restrict to same sample (optional) ---
  if (restrict_to_sample) {
    if (!("SampleID" %in% names(cell_meta))) {
      stop("restrict_to_sample=TRUE but 'SampleID' not found in cell_meta.")
    }
    coi_sample <- cell_meta %>% filter(cell == coi) %>% pull(SampleID)
    meta_nn <- cell_meta %>% filter(SampleID == coi_sample)
  } else {
    meta_nn <- cell_meta
  }
  
  # --- Build KNN set (COI + neighbors) ---
  xy <- meta_nn %>% select(CenterX_global_px, CenterY_global_px) %>% as.matrix()
  rownames(xy) <- meta_nn$cell
  if (!(coi %in% rownames(xy))) stop("COI is not in neighborhood search set (check SampleID restriction).")
  
  k_eff <- min(k, nrow(xy)) # guard if k > n
  nn <- FNN::get.knnx(data = xy, query = xy, k = k_eff)
  rownames(nn$nn.index) <- rownames(xy)
  neighboring_cells <- rownames(nn$nn.index)[nn$nn.index[coi,]]
  neighboring_cells <- unique(c(neighboring_cells, coi))
  
  # --- Prep plotting df: polygons + selected meta ---
  # Choose available cluster var gracefully
  cluster_present <- cluster_var %in% names(cell_meta)
  meta_cols <- c("cell","VirusPos","SIVsum","SampleID")
  if (cluster_present) meta_cols <- c(meta_cols, cluster_var)
  
  meta4join <- cell_meta %>% select(all_of(intersect(meta_cols, names(cell_meta))))
  
  plot.df <- cell_polygons %>%
    filter(cell %in% neighboring_cells) %>%
    left_join(meta4join, by = "cell")
  
  # --- Build scale bar if function provided ---
  scale_bar <- NULL
  if (!is.null(scale_bar_fun)) {
    scale_bar <- scale_bar_fun(
      x_vals = plot.df$x_global_px,
      y_vals = plot.df$y_global_px
    )
  }
  
  # --- Aesthetics: fill mapping ---
  fill_mapping <- switch(
    fill_by,
    "VirusPos" = aes(fill = .data[["VirusPos"]]),
    "cluster"  = {
      if (!cluster_present) warning("cluster_var not found; falling back to VirusPos.")
      aes(fill = .data[[if (cluster_present) cluster_var else "VirusPos"]])
    },
    "SIVsum"   = aes(fill = .data[["SIVsum"]])
  )
  
  # --- Base plot ---
  p <- ggplot(
    plot.df,
    aes(x = x_global_px, y = y_global_px, group = cell,
        alpha = (cell == coi), color = (cell == coi))
  ) +
    fill_mapping +
    geom_polygon(linewidth = linewidth) +
    scale_alpha_manual(values = c("TRUE" = alpha_coi, "FALSE" = alpha_bg)) +
    scale_color_manual(values = c("TRUE" = outline_coi, "FALSE" = outline_bg)) +
    coord_fixed() +
    theme_void() +
    guides(alpha = "none", color = "none")
  
  if (fill_by == "SIVsum") {
    p <- p + scale_fill_viridis_c(option = "C")
  } else {
    p <- p + labs(fill = if (fill_by == "cluster") "Cell type" else "Viral status")
  }
  
  if (!is.null(scale_bar)) {
    p <- p + scale_bar$bg + scale_bar$rect + scale_bar$label
  }
  
  # return plot (and useful bits invisibly)
  attr(p, "coi") <- coi
  attr(p, "neighbors") <- neighboring_cells
  return(p)
}


#' Plot a gene- or feature-specific cell neighborhood
#'
#' Extends \code{\link{neighborhood_plot}} to support visualization of local
#' gene expression patterns in addition to viral status, cluster identity, or
#' SIV summary signal.
#'
#' @param cell_meta Data frame of per-cell metadata including
#'   `cell`, `CenterX_global_px`, and `CenterY_global_px`.
#' @param cell_polygons Data frame of per-cell polygon coordinates including
#'   `cell`, `x_global_px`, and `y_global_px`.
#' @param expr_matrix Optional numeric expression matrix (cells × genes), required
#'   when using `fill_by = "gene"`, or when `SIVsum` is not in `cell_meta`.
#' @param coi Character. Cell ID for the cell of interest; if `NULL`, defaults
#'   to the cell with the highest `SIVsum`.
#' @param k Integer. Number of nearest neighbors to include (default: `50`).
#' @param restrict_to_sample Logical. If `TRUE`, restricts to neighbors within
#'   the same `SampleID` (default: `TRUE`).
#' @param cluster_var Character. Column defining cell type or cluster identity.
#' @param fill_by Character. Variable used to color polygons:
#'   `"VirusPos"`, `"cluster"`, `"SIVsum"`, or `"gene"`.
#' @param gene Character. Target gene to visualize when `fill_by = "gene"`.
#' @param gene_transform Character. Transformation applied to gene expression:
#'   `"none"`, `"log1p"`, or `"sqrt"`.
#' @param alpha_coi,alpha_bg Numeric transparency values for COI and neighbors.
#' @param outline_coi,outline_bg Colors for polygon outlines.
#' @param linewidth Numeric. Border line width for polygons.
#' @param scale_bar_fun Function providing scale bar annotations (default:
#'   `make_scale_bar_r`); use `NULL` to skip.
#' @param seed Integer. Random seed for reproducibility.
#'
#' @return A `ggplot2` object showing neighborhood polygons colored by the
#'   specified variable. The plot carries attributes:
#'   \itemize{
#'     \item `"coi"` — the cell of interest.
#'     \item `"neighbors"` — vector of neighboring cell IDs.
#'     \item `"gene"` — gene name when `fill_by = "gene"`.
#'   }
#'
#' @details
#' This function provides flexible visualization of local spatial environments
#' for a given cell. When `fill_by = "gene"`, expression values are retrieved
#' from `expr_matrix` and optionally transformed for dynamic range compression.
#' Color scales adapt automatically (viridis continuous for numeric values).
#'
#' @examples
#' p_gene <- neighborhood_plot2(
#'   cell_meta, cell_polygons, expr_matrix,
#'   fill_by = "gene", gene = "CD3E", gene_transform = "log1p"
#' )
#' print(p_gene)
#'
#' @seealso [neighborhood_plot()], [make_scale_bar_r()]
#' @export
neighborhood_plot2 <- function(
    cell_meta,
    cell_polygons,
    expr_matrix = NULL,                 # needed if SIVsum or gene not already in cell_meta
    coi = NULL,                         # e.g., "c_1_100_555"; if NULL choose max SIVsum
    k = 50,                             # COI + (k-1) neighbors
    restrict_to_sample = TRUE,          # only neighbors from same SampleID
    cluster_var = "RNA_Pre.Napari.v3.0.071825_Cell.Typing.InSituType.2_1_clusters",
    fill_by = c("VirusPos","cluster","SIVsum","gene"),  # <- added "gene"
    gene = NULL,                        # <- target gene name for fill_by = "gene"
    gene_transform = c("none","log1p","sqrt"),
    alpha_coi = 1, alpha_bg = 0.6,
    outline_coi = "black", outline_bg = "grey40",
    linewidth = 0.25,
    scale_bar_fun = make_scale_bar_r,   # your helper; set to NULL to skip
    seed = 2025
) {
  stopifnot(requireNamespace("dplyr", quietly = TRUE),
            requireNamespace("FNN",   quietly = TRUE),
            requireNamespace("ggplot2", quietly = TRUE))
  fill_by <- match.arg(fill_by)
  gene_transform <- match.arg(gene_transform)
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
  })
  
  # --- Required columns ---
  req_meta <- c("cell","CenterX_global_px","CenterY_global_px")
  miss_meta <- setdiff(req_meta, names(cell_meta))
  if (length(miss_meta)) stop("cell_meta is missing: ", paste(miss_meta, collapse=", "))
  req_poly <- c("cell","x_global_px","y_global_px")
  miss_poly <- setdiff(req_poly, names(cell_polygons))
  if (length(miss_poly)) stop("cell_polygons is missing: ", paste(miss_poly, collapse=", "))
  
  # --- Ensure SIVsum/VirusPos in meta if needed ---
  if (fill_by %in% c("VirusPos","SIVsum") && !("SIVsum" %in% names(cell_meta))) {
    if (is.null(expr_matrix) || !("SIVsum" %in% colnames(expr_matrix))) {
      stop("SIVsum not found in cell_meta and expr_matrix either missing or lacks 'SIVsum'.")
    }
    if (!all(cell_meta$cell %in% rownames(expr_matrix))) {
      warning("Some cell_meta$cell are not in expr_matrix rownames; SIVsum becomes NA for them.")
    }
    cell_meta$SIVsum <- expr_matrix[cell_meta$cell, "SIVsum"]
  }
  if (fill_by %in% c("VirusPos","SIVsum") && !("VirusPos" %in% names(cell_meta))) {
    cell_meta$VirusPos <- ifelse(is.na(cell_meta$SIVsum), NA_character_,
                                 ifelse(cell_meta$SIVsum > 0, "SIVpos", "SIVneg"))
  }
  
  # --- COI selection ---
  if (is.null(coi)) {
    coi <- cell_meta %>% arrange(dplyr::desc(SIVsum)) %>% slice(1) %>% pull(cell)
  } else if (!(coi %in% cell_meta$cell)) {
    stop("Provided 'coi' not found in cell_meta$cell.")
  }
  
  # --- Restrict to same sample (optional) ---
  meta_nn <- if (restrict_to_sample) {
    if (!("SampleID" %in% names(cell_meta))) stop("restrict_to_sample=TRUE but 'SampleID' not in cell_meta.")
    coi_sample <- cell_meta %>% filter(cell == coi) %>% pull(SampleID)
    cell_meta %>% filter(SampleID == coi_sample)
  } else cell_meta
  
  # --- KNN neighborhood ---
  xy <- meta_nn %>% select(CenterX_global_px, CenterY_global_px) %>% as.matrix()
  rownames(xy) <- meta_nn$cell
  if (!(coi %in% rownames(xy))) stop("COI not found in neighborhood search set (check SampleID restriction).")
  k_eff <- min(k, nrow(xy))
  nn <- FNN::get.knnx(data = xy, query = xy, k = k_eff)
  rownames(nn$nn.index) <- rownames(xy)
  neighboring_cells <- rownames(nn$nn.index)[nn$nn.index[coi,]]
  neighboring_cells <- unique(c(neighboring_cells, coi))
  
  # --- Prepare plotting df ---
  cluster_present <- cluster_var %in% names(cell_meta)
  meta_cols <- c("cell","VirusPos","SIVsum","SampleID")
  if (cluster_present) meta_cols <- c(meta_cols, cluster_var)
  meta4join <- cell_meta %>% select(all_of(intersect(meta_cols, names(cell_meta))))
  
  plot.df <- cell_polygons %>%
    filter(cell %in% neighboring_cells) %>%
    left_join(meta4join, by = "cell")
  
  # --- Add gene expression if requested ---
  if (fill_by == "gene") {
    if (is.null(expr_matrix)) stop("fill_by='gene' requires expr_matrix.")
    if (is.null(gene)) stop("Please provide 'gene' when fill_by = 'gene'.")
    # exact match or case-insensitive fallback
    gcol <- if (gene %in% colnames(expr_matrix)) {
      gene
    } else {
      ci <- match(tolower(gene), tolower(colnames(expr_matrix)))
      if (is.na(ci)) stop("Gene '", gene, "' not found in expr_matrix (checked case-insensitively).")
      colnames(expr_matrix)[ci]
    }
    # align and pull values only for neighbors (faster)
    if (!all(neighboring_cells %in% rownames(expr_matrix))) {
      warning("Some neighboring cells missing from expr_matrix; gene values will be NA.")
    }
    gene_vec <- expr_matrix[neighboring_cells, gcol, drop = TRUE]
    gene_df <- data.frame(cell = neighboring_cells, .GeneExpr = as.numeric(gene_vec))
    plot.df <- plot.df %>% left_join(gene_df, by = "cell")
    
    # transform if requested
    if (gene_transform == "log1p") plot.df$.GeneExpr <- log1p(plot.df$.GeneExpr)
    if (gene_transform == "sqrt")  plot.df$.GeneExpr <- sqrt(plot.df$.GeneExpr)
  }
  
  # --- Scale bar ---
  scale_bar <- NULL
  if (!is.null(scale_bar_fun)) {
    scale_bar <- scale_bar_fun(x_vals = plot.df$x_global_px, y_vals = plot.df$y_global_px)
  }
  
  # --- Fill mapping ---
  fill_mapping <- switch(
    fill_by,
    "VirusPos" = aes(fill = .data[["VirusPos"]]),
    "cluster"  = {
      if (!cluster_present) warning("cluster_var not found; falling back to VirusPos.")
      aes(fill = .data[[if (cluster_present) cluster_var else "VirusPos"]])
    },
    "SIVsum"   = aes(fill = .data[["SIVsum"]]),
    "gene"     = aes(fill = .data[[".GeneExpr"]])
  )
  
  # --- Plot ---
  p <- ggplot(
    plot.df,
    aes(x = x_global_px, y = y_global_px, group = cell,
        alpha = (cell == coi), color = (cell == coi))
  ) +
    fill_mapping +
    geom_polygon(linewidth = linewidth) +
    scale_alpha_manual(values = c("TRUE" = alpha_coi, "FALSE" = alpha_bg)) +
    scale_color_manual(values = c("TRUE" = outline_coi, "FALSE" = outline_bg)) +
    coord_fixed() +
    theme_void() +
    guides(alpha = "none", color = "none")
  
  if (fill_by %in% c("SIVsum","gene")) {
    p <- p + scale_fill_viridis_c(option = "C") +
      labs(fill = if (fill_by == "gene") paste0(gene, if (gene_transform!="none") paste0(" (", gene_transform, ")") else "")
           else "SIVsum")
  } else {
    p <- p + labs(fill = if (fill_by == "cluster") "Cell type" else "Viral status")
  }
  
  if (!is.null(scale_bar)) {
    p <- p + scale_bar$bg + scale_bar$rect + scale_bar$label
  }
  
  attr(p, "coi") <- coi
  attr(p, "neighbors") <- neighboring_cells
  if (fill_by == "gene") attr(p, "gene") <- gene
  return(p)
}









# ---- ROI helpers ------------------------------------------------------

# Select non-overlapping hotspot centers among SIV_Plausible cells
#' Select Non-Overlapping Hotspot Seed Centers
#'
#' Identifies spatially separated seed points (typically SIV-plausible cells)
#' used to define regions of interest (ROIs). The function ranks candidate cells
#' by a scoring column (e.g., infection score), then iteratively selects the
#' highest-scoring cells while enforcing a minimum spatial separation between
#' chosen seeds.
#'
#' @param meta A data frame of per-cell metadata containing spatial coordinates
#'   and scoring information. Must include:
#'   \describe{
#'     \item{\code{x_slide_mm}}{Numeric x-coordinate (in mm).}
#'     \item{\code{y_slide_mm}}{Numeric y-coordinate (in mm).}
#'   }
#' @param score_col Character string specifying the column name containing the
#'   numeric score used to rank potential seeds (e.g., \code{"SIV_UCell"}).
#' @param flag_col Character string specifying the logical column indicating
#'   which cells are eligible to be considered as seeds (e.g., \code{"SIV_Plausible"}).
#' @param r_mm Numeric. Base spatial scale (in mm) used for ROI construction; also
#'   serves as the default for \code{min_sep} if not explicitly provided.
#' @param max_seeds Integer. Maximum number of seeds to select (default: \code{Inf}).
#' @param min_sep Numeric. Minimum allowed Euclidean distance (in mm) between
#'   two selected seeds (default: \code{r_mm}).
#' @param verbose Logical. If \code{TRUE}, prints per-iteration diagnostic output
#'   showing seed coordinates and current distance thresholds (default: \code{FALSE}).
#'
#' @return A data frame of selected seed cells (subset of \code{meta}) containing
#'   all original metadata columns for each chosen seed. The number of returned rows
#'   corresponds to the number of independent hotspots identified.
#'
#' @details
#' The algorithm proceeds as follows:
#' \enumerate{
#'   \item Filter candidate cells where \code{flag_col == TRUE}.
#'   \item Sort candidates in descending order of \code{score_col}.
#'   \item Iteratively select the highest-scoring cell, then remove all remaining
#'         candidates within \code{min_sep} of that seed.
#' }
#'
#' This ensures a spatially thinned set of non-overlapping hotspots, suitable for
#' defining independent ROI centers.
#'
#' The function uses a soft exclusion criterion (\code{0.9 * min_sep}) to prevent
#' over-thinning in dense areas, balancing independence with coverage.
#'
#' @examples
#' meta <- data.frame(
#'   x_slide_mm = runif(1000, 0, 2),
#'   y_slide_mm = runif(1000, 0, 2),
#'   SIV_UCell = runif(1000),
#'   SIV_Plausible = sample(c(TRUE, FALSE), 1000, TRUE)
#' )
#' seeds <- select_hotspot_seeds(meta, score_col = "SIV_UCell",
#'                               flag_col = "SIV_Plausible",
#'                               r_mm = 0.05)
#' head(seeds)
#'
#' @seealso
#' \code{\link{build_rois_fast}} and \code{\link{make_control_centers}} for
#' downstream ROI generation, and \code{\link{summarize_rois}} for ROI summarization.
#'
#' @export
select_hotspot_seeds <- function(meta, score_col, flag_col, r_mm, 
                                 max_seeds = Inf, min_sep = NULL, verbose = FALSE) {
  if (is.null(min_sep)) min_sep <- r_mm
  stopifnot(all(c("x_slide_mm", "y_slide_mm") %in% names(meta)))
  
  df <- meta
  df$idx <- seq_len(nrow(df))
  df <- df[df[[flag_col]] %in% TRUE, , drop = FALSE]
  if (!nrow(df)) return(df[0, ])
  
  df$x_slide_mm <- as.numeric(df$x_slide_mm)
  df$y_slide_mm <- as.numeric(df$y_slide_mm)
  df <- df[order(df[[score_col]], decreasing = TRUE), ]
  picked_idx <- integer(0)
  
  while (nrow(df) > 0 && length(picked_idx) < max_seeds) {
    seed <- df[1, ]
    picked_idx <- c(picked_idx, seed$idx)
    
    x_num <- df$x_slide_mm
    y_num <- df$y_slide_mm
    sx <- seed$x_slide_mm
    sy <- seed$y_slide_mm
    
    dx <- x_num - sx
    dy <- y_num - sy
    dist2 <- dx * dx + dy * dy
    
    if (verbose)
      cat("Iter:", length(picked_idx), 
          "seed:", round(sx, 3), round(sy, 3),
          "min dist² =", min(dist2), "thresh² =", min_sep^2, "\n")
    
    # softer exclusion (previous working behavior)
    keep <- dist2 > ((min_sep * 0.9)^2)
    df <- df[keep, , drop = FALSE]
  }
  
  out <- meta[seq_len(nrow(meta)) %in% picked_idx, , drop = FALSE]
  message("Picked ", nrow(out), " seeds from ", sum(meta[[flag_col]], na.rm = TRUE), 
          " plausible cells (min_sep = ", min_sep, ").")
  out
}

# Build ROIs (one circle per center), returns long table (one row per cell in an ROI)
#' Build Regions of Interest (ROIs) Around Specified Centers
#'
#' Constructs circular regions of interest (ROIs) by selecting all cells within a
#' specified radius (\code{r_mm}) around given center coordinates. Each ROI is labeled
#' and returned as a long-format data frame suitable for downstream summarization
#' and visualization.
#'
#' @param meta A data frame of per-cell metadata containing spatial coordinates and
#'   any additional cell-level annotations. Must include:
#'   \describe{
#'     \item{\code{x_slide_mm}}{Numeric x-coordinate (in mm).}
#'     \item{\code{y_slide_mm}}{Numeric y-coordinate (in mm).}
#'   }
#' @param centers A data frame of ROI center coordinates with columns
#'   \code{x_slide_mm} and \code{y_slide_mm}. Each row represents one ROI center.
#' @param r_mm Numeric. Radius (in millimeters) defining the circular region around
#'   each center to include cells (default: \code{0.05}).
#' @param roi_prefix Character string defining a prefix for ROI identifiers and types
#'   (e.g., \code{"case"} or \code{"control"}; default: \code{"case"}).
#'
#' @return A data frame where each row corresponds to a single cell included in one ROI.
#'   The output retains all columns from \code{meta} and adds:
#'   \describe{
#'     \item{\code{roi_id}}{Unique ROI identifier (e.g., \code{"case_1"}).}
#'     \item{\code{roi_cx}, \code{roi_cy}}{ROI center coordinates.}
#'     \item{\code{roi_radius_mm}}{Numeric radius used to define the ROI.}
#'     \item{\code{roi_type}}{ROI label defined by \code{roi_prefix}.}
#'   }
#'
#' @details
#' This function performs a straightforward per-center Euclidean distance check to
#' include all cells within a circular region of radius \code{r_mm}. While slower
#' than \code{\link{build_rois_fast}} for large datasets, it is easier to interpret
#' and useful for validation or small-scale analyses.
#'
#' If no cells fall within the ROI for a given center, that ROI is skipped.
#' An empty data frame is returned if no centers produce valid ROIs.
#'
#' @examples
#' meta <- data.frame(
#'   x_slide_mm = runif(1000, 0, 2),
#'   y_slide_mm = runif(1000, 0, 2)
#' )
#' centers <- data.frame(
#'   x_slide_mm = c(0.5, 1.5),
#'   y_slide_mm = c(0.5, 1.5)
#' )
#' rois <- build_rois(meta, centers, r_mm = 0.1, roi_prefix = "case")
#' head(rois)
#'
#' @seealso
#' \code{\link{build_rois_fast}} for a faster nearest-neighbor version, and
#' \code{\link{summarize_rois}} for ROI-level aggregation.
#'
#' @importFrom dplyr bind_rows
#' @export
build_rois <- function(meta,
                       centers,          # data.frame with x_slide_mm, y_slide_mm, and an id
                       r_mm = .05,
                       roi_prefix = "case") {
  stopifnot(all(c("x_slide_mm","y_slide_mm") %in% names(meta)))
  if (!nrow(centers)) return(data.frame()[0, ])
  
  out_ls <- vector("list", nrow(centers))
  for (i in seq_len(nrow(centers))) {
    cx <- centers$x_slide_mm[i]
    cy <- centers$y_slide_mm[i]
    rid <- paste0(roi_prefix, "_", i)
    
    dx <- meta$x_slide_mm - cx
    dy <- meta$y_slide_mm - cy
    in_roi <- (dx*dx + dy*dy) <= (r_mm^2)
    
    if (any(in_roi)) {
      tmp <- meta[in_roi, , drop = FALSE]
      tmp$roi_id <- rid
      tmp$roi_cx <- cx
      tmp$roi_cy <- cy
      tmp$roi_radius_mm <- r_mm
      tmp$roi_type <- roi_prefix
      out_ls[[i]] <- tmp
    }
  }
  dplyr::bind_rows(out_ls)
}


#' Build ROIs Around Seed Centers (Fast Approximation)
#'
#' Efficiently constructs regions of interest (ROIs) by identifying all cells within
#' a specified radius (\code{r_mm}) around each provided seed (center) coordinate.
#' This implementation uses fast nearest-neighbor search via the \pkg{FNN} package
#' for scalable performance on large datasets.
#'
#' @param meta A data frame of per-cell metadata containing spatial coordinates and
#'   other cell-level annotations. Must include:
#'   \describe{
#'     \item{\code{x_slide_mm}}{Numeric x-coordinate (in mm).}
#'     \item{\code{y_slide_mm}}{Numeric y-coordinate (in mm).}
#'   }
#' @param seed_df A data frame of ROI seed centers, with columns
#'   \code{x_slide_mm} and \code{y_slide_mm}, specifying the ROI centers to be expanded.
#' @param r_mm Numeric. ROI radius in millimeters used to define the circular
#'   neighborhood for each ROI (default: \code{0.05}).
#' @param label Character string defining the ROI type label to assign
#'   (e.g., \code{"case"} or \code{"control"}; default: \code{"case"}).
#'
#' @return A data frame in “long” format containing all cells that fall within at least
#'   one ROI. Each row corresponds to a single cell and includes:
#'   \describe{
#'     \item{\code{roi_id}}{Unique ROI identifier (e.g., \code{"case_1"}).}
#'     \item{\code{roi_centerx}, \code{roi_centery}}{Coordinates of the ROI center.}
#'     \item{\code{roi_type}}{ROI label (e.g., \code{"case"} or \code{"control"}).}
#'     \item{\code{roi_radiusmm}}{Numeric radius used to define the ROI.}
#'   }
#'   The output also retains all columns from the input \code{meta} data frame
#'   corresponding to cells included within an ROI.
#'
#' @details
#' This function leverages \pkg{FNN}’s \code{\link[FNN]{get.knnx}} to efficiently
#' compute nearest-neighbor distances between ROI centers and all cell coordinates.
#' Each cell within \code{r_mm} of a seed center is included in that ROI.
#'
#' It provides a significant speed-up compared to pure R distance loops, particularly
#' for large spatial datasets typical of single-cell or spatial transcriptomic analyses.
#'
#' A simple diagnostic scatter plot is generated by default, showing all cell
#' coordinates (gray) and ROI seed centers (red).
#'
#' @examples
#' meta <- data.frame(
#'   x_slide_mm = runif(1000, 0, 2),
#'   y_slide_mm = runif(1000, 0, 2)
#' )
#' seeds <- data.frame(
#'   x_slide_mm = c(0.5, 1.5),
#'   y_slide_mm = c(0.5, 1.5)
#' )
#' rois <- build_rois_fast(meta, seed_df = seeds, r_mm = 0.1, label = "case")
#' head(rois)
#'
#' @seealso
#' \code{\link{make_control_centers}}, \code{\link{select_hotspot_seeds}},
#' and \code{\link{summarize_rois}} for related ROI generation and analysis steps.
#'
#' @importFrom FNN get.knnx
#' @export
build_rois_fast <- function(meta,
                            seed_df,
                            r_mm = .05,
                            label = "case") {
  # seed_df = seeds
  stopifnot(requireNamespace("FNN", quietly = TRUE))
  if (!nrow(seed_df)) return(meta[0, , drop = FALSE])
  
  # Drop NA coords
  meta <- meta[!is.na(meta$x_slide_mm) & !is.na(meta$y_slide_mm), , drop = FALSE]
  seed_df <- seed_df[!is.na(seed_df$x_slide_mm) & !is.na(seed_df$y_slide_mm), , drop = FALSE]
  if (!nrow(meta) || !nrow(seed_df)) return(meta[0, , drop = FALSE])
  
  coords_all  <- as.matrix(meta[, c("x_slide_mm", "y_slide_mm")])
  coords_seed <- as.matrix(seed_df[, c("x_slide_mm", "y_slide_mm")])
  
  plot(coords_all, pch=20, col = "grey")
  points(coords_seed, pch=20, col = "red")
  
  nn <- FNN::get.knnx(data = coords_all, query = coords_seed, k = nrow(coords_all))
  
  # nn$nn.dist[1:10,1:10]
  # dim(nn$nn.dist)
  # 
  # print(summary(as.vector(nn$nn.dist)))
  
  roi_list <- lapply(seq_len(nrow(coords_seed)), function(i) {
    dists <- nn$nn.dist[i, ]
    keep_idx <- which(dists <= r_mm)
    if (!length(keep_idx)) return(NULL)
    
    roi_cells <- meta[keep_idx, , drop = FALSE]
    # Emit the exact column names downstream code expects
    roi_cells$roi_id       <- paste0(label, "_", i)
    roi_cells$roi_centerx  <- coords_seed[i, 1]
    roi_cells$roi_centery  <- coords_seed[i, 2]
    roi_cells$roi_type     <- label
    roi_cells$roi_radiusmm <- r_mm
    roi_cells
  })
  
  do.call(rbind, roi_list)
}


# Find control centers matched on density and far from plausible positives
# Returns a centers table similar to 'centers'
#' Generate Matched Control ROI Centers
#'
#' Identifies control ROI centers that are spatially distant from SIV-plausible cells
#' but matched in local cell density to corresponding case ROIs. Each case ROI center
#' is used to select one or more nearby control centers that have similar neighborhood
#' density, ensuring comparable sampling for downstream case–control comparisons.
#'
#' @param meta A data frame of per-cell metadata with:
#'   \describe{
#'     \item{x_slide_mm}{Numeric x-coordinate (in mm).}
#'     \item{y_slide_mm}{Numeric y-coordinate (in mm).}
#'   }
#' @param case_centers A data frame of case ROI centers with columns:
#'   \describe{
#'     \item{x_slide_mm}{X-coordinate of ROI center (mm).}
#'     \item{y_slide_mm}{Y-coordinate of ROI center (mm).}
#'   }
#' @param r_mm Numeric. ROI radius in millimeters; defines the neighborhood
#'   used for density comparison and exclusion of nearby SIV-plausible regions.
#' @param n_controls Integer. Number of matched control centers to generate per case
#'   center (default: 2).
#' @param density_tol Numeric fraction (default: 0.20). Tolerance for acceptable
#'   differences in local cell density between case and control ROIs. For example,
#'   \code{0.20} allows ±20% variation in the number of cells within \code{r_mm}.
#' @param plausible_col Character string specifying the column name in \code{meta}
#'   that indicates SIV-plausible cells (default: \code{"SIV_Plausible"}).
#' @param attempts_per_case Integer (default: 500). Maximum number of random attempts
#'   to find suitable control centers per case.
#' @param seed Integer random seed (default: 2025) for reproducibility.
#' @param min_sep Numeric. Minimum separation (in mm) between control centers to
#'   prevent spatial overlap (default: \code{r_mm}).
#'
#' @return A data frame of control ROI centers with:
#'   \describe{
#'     \item{x_slide_mm, y_slide_mm}{Coordinates of each control ROI center.}
#'     \item{case_id}{ID of the matched case ROI.}
#'     \item{control_rank}{Numeric rank of control (1..n_controls).}
#'     \item{control_id}{Unique control ROI identifier.}
#'   }
#'
#' @details
#' The algorithm randomly samples candidate positions from cells that are not
#' marked as SIV-plausible and are farther than \code{r_mm} from any plausible
#' cell. For each case center, it iteratively selects control centers whose local
#' cell density (number of cells within \code{r_mm}) falls within the tolerance
#' range of the corresponding case ROI. Candidates too close to an existing control
#' are excluded using the \code{min_sep} criterion.
#'
#' If no valid control centers are found for a given case, the result will return
#' an empty data frame.
#'
#' @examples
#' meta <- data.frame(
#'   x_slide_mm = runif(500, 0, 2),
#'   y_slide_mm = runif(500, 0, 2),
#'   SIV_Plausible = sample(c(TRUE, FALSE), 500, TRUE)
#' )
#' case_centers <- data.frame(
#'   x_slide_mm = runif(5, 0.2, 1.8),
#'   y_slide_mm = runif(5, 0.2, 1.8)
#' )
#' ctrs <- make_control_centers(meta, case_centers, r_mm = 0.1, n_controls = 2)
#' head(ctrs)
#'
#' @seealso
#' \code{\link{select_hotspot_seeds}}, \code{\link{build_rois_fast}},
#' and \code{\link{summarize_rois}} for related ROI generation and summarization steps.
#'
#' @export
make_control_centers <- function(meta, case_centers, r_mm = 0.05, n_controls = 2,
                                 density_tol = 0.20, plausible_col = "SIV_Plausible",
                                 attempts_per_case = 500, seed = 2025, min_sep = 0.05) {
  stopifnot(all(c("x_slide_mm","y_slide_mm") %in% names(meta)))
  set.seed(seed)
  
  plausible <- meta[[plausible_col]] %in% TRUE
  if (any(plausible)) {
    px <- meta$x_slide_mm[plausible]
    py <- meta$y_slide_mm[plausible]
    mask_far <- rep(TRUE, nrow(meta))
    for (j in seq_along(px)) {
      dx <- meta$x_slide_mm - px[j]
      dy <- meta$y_slide_mm - py[j]
      mask_far[(dx*dx + dy*dy) <= (r_mm^2)] <- FALSE
    }
    cand <- meta[mask_far, , drop = FALSE]
  } else {
    cand <- meta
  }
  
  out <- vector("list", nrow(case_centers))
  for (i in seq_len(nrow(case_centers))) {
    cx <- case_centers$x_slide_mm[i]
    cy <- case_centers$y_slide_mm[i]
    
    n_case <- sum((meta$x_slide_mm - cx)^2 + (meta$y_slide_mm - cy)^2 <= (r_mm^2))
    target_min <- floor((1 - density_tol) * n_case)
    target_max <- ceiling((1 + density_tol) * n_case)
    
    picked <- list()
    attempts <- 0
    while (length(picked) < n_controls && attempts < attempts_per_case && nrow(cand) > 0) {
      attempts <- attempts + 1
      j <- sample.int(nrow(cand), 1)
      tx <- cand$x_slide_mm[j]
      ty <- cand$y_slide_mm[j]
      
      n_here <- sum((meta$x_slide_mm - tx)^2 + (meta$y_slide_mm - ty)^2 <= (r_mm^2))
      if (n_here >= target_min && n_here <= target_max) {
        picked[[length(picked) + 1]] <- cand[j, c("x_slide_mm","y_slide_mm")]
        cand <- cand[(cand$x_slide_mm - tx)^2 + (cand$y_slide_mm - ty)^2 > (min_sep^2), , drop = FALSE]
      }
    }
    
    if (length(picked)) {
      cc <- dplyr::bind_rows(picked)
      cc$case_id <- paste0("case_", i)
      cc$control_rank <- seq_len(nrow(cc))
      out[[i]] <- cc
    }
  }
  
  ctrs <- dplyr::bind_rows(out)
  if (!nrow(ctrs)) return(ctrs)
  ctrs$control_id <- paste0(ctrs$case_id, "_ctrl_", ctrs$control_rank)
  ctrs
}


# Summarize each ROI to one row (fractions and counts)
#' Summarize ROI-Level Cell Composition and SIV Metrics
#'
#' Aggregates per-cell metadata within each region of interest (ROI) to produce
#' ROI-level summaries of cell counts, infection-related fractions, and cell-type
#' composition. This function serves as the core summarization step in the
#' \pkg{CellNicheR} spatial ROI analysis pipeline.
#'
#' @param roi_long A data frame containing all cells assigned to one or more ROIs.
#'   Must include at least the columns:
#'   \describe{
#'     \item{\code{roi_id}}{Unique ROI identifier.}
#'     \item{\code{roi_type}}{ROI label (e.g., \code{"case"} or \code{"control"}).}
#'     \item{\code{SampleID}}{Sample or section identifier.}
#'     \item{\code{x_slide_mm}, \code{y_slide_mm}}{Optional; per-cell spatial coordinates.}
#'   }
#'   Additional columns specified by the arguments below must also be present.
#'
#' @param celltype_col Character string giving the name of the column in \code{roi_long}
#'   that defines cell types or clusters.
#' @param positive_col Character string specifying the column indicating SIV-positive
#'   cells (default: \code{"SIV_Positive"}).
#' @param plausible_col Character string specifying the column indicating
#'   SIV-plausible (susceptible) cells (default: \code{"SIV_Plausible"}).
#' @param score_col Character string specifying the column containing per-cell
#'   SIV activity or infection scores (default: \code{"SIV_UCell"}).
#'
#' @return A data frame with one row per ROI, containing:
#'   \describe{
#'     \item{\code{roi_id}, \code{roi_type}, \code{SampleID}}{ROI identifiers.}
#'     \item{\code{n_cells}}{Number of cells within the ROI.}
#'     \item{\code{frac_pos}}{Fraction of cells that are SIV-positive.}
#'     \item{\code{frac_plaus}}{Fraction of cells that are SIV-plausible.}
#'     \item{\code{mean_score}, \code{median_score}}{Mean and median SIV score per ROI.}
#'     \item{Cell-type columns}{Additional columns representing the proportion
#'       of each cell type (from \code{celltype_col}) within the ROI, as wide-format
#'       fractional values that sum to 1 per ROI.}
#'   }
#'
#' @details
#' The function first computes infection-related summaries for each ROI using
#' the specified \code{positive_col}, \code{plausible_col}, and \code{score_col}.
#' It then tabulates per-ROI cell-type proportions and merges these into a single
#' summary table suitable for downstream visualization or statistical analysis.
#'
#' @examples
#' roi_long <- data.frame(
#'   roi_id = rep(paste0("case_", 1:3), each = 10),
#'   roi_type = rep("case", 30),
#'   SampleID = "LN1",
#'   SIV_Positive = rbinom(30, 1, 0.2),
#'   SIV_Plausible = rbinom(30, 1, 0.5),
#'   SIV_UCell = runif(30),
#'   CellType = sample(c("T.cell", "B.cell", "Macrophage"), 30, TRUE)
#' )
#' summarize_rois(roi_long, celltype_col = "CellType")
#'
#' @seealso
#' \code{\link{run_roi_pipeline_one}}, \code{\link{build_rois_fast}},
#' and \code{\link{compute_overlap_stats}} for related ROI processing functions.
#'
#' @export
summarize_rois <- function(roi_long,
                           celltype_col,
                           positive_col = "SIV_Positive",
                           plausible_col = "SIV_Plausible",
                           score_col = "SIV_UCell") {
  
  stopifnot(all(c("roi_id","roi_type") %in% names(roi_long)))
  ct <- celltype_col
  
  roi_long %>%
    dplyr::group_by(roi_id, roi_type, SampleID, .add = TRUE) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      frac_pos = mean(.data[[positive_col]] %in% TRUE, na.rm = TRUE),
      frac_plaus = mean(.data[[plausible_col]] %in% TRUE, na.rm = TRUE),
      mean_score = mean(.data[[score_col]], na.rm = TRUE),
      median_score = median(.data[[score_col]], na.rm = TRUE),
      .groups = "drop"
    ) -> base
  
  # cell-type composition wide (percentage per type)
  comp <- roi_long %>%
    dplyr::count(roi_id, .data[[ct]], name = "n") %>%
    dplyr::group_by(roi_id) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    tidyr::pivot_wider(id_cols = roi_id,
                       names_from = .data[[ct]],
                       values_from = prop,
                       values_fill = 0)
  
  dplyr::left_join(base, comp, by = "roi_id")
}


#' Compute Overlap Statistics Between ROI Centers
#'
#' Calculates pairwise spatial overlap metrics between regions of interest (ROIs),
#' given their center coordinates and a specified ROI radius. The function
#' quantifies how often circular ROIs overlap, both overall and separated by ROI type
#' (e.g., case vs. control). It is useful for assessing redundancy or spatial crowding
#' among ROI placements in \pkg{CellNicheR} analyses.
#'
#' @param centers A data frame containing ROI metadata with at least the columns:
#'   \describe{
#'     \item{\code{roi_centerx}}{Numeric x-coordinate of each ROI center (in mm).}
#'     \item{\code{roi_centery}}{Numeric y-coordinate of each ROI center (in mm).}
#'     \item{\code{roi_type}}{Character label identifying ROI type, typically
#'       \code{"case"} or \code{"control"}.}
#'   }
#' @param r Numeric scalar. Radius (in the same units as coordinates, typically mm)
#'   used to define ROI circular footprints.
#'
#' @return A named list of overlap statistics:
#'   \describe{
#'     \item{\code{pairs_total}}{Total number of unique ROI pairs evaluated.}
#'     \item{\code{pairs_overlapping}}{Count of ROI pairs whose centers are closer
#'       than twice the radius (i.e., overlapping circles).}
#'     \item{\code{frac_pairs_overlapping}}{Fraction of all pairs that overlap.}
#'     \item{\code{frac_pairs_case_case}}{Fraction of overlapping pairs among
#'       case–case ROI pairs.}
#'     \item{\code{frac_pairs_ctrl_ctrl}}{Fraction of overlapping pairs among
#'       control–control ROI pairs.}
#'     \item{\code{frac_pairs_between}}{Fraction of overlapping pairs between
#'       case and control ROIs.}
#'     \item{\code{frac_circles_that_overlap_any}}{Proportion of individual ROIs
#'       that overlap at least one other ROI.}
#'   }
#'
#' @details
#' The function computes Euclidean pairwise distances between ROI centers and
#' identifies overlaps where the center-to-center distance is less than
#' \eqn{2 * r}. It uses only the upper triangle of the distance matrix to avoid
#' duplicate pair counting.
#'
#' @examples
#' centers <- data.frame(
#'   roi_centerx = runif(10, 0, 1),
#'   roi_centery = runif(10, 0, 1),
#'   roi_type = rep(c("case", "control"), each = 5)
#' )
#' compute_overlap_stats(centers, r = 0.1)
#'
#' @export
compute_overlap_stats <- function(centers, r){
  mat <- as.matrix(dist(centers[, c("roi_centerx","roi_centery")]))
  n <- nrow(mat)
  # upper triangle pairs only
  ut <- upper.tri(mat)
  overlap_pair <- mat[ut] < (2*r)
  
  # within/between by roi_type
  type <- centers$roi_type
  A <- outer(type, type, FUN = function(a,b) paste(a,b,sep="|"))
  within_case <- (A == "case|case") & ut
  within_ctrl <- (A == "control|control") & ut
  between    <- (A %in% c("case|control","control|case")) & ut
  
  list(
    pairs_total = sum(ut),
    pairs_overlapping = sum(overlap_pair),
    frac_pairs_overlapping = mean(overlap_pair),
    frac_pairs_case_case = mean(mat[within_case] < (2*r)),
    frac_pairs_ctrl_ctrl = mean(mat[within_ctrl] < (2*r)),
    frac_pairs_between   = mean(mat[between] < (2*r)),
    frac_circles_that_overlap_any = mean(apply(mat < (2*r) & !diag(n), 1, any))
  )
}