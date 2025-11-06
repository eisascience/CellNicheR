### Thresholding and metadata augmentation -------------------------
#' Add SIV infection flags based on UCell score threshold
#'
#' This function computes a percentile-based threshold on the SIV-associated UCell score
#' (or any continuous infection-related score) and annotates cells as SIV-positive
#' if they exceed that threshold. The cutoff is stored in the Seurat object's
#' metadata and also recorded in the `@misc` slot for reference.
#'
#' @param serobj A Seurat object containing per-cell metadata and the specified score column.
#' @param score_col Character. The metadata column name containing the numeric SIV-related score
#'   (default: `"SIV_UCell"`).
#' @param cutoff_quantile Numeric between 0 and 1. The quantile of the score distribution used
#'   as the positivity threshold (default: `0.99` for 99th percentile).
#'
#' @return A Seurat object with two additions:
#'   \itemize{
#'     \item A new logical column `SIV_Positive` in the metadata, indicating cells above threshold.
#'     \item The computed numeric threshold stored in `serobj@misc$SIV_threshold`.
#'   }
#'
#' @details
#' This function is designed for CosMx or other spatial transcriptomic datasets
#' where a custom SIV gene signature (e.g., gag, pol, env) is quantified per cell.
#' The percentile cutoff can be tuned depending on the noise characteristics of
#' each slide or experimental condition.
#'
#' @examples
#' serobj <- add_SIV_flags(serobj, score_col = "SIV_UCell", cutoff_quantile = 0.99)
#' head(serobj@meta.data$SIV_Positive)
#'
#' @export
add_SIV_flags <- function(serobj, score_col = "SIV_UCell", cutoff_quantile = 0.99) {
  thr <- quantile(serobj[[score_col]][, 1], cutoff_quantile, na.rm = TRUE)
  serobj$SIV_Positive <- serobj[[score_col]][, 1] >= thr
  serobj@misc$SIV_threshold <- thr
  message("SIV threshold (", cutoff_quantile * 100, "th percentile): ", round(thr, 4))
  return(serobj)
}

### 2. Define biologically plausible target cell types ----------------
#' Annotate biologically plausible SIV-positive cells
#'
#' This function filters the SIV-positive cells to those belonging to cell types
#' known or expected to support SIV infection (for example, CD4 T cells or myeloid cells).
#' It creates a new logical metadata column indicating which SIV-positive cells
#' are biologically plausible based on their assigned cell type.
#'
#' @param serobj A Seurat object containing per-cell metadata, including a logical
#'   `SIV_Positive` column and a cell-type annotation column.
#' @param celltype_col Character. The metadata column name containing cell type labels
#'   (for example, `"RNA_Pre.Napari.v3.0.071825_Cell.Typing.InSituType.2_1_clusters"`).
#' @param susceptible_types Character vector of cell types considered biologically
#'   plausible SIV targets. Typical examples include `"T.cell.CD4"`, `"Macrophage"`,
#'   `"Monocyte"`, `"Conventional.dendritic.cell"`, etc.
#'
#' @return A Seurat object with an added logical column `SIV_Plausible` in the metadata.
#'   Cells are marked `TRUE` if they are SIV-positive *and* belong to one of the
#'   specified susceptible cell types.
#'
#' @details
#' This filtering step is useful for distinguishing likely true SIV infections from
#' false positives caused by ambient RNA, spillover, or technical artifacts.
#' The resulting flag can be used for focused spatial analysis and infection
#' hotspot quantification.
#'
#' @examples
#' susceptible <- c("T.cell.CD4", "T.cell.regulatory", "Macrophage", "Monocyte")
#' serobj <- add_SIV_plausible(serobj,
#'                             celltype_col = "RNA_Pre.Napari.v3.0.071825_Cell.Typing.InSituType.2_1_clusters",
#'                             susceptible_types = susceptible)
#' table(serobj$SIV_Plausible)
#'
#' @export
add_SIV_plausible <- function(serobj, celltype_col, susceptible_types) {
  serobj$SIV_Plausible <- serobj$SIV_Positive &
    serobj[[celltype_col]][, 1] %in% susceptible_types
  return(serobj)
}

### 3. Spatial plotting ------------------------------------------------
#' Plot spatial distribution of SIV-associated score
#'
#' Generates a spatial scatter plot showing the per-cell distribution of an
#' SIV-related continuous score (for example, UCell score) across the tissue.
#' Each point represents a single cell positioned by its CosMx slide coordinates
#' and colored according to the specified score column.
#'
#' @param serobj A Seurat object containing spatial metadata with `x_slide_mm` and
#'   `y_slide_mm` coordinates, as well as the specified score column in the metadata.
#' @param color_col Character. The metadata column containing the continuous SIV-related
#'   score to visualize (default: `"SIV_UCell"`).
#' @param palette Character vector of colors used for the gradient scale
#'   (default: `c("white", "grey", "lightblue", "gold", "red")`).
#' @param title Character. Plot title (default: `"SIV_UCell spatial distribution"`).
#'
#' @return A `ggplot2` object showing spatially resolved expression or infection
#'   score values across all cells in the sample.
#'
#' @details
#' This plot provides an overview of potential viral hotspots and regions of
#' elevated SIV-associated transcriptional signal. The color palette can be
#' adjusted to emphasize low or high expression regions.
#'
#' @examples
#' p_spatial <- plot_spatial_SIV(serobj,
#'                               color_col = "SIV_UCell",
#'                               palette = c("white", "grey", "lightblue", "gold", "red"),
#'                               title = "SIV_UCell spatial distribution")
#' print(p_spatial)
#'
#' @export
plot_spatial_SIV <- function(serobj, color_col = "SIV_UCell",
                             palette = c("white", "grey", "lightblue", "gold", "red"),
                             title = "SIV_UCell spatial distribution") {
  ggplot(serobj@meta.data, aes(x = x_slide_mm, y = y_slide_mm, color = !!sym(color_col))) +
    geom_point(size = 0.3) +
    scale_color_gradientn(colours = palette, name = color_col) +
    coord_fixed() +
    theme_void(base_size = 14) +
    ggtitle(title)
}

#' Plot binary spatial map of SIV-positive cells
#'
#' Creates a spatial scatter plot distinguishing SIV-positive and SIV-negative
#' cells across the tissue based on a logical metadata column. Each point
#' represents a single cell positioned using CosMx slide coordinates, colored
#' red for positive and grey for negative.
#'
#' @param serobj A Seurat object containing spatial metadata with `x_slide_mm` and
#'   `y_slide_mm` coordinates, and a logical column indicating positive cells.
#' @param flag_col Character. The metadata column name indicating SIV-positive
#'   cells (default: `"SIV_Positive"`).
#' @param legend_title Character. Legend title for the positive population
#'   (default: `"SIV+"`).
#'
#' @return A `ggplot2` object showing the spatial distribution of SIV-positive
#'   cells (red) overlaid on all other cells (grey). The total number of
#'   positive cells (`N = ...`) is displayed in the plot subtitle.
#'
#' @details
#' This plot highlights the spatial arrangement of cells identified as SIV-positive
#' by thresholding or scoring methods (e.g., via `add_SIV_flags`). It can be used
#' to visually assess clustering or hotspot localization of infected cells.
#'
#' @examples
#' p_bin <- plot_spatial_binary(serobj,
#'                              flag_col = "SIV_Positive",
#'                              legend_title = "SIV+")
#' print(p_bin)
#'
#' @export
plot_spatial_binary <- function(serobj,
                                flag_col = "SIV_Positive",
                                legend_title = "SIV+") {
  meta_df <- serobj@meta.data
  meta_df <- meta_df[order(meta_df[[flag_col]]), ]
  
  # Count total TRUE cells
  n_pos <- sum(meta_df[[flag_col]], na.rm = TRUE)
  
  ggplot(meta_df, aes(x = x_slide_mm, y = y_slide_mm, color = !!sym(flag_col))) +
    geom_point(size = 0.3) +
    scale_color_manual(values = c("grey80", "red"), name = legend_title) +
    coord_fixed() +
    theme_void(base_size = 14) +
    ggtitle(
      paste0(legend_title, " vs others"),
      subtitle = paste("N =", formatC(n_pos, format = "d", big.mark = ","),
                       "cells")
    )
}

#' Plot spatial distribution of biologically plausible SIV-positive cells
#'
#' Generates a spatial scatter plot showing the subset of SIV-positive cells that are
#' biologically plausible infection targets (for example, CD4 T cells, macrophages,
#' or dendritic cells). Each point represents a cell, positioned by CosMx slide
#' coordinates, colored red for plausible positives and grey for others.
#'
#' @param serobj A Seurat object containing spatial metadata with `x_slide_mm` and
#'   `y_slide_mm` coordinates, and a logical column `SIV_Plausible` indicating
#'   biologically plausible SIV-positive cells.
#'
#' @return A `ggplot2` object showing plausible SIV-positive cells (red) overlaid on
#'   all other cells (grey). The total count of plausible positives (`N = ...`) is
#'   displayed in the plot subtitle.
#'
#' @details
#' This visualization focuses on SIV-positive cells that belong to cell types known
#' to be permissive to infection. It is typically used after calling
#' `add_SIV_plausible()` to refine the SIV-positive subset based on biological
#' plausibility. It can reveal tissue regions enriched in likely infected cell types.
#'
#' @examples
#' p_plaus <- plot_spatial_plausible(serobj)
#' print(p_plaus)
#'
#' @export
plot_spatial_plausible <- function(serobj) {
  meta_df <- serobj@meta.data
  meta_df <- meta_df[order(meta_df$SIV_Plausible), ]
  
  # Count total plausible positives
  n_plaus <- sum(meta_df$SIV_Plausible, na.rm = TRUE)
  
  ggplot(meta_df, aes(x = x_slide_mm, y = y_slide_mm, color = SIV_Plausible)) +
    geom_point(size = 0.3) +
    scale_color_manual(values = c("grey80", "red"), name = "Plausible SIV+") +
    coord_fixed() +
    theme_void(base_size = 14) +
    ggtitle(
      "Filtered SIV+ cells",
      subtitle = paste("N =", formatC(n_plaus, format = "d", big.mark = ","),
                       "cells")
    )
}

### 4. Histogram of SIV_UCell distribution ----------------------------
#' Plot histogram of SIV-associated score with cutoff annotation
#'
#' Creates a histogram showing the distribution of an SIV-related score (for example,
#' UCell score) across all cells, using the same color palette as the spatial plots.
#' The chosen threshold (usually the 95th or 99th percentile) is indicated by a dashed
#' vertical line and labeled directly on the plot.
#'
#' @param serobj A Seurat object containing the specified score column in its metadata.
#'   The threshold value is optionally retrieved from `serobj@misc$SIV_threshold`.
#' @param score_col Character. The metadata column containing the numeric SIV-related
#'   score to visualize (default: `"SIV_UCell"`).
#' @param palette Character vector specifying the color gradient used to fill histogram
#'   bars (default: `c("white", "grey", "lightblue", "gold", "red")`).
#' @param cutoff Numeric. Optional score threshold value to display as a dashed
#'   vertical line. If `NULL`, the function automatically uses the threshold stored
#'   in `serobj@misc$SIV_threshold`.
#'
#' @return A `ggplot2` object showing the histogram of the score distribution with
#'   the selected threshold annotated.
#'
#' @details
#' This plot is typically used alongside `plot_spatial_SIV()` to visualize how the
#' chosen threshold partitions the population into SIV-positive and SIV-negative cells.
#' The color palette matches that used in spatial plots to ensure visual consistency.
#'
#' @examples
#' p_hist <- plot_SIV_histogram(serobj,
#'                              score_col = "SIV_UCell",
#'                              palette = c("white", "grey", "lightblue", "gold", "red"),
#'                              cutoff = serobj@misc$SIV_threshold)
#' print(p_hist)
#'
#' @export
plot_SIV_histogram <- function(serobj, score_col = "SIV_UCell",
                               palette = c("white", "grey", "lightblue", "gold", "red"),
                               cutoff = NULL) {
  if (is.null(cutoff)) cutoff <- serobj@misc$SIV_threshold
  ggplot(serobj@meta.data, aes(x = !!sym(score_col))) +
    geom_histogram(aes(fill = ..x..), bins = 80, color = "black", size = 0.1) +
    scale_fill_gradientn(colours = palette, name = score_col) +
    geom_vline(xintercept = cutoff, color = "black",
               linetype = "dashed", linewidth = 0.8) +
    annotate("text", x = cutoff, y = Inf,
             label = paste0("cutoff (", round(cutoff, 3), ")"),
             vjust = 2, hjust = 1.1, size = 3.5) +
    theme_classic(base_size = 14) +
    xlab(score_col) + ylab("Cell count") +
    ggtitle(paste("Distribution of", score_col))
}

### 5. Cell-type summary plotting -------------------------------------
#' Summarize and visualize SIV-positive cells by cell type
#'
#' Generates a bar plot summarizing the fraction (or count) of SIV-positive cells
#' across annotated cell types. The plot can optionally display percentage labels
#' on each bar, showing the proportion of positive cells within each cell type.
#'
#' @param serobj A Seurat object containing per-cell metadata, including a cell type
#'   annotation column and a logical flag column (for example, `SIV_Positive` or
#'   `SIV_Plausible`).
#' @param celltype_col Character. The metadata column containing cell type labels.
#' @param flag_col Character. The metadata column indicating SIV-positive cells
#'   (default: `"SIV_Positive"`).
#' @param label_percent Logical. Whether to annotate bars with the percentage of
#'   positive cells (default: `TRUE`).
#' @param label_title Character. Legend title for the positive flag
#'   (default: `"SIV+"`).
#'
#' @return A `ggplot2` object showing the proportion of SIV-positive cells within
#'   each annotated cell type. Cell types are ordered by decreasing proportion
#'   of positives.
#'
#' @details
#' This function computes, for each cell type, the proportion of cells that are
#' positive according to the specified logical flag column. It is useful for
#' identifying which cell types are most enriched for viral signal. When used
#' with `SIV_Plausible`, the plot highlights biologically plausible infection
#' patterns across permissive cell types.
#'
#' @examples
#' # Summarize all SIV-positive cells by cell type
#' p_summary_all <- summarize_SIV_by_celltype(
#'   serobj,
#'   celltype_col = "RNA_Pre.Napari.v3.0.071825_Cell.Typing.InSituType.2_1_clusters",
#'   flag_col = "SIV_Positive"
#' )
#'
#' # Summarize only biologically plausible SIV-positive cells
#' p_summary_plaus <- summarize_SIV_by_celltype(
#'   serobj,
#'   celltype_col = "RNA_Pre.Napari.v3.0.071825_Cell.Typing.InSituType.2_1_clusters",
#'   flag_col = "SIV_Plausible",
#'   label_title = "Plausible SIV+"
#' )
#'
#' @export
summarize_SIV_by_celltype <- function(serobj, celltype_col,
                                      flag_col = "SIV_Positive",
                                      label_percent = TRUE,
                                      label_title = "SIV+") {
  
  # Extract the metadata for consistent length
  meta <- serobj@meta.data
  
  # Build the table
  df <- as.data.frame(table(
    CellType = meta[[celltype_col]],
    Flag = meta[[flag_col]]
  ))
  
  # Compute proportions
  df <- df %>%
    group_by(CellType) %>%
    mutate(Prop = Freq / sum(Freq))
  
  # Order by decreasing fraction of TRUE
  order_df <- df %>%
    filter(Flag == TRUE) %>%
    arrange(desc(Prop))
  
  df$CellType <- factor(df$CellType, levels = order_df$CellType)
  
  # Create label data frame
  df_label <- df %>%
    filter(Flag == TRUE) %>%
    mutate(percent = round(Prop * 100, 1))
  
  # Plot
  p <- ggplot(df, aes(x = CellType, y = Prop, fill = Flag)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
    scale_fill_manual(values = c("grey80", "red"), name = label_title) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    coord_flip() +
    labs(x = "Cell type", y = "Proportion of cells",
         title = paste(label_title, "cells by cell types")) +
    theme_classic(base_size = 14) +
    theme(axis.text.y = element_text(size = 10))
  
  if (label_percent) {
    p <- p + geom_text(
      data = df_label,
      aes(label = paste0(percent, "%")),
      position = position_stack(vjust = 0.5),
      size = 3
    )
  }
  
  return(p)
}




#' Generate combined spatial and summary visualization of SIV infection
#'
#' Creates a composite figure integrating spatial, histogram, and cell-type
#' summary plots for SIV-associated signals within a Seurat object. This function
#' applies the full analysis pipeline: thresholding based on SIV score,
#' filtering for biologically plausible infection targets, and visualizing
#' both spatial distributions and per-cell-type summaries.
#'
#' @param serobj A Seurat object containing spatial metadata (`x_slide_mm`,
#'   `y_slide_mm`) and an SIV-related score column (for example, `SIV_UCell`).
#' @param celltype_col Character. The metadata column specifying per-cell type
#'   annotations used to stratify SIV-positive counts.
#' @param susceptible_types Character vector of biologically plausible
#'   SIV-permissive cell types (for example, `"T.cell.CD4"`, `"Macrophage"`,
#'   `"Monocyte"`, `"Conventional.dendritic.cell"`).
#' @param cutoff_quantile Numeric between 0 and 1. Quantile used to determine the
#'   threshold for SIV positivity (default: `0.99` for 99th percentile).
#'
#' @return A `patchwork` composite object combining six panels:
#'   \itemize{
#'     \item Top row:
#'       \enumerate{
#'         \item Continuous SIV score spatial map (`plot_spatial_SIV`)
#'         \item Histogram of SIV score with threshold (`plot_SIV_histogram`)
#'         \item Binary SIV+ vs others map (`plot_spatial_binary`)
#'       }
#'     \item Bottom row:
#'       \enumerate{
#'         \item Spatial map of biologically plausible SIV+ cells (`plot_spatial_plausible`)
#'         \item Bar plot of SIV+ fraction by cell type (`summarize_SIV_by_celltype`)
#'         \item Bar plot of plausible SIV+ fraction by cell type (`summarize_SIV_by_celltype`)
#'       }
#'   }
#'
#' @details
#' This high-level function runs the full visualization workflow for a given
#' CosMx-derived Seurat object. It performs thresholding (via `add_SIV_flags`),
#' filters plausible infected cell types (`add_SIV_plausible`), and combines
#' outputs from all lower-level plotting functions into a single summary figure.
#' The output can be saved directly using `ggsave()` for slide- or sample-level reports.
#'
#' @examples
#' combined_plot <- run_SIV_visualization(
#'   serobj = SerObjLS$PerLN_34100_LN1,
#'   celltype_col = "RNA_Pre.Napari.v3.0.071825_Cell.Typing.InSituType.2_1_clusters",
#'   susceptible_types = c("T.cell.CD4", "T.cell.regulatory",
#'                         "Macrophage", "Monocyte",
#'                         "Conventional.dendritic.cell",
#'                         "Plasmacytoid.dendritic.cell"),
#'   cutoff_quantile = 0.99
#' )
#' print(combined_plot)
#'
#' @export
run_SIV_visualization <- function(serobj, celltype_col,
                                  susceptible_types,
                                  cutoff_quantile = 0.99) {
  
  serobj <- add_SIV_flags(serobj, "SIV_UCell", cutoff_quantile)
  serobj <- add_SIV_plausible(serobj, celltype_col, susceptible_types)
  
  p_spatial <- plot_spatial_SIV(serobj)
  p_hist <- plot_SIV_histogram(serobj) + NoLegend()
  p_binary <- plot_spatial_binary(serobj)
  p_plaus <- plot_spatial_plausible(serobj) + NoLegend()
  p_summary_all <- summarize_SIV_by_celltype(serobj, celltype_col, "SIV_Positive")
  p_summary_plaus <- summarize_SIV_by_celltype(serobj, celltype_col,
                                               "SIV_Plausible",
                                               label_title = "Plausible SIV+") + NoLegend()
  
  (p_spatial | p_hist | p_binary) / ( p_plaus | (p_summary_all / p_summary_plaus))
}







run_roi_pipeline_one <- function(serobj, celltype_col,
                                 r_mm = 0.05, max_seeds = Inf,
                                 n_controls = 2, density_tol = 0.20, min_sep = 0.05) {
  meta <- serobj@meta.data
  meta$x_slide_mm <- as.numeric(meta$x_slide_mm)
  meta$y_slide_mm <- as.numeric(meta$y_slide_mm)
  
  if (!("LN_class" %in% names(meta))) {
    meta$LN_class <- dplyr::case_when(
      grepl("Mes", meta$SampleID, ignore.case = TRUE) ~ "Mes",
      grepl("Per", meta$SampleID, ignore.case = TRUE) ~ "Per",
      TRUE ~ "Other"
    )
  }
  
  seeds <- select_hotspot_seeds(meta, "SIV_UCell", "SIV_Plausible",
                                r_mm = r_mm, max_seeds = max_seeds, min_sep = min_sep)
  if (!nrow(seeds)) {
    message("No plausible SIV+ seeds found; returning empty results.")
    return(list(roi_cells = meta[0, ], roi_summary = meta[0, ]))
  }
  
  message("Building ROIs around ", nrow(seeds), " SIV+ seeds...")
  rois_case <- build_rois_fast(meta, seeds, r_mm = r_mm, label = "case")
  message("Built ", length(unique(rois_case$roi_id)), " case ROIs.")
  
  ctrl_centers <- make_control_centers(meta, seeds[, c("x_slide_mm","y_slide_mm"), drop = FALSE],
                                       r_mm = r_mm, n_controls = n_controls,
                                       density_tol = density_tol, plausible_col = "SIV_Plausible")
  
  if (nrow(ctrl_centers)) {
    rois_ctrl <- build_rois_fast(meta, ctrl_centers, r_mm = r_mm, label = "control")
    message("Built ", length(unique(rois_ctrl$roi_id)), " control ROIs.")
  } else {
    message("No valid control centers found.")
    rois_ctrl <- meta[0, ]
  }
  
  roi_long <- dplyr::bind_rows(rois_case, rois_ctrl)
  names(roi_long) <- make.unique(names(roi_long), sep = "__dup")
  
  roi_sum <- summarize_rois(roi_long, celltype_col)
  message("ROI pipeline completed: ", nrow(roi_sum), " summarized entries.")
  list(roi_cells = roi_long, roi_summary = roi_sum)
}



