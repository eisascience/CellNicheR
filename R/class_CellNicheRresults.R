#' Construct a CellNicheRresults Object
#'
#' Create a standardized results container for ROI-based spatial niche analysis.
#'
#' @param roi_cells A data frame containing per-cell ROI assignments and metadata.
#' @param roi_summary A data frame containing summarized per-ROI statistics.
#' @param meta (optional) The metadata table used as input, if available.
#' @param params (optional) A list of pipeline parameters (r_mm, n_controls, etc.).
#'
#' @return An object of class \code{"CellNicheRresults"}.
#' @export
new_CellNicheRresults <- function(roi_cells, roi_summary, meta = NULL, params = list()) {
  stopifnot(is.data.frame(roi_cells), is.data.frame(roi_summary))
  structure(
    list(
      roi_cells = roi_cells,
      roi_summary = roi_summary,
      meta = meta,
      params = params
    ),
    class = "CellNicheRresults"
  )
}


#' @export
print.CellNicheRresults <- function(x, ...) {
  cat("CellNicheRresults object\n")
  cat(" ├─", nrow(x$roi_cells), "cells across",
      length(unique(x$roi_cells$roi_id)), "ROIs\n")
  cat(" ├─ ROI types:", paste(unique(x$roi_cells$roi_type), collapse = ", "), "\n")
  if (!is.null(x$params)) {
    cat(" └─ Parameters:",
        paste(names(x$params), "=", unlist(x$params)), "\n")
  }
  invisible(x)
}