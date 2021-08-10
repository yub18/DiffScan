#' Initialize DiffScan
#'
#' Initialize with optional quality control.
#'
#' Set \code{qc_ctrl = FALSE} to skip QC. Otherwise, \code{qc_ctrl} should be a list
#' with components \code{n_trim_tail} and \code{bad_sites}. \code{n_trim_tail} nucleotide
#' positions in the tail of each transcript will be trimmed. Nucleotide positions in
#' \code{bad_sites} will be removed.
#'
#' Transcripts are cut into segments with length not exceeding 100 nt.
#'
#' @param r A list of reactivities of two conditions.
#' @param qc_ctrl Quality control settings. See "Details".
#' @param ncores Number of cores for parallel computation. Automatically determined by default.
#'
#' @return Initialized reactivities.
#' @export
init <- function(r, qc_ctrl = FALSE, ncores = "auto") {
  library(pipeR)
  if(.Platform$OS.type == 'windows') ncores <- 1
  else {
    if (ncores == "auto") ncores <- parallel::detectCores() - 1
  }
  cut_into_seg <- function(r_rna_cut, nseg) {
    if (length(r_rna_cut) == 0) {
      return(NULL)
    }
    segs <- Reduce(c, r_rna_cut %>>% rlist::list.map({
      r_seg <- dplyr::as_tibble(.)
      abs_start_seg <- strsplit(.name, "_to_")[[1]][1] %>>% as.integer()
      intvs <- split(seq_len(nrow(r_seg)), ceiling(seq_len(nrow(r_seg)) / nseg))
      names(intvs) <- intvs %>>% rlist::list.mapv({
        abs_start <- .[1] + abs_start_seg - 1
        abs_stop <- .[length(.)] + abs_start_seg - 1
        paste0(abs_start, "_to_", abs_stop)
      })
      intvs %>>% rlist::list.exclude(length(.) < 5) %>>% rlist::list.map(r_seg[., ])
    }))
    if (length(segs) == 0) {
      return(NULL)
    }
    segs
  }

  message(length(r), " rnas in total")
  nseg <- 100

  if (identical(FALSE, qc_ctrl)) {
    message("skip quality control")
    r <- parallel::mcmapply(
      FUN = function(r_rna) {
        r_rna_cut <- list()
        r_rna_cut[[paste0(1, "_to_", nrow(r_rna))]] <- r_rna
        cut_into_seg(r_rna_cut = r_rna_cut, nseg = nseg)
      },
      r_rna = r, SIMPLIFY = F, mc.cores = ncores
    )
    return(r)
  }

  if (!("n_trim_tail" %in% names(qc_ctrl))) qc_ctrl$n_trim_tail <- 0
  if (!("bad_sites" %in% names(qc_ctrl))) qc_ctrl$bad_sites <- list()
  cut_bad_sites <- function(r_rna, bad_sites_rna) {
    if (length(bad_sites_rna) == 0) {
      new <- list(r_rna)
      names(new) <- paste0(1, "_to_", nrow(r_rna))
      return(new)
    }
    whether_site_cut <- (seq_len(nrow(r_rna)) %in% bad_sites_rna)
    if (all(whether_site_cut)) {
      return(NULL)
    }
    rle_site_cut <- rle(whether_site_cut)
    r_cut <- list()
    pos <- 1
    for (k in seq_along(rle_site_cut$values)) {
      if ((!rle_site_cut$values[k]) & (rle_site_cut$lengths[k] > 1)) {
        reach <- pos + rle_site_cut$lengths[k] - 1
        r_cut[[paste0(pos, "_to_", reach)]] <- r_rna[pos:reach, ]
      }
      pos <- pos + rle_site_cut$lengths[k]
    }
    if (length(r_cut) == 0) {
      return(NULL)
    }
    r_cut
  }
  r <- parallel::mcmapply(
    FUN = function(r_rna, bad_sites_rna) {
      r_rna <- r_rna[seq_len(nrow(r_rna) - qc_ctrl$n_trim_tail), ]
      r_rna_cut <- cut_bad_sites(r_rna = r_rna, bad_sites_rna = bad_sites_rna)
      cut_into_seg(r_rna_cut = r_rna_cut, nseg = nseg)
    },
    r_rna = r, bad_sites_rna = qc_ctrl$bad_sites[names(r)], SIMPLIFY = F, mc.cores = ncores
  ) %>>% rlist::list.exclude(identical(NULL, .))

  message(length(r), " rnas kept after QC")
  return(r)
}
