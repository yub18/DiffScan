#' Return the number of nucleotide positions in SVRs of a transcript
#'
#' @export
nSVR <- function(SVR) {
  sum(SVR$Stop - SVR$Start + 1)
}

#' Return nucleotide positions in SVRs of a transcript
#'
#' @export
vSVR <- function(SVR) {
  plyr::alply(SVR, 1, function(x) x$Start:x$Stop) %>>% unlist()
}


#' Return the number of nucleotide positions of a transcript that are scanned
#'
#' @export
nscan <- function(scan) {
  sum(scan[, 2] - scan[, 1] + 1)
}

#' Return nucleotide positions of a transcript that are scanned
#'
#' @export
vscan <- function(scan) {
  plyr::alply(scan, 1, function(x) x[1]:x[2]) %>>% unlist()
}

#' Convert rcut of a transcript into a tibble
#'
#' @param r_cut Reactivities that are cut into a list.
#'
#' @return A tibble with columns \code{pos} indicating position.
#' @export
rcut2rf <- function(r_cut) {
  r_cut %>>% rlist::list.map({
    start <- strsplit(.name, "_to_")[[1]][1] %>>% as.integer()
    dplyr::tibble(pos = start + seq_len(nrow(.)) - 1, .)
  }) %>>% rlist::list.rbind()
}

#' Convert SVRs of a transcript into a status vector
#'
#' The status vector indicates whether each nucleotide position is in SVRs.
#'
#' @param n Length of the transcript.
#'
#' @export
SVR2diff <- function(SVR, n) {
  if (nrow(SVR) == 0) {
    return(rep(F, n))
  }
  (1:n) %in% unlist(plyr::alply(SVR, 1, function(x) x$Start:x$Stop))
}


#' Reverse \code{SVR2diff}
#'
#' Reverse \code{\link{SVR2diff}}
#'
#' @param len_merge Merge SVRs that are no more \code{len_merge}'nt far.
#'
#' @export
diff2SVR <- function(diff, len_merge = 0) {
  if ((len_merge %% 1 != 0) | (len_merge < 0)) stop("len_merge must be positive integer.")
  x = which(diff)
  return(data.frame(
    Start = c(x[1], x[which(diff(x) > 1 + len_merge) + 1]),
    Stop = c(x[which(diff(x) > 1 + len_merge)], tail(x, 1))
  ))
}
