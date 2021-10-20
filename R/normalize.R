#' Normalize reactivities
#'
#' Remove systematic bias before differential analysis. Skip this step
#' if the user is convinced of the comparability between reactivities.
#'
#' Users should inspect normalized reactivities carefully.
#'
#' @param r Initialized reactivities from \code{\link{init}}.
#' @param fake0 A logical variable indicating whether there are a substantial number
#' of fake 0's in reactivities (\code{FALSE} by default). Fake 0's arises when SP
#' experiments fail to measure the structure of part of the nucleotide positions.
#' @param plotfile File name of plot of details. Set to \code{FALSE} (default) to skip.
#' @param ncores Number of cores for parallel computation. Automatically determined by default.
#' @param verbose A logical variable indicating whether to print details (\code{TRUE} by default).
#'
#' @return A list containing the following components:
#' \item{\code{r}}{Normalized reactivities.}
#' \item{\code{fit}}{The fitted model for the normalization.}
#' @export
#'
normalize <- function(r, fake0 = FALSE, plotfile = FALSE, ncores = "auto", verbose = TRUE) {
  if(.Platform$OS.type == 'windows') ncores <- 1
  else {
    if (ncores == "auto") ncores <- parallel::detectCores() - 1
  }

  message("merge transcripts")
  new_r <- parallel::mcmapply(
    FUN = function(r_rna) {
      r_rna <- rlist::list.rbind(r_rna)
      if (any(is.na(r_rna)) | any(is.infinite(unlist(r_rna)))) {
        message("perform quality control in init first")
        stop()
      }
      r_rna
    },
    r_rna = r, SIMPLIFY = F, mc.cores = ncores
  ) %>>% rlist::list.rbind()

  trace_bindrow <- function(x, type = c("trace_rna", "trace_seg")) {
    if (type == "trace_rna") {
      nrows <- x %>>% rlist::list.mapv(sum(sapply(., nrow)))
    } else {
      nrows <- sapply(x, nrow)
    }
    marks <- c(0, cumsum(nrows))
    id_row <- lapply(2:length(marks), function(i) (marks[i - 1] + 1):(marks[i]))
    names(id_row) <- names(x)
    id_row
  }
  intv_to_traceback_rna <- trace_bindrow(x = r, type = "trace_rna")
  intv_to_traceback_seg <- parallel::mcmapply(
    FUN = function(r_rna_cut) trace_bindrow(x = r_rna_cut, type = "trace_seg"),
    r_rna_cut = r, SIMPLIFY = F, mc.cores = ncores
  )
  rm(r)

  helper <- function(rA, rB) {
    len <- nrow(rA)
    if (!identical(FALSE, plotfile)) {
      plot_helper <- function(A, B, eps, main) {
        graphics::pairs(cbind(logA = log(A + eps), logB = log(B + eps)), lower.panel = NULL, main = main)
      }
      grDevices::pdf(plotfile)
      plot_id <- sample(1:len, size = min(len, 1e3))
      plot_helper(A = rA[plot_id, ], B = rB[plot_id, ], eps = 1 - min(unlist(rA), unlist(rB)), main = "Raw")
    }

    if (verbose) message("normalize to [0,1]")
    extreme_message <- "Extreme reactivities. Perform corresponding quality control first."
    pre_norm <- function(raw) {
      low <- which(raw <= stats::quantile(raw, 0.05))
      high <- which(raw >= stats::quantile(raw, 0.95))
      mid <- setdiff(seq_along(raw), union(low, high))
      if (length(mid) < 5) stop(extreme_message)
      r <- c()
      r[low] <- 0
      r[high] <- 1
      r[mid] <- raw[mid] - min(raw[mid])
      if (max(r[mid]) == 0) stop(extreme_message)
      r[mid] <- r[mid] / max(r[mid])
      r
    }
    rA <- apply(rA, 2, pre_norm)
    rB <- apply(rB, 2, pre_norm)
    if (!identical(FALSE, plotfile)) {
      plot_helper(A = rA[plot_id, ], B = rB[plot_id, ], eps = 1, main = "Scale to 0~1")
    }

    if (verbose) message("quantile normalization")
    qtnorm <- function(X) {
      if (fake0) {
        X_norm <- X
        set <- which(apply(X, 1, function(z) all(z != 0)))
        if (length(set) < 2) stop(extreme_message)
        X_norm[set, ] <- preprocessCore::normalize.quantiles(X[set, ])
        set1 <- setdiff(seq_len(nrow(X)), set)
        if (length(set1) < 2) stop("Set fake0=FALSE")
        for (j in seq_len(ncol(X))) {
          setj <- intersect(set1, which(X[, j] != 0))
          approxs <- stats::approx(
            x = X[set, j], y = X_norm[set, j], xout = X[setj, j]
          )$y
          X_norm[setj, j] <- ifelse(is.na(approxs), X_norm[setj, j], approxs)
        }
      } else {
        X_norm <- preprocessCore::normalize.quantiles(X)
      }
      colnames(X_norm) <- colnames(X)
      X_norm
    }
    if (ncol(rA) > 1) rA <- qtnorm(rA)
    if (ncol(rB) > 1) rB <- qtnorm(rB)
    if ((ncol(rA) > 1 | ncol(rB) > 1) & !identical(FALSE, plotfile)) {
      plot_helper(A = rA[plot_id, ], B = rB[plot_id, ], eps = 1, main = "Quantile normalization")
    }

    meanA <- rowMeans(rA)
    meanB <- rowMeans(rB)
    if (verbose) {
      message("meanA:")
      print(summary(meanA))
      message("meanB:")
      print(summary(meanB))
    }

    if (verbose) message("determine invariant set")
    extract_inv_ss <- function(r_rep, upper_qt) {
      which(r_rep >= stats::quantile(r_rep[(r_rep != 1) & (r_rep != 0)],
        probs = upper_qt, na.rm = T
      ))
    }
    extract_inv_ds <- function(r_rep, lower_qt) {
      which(r_rep <= stats::quantile(r_rep[(r_rep != 1) & (r_rep != 0)],
        probs = lower_qt, na.rm = T
      ))
    }
    extreme_sites <- if (!fake0) {
      intersect(which(meanA == 0 | meanA == 1), which(meanB == 0 | meanB == 1))
    } else {
      union(which(meanA == 0 | meanB == 0), which(meanA == 1 & meanB == 1))
    }
    grids <- data.frame(
      lower = rep(seq(from = 0.05, to = 0.4, by = 0.05), each = 8),
      upper = rep(seq(from = 0.95, to = 0.6, by = -0.05), times = 8)
    )
    res <- parallel::mcmapply(
      FUN = function(lower_qt, upper_qt) {
        inv_ss <- Reduce(intersect, c(
          plyr::alply(rA, 2, extract_inv_ss, upper_qt = upper_qt),
          plyr::alply(rB, 2, extract_inv_ss, upper_qt = upper_qt)
        ))
        inv_ds <- Reduce(intersect, c(
          plyr::alply(rA, 2, extract_inv_ds, lower_qt = lower_qt),
          plyr::alply(rB, 2, extract_inv_ds, lower_qt = lower_qt)
        ))
        inv_ss <- setdiff(inv_ss, extreme_sites)
        inv_ds <- setdiff(inv_ds, extreme_sites)
        inv <- union(inv_ss, inv_ds)
        if (length(inv) >= max(10, len / 20)) {
          cor_inv <- tryCatch(
            stats::cor(meanA[inv], meanB[inv], method = "spearman"),
            warning = function(w) -2
          )
        } else {
          cor_inv <- -2
        }
        list(
          cor_inv = cor_inv, inv = inv, lower_qt = lower_qt,
          upper_qt = upper_qt, inv_ss = inv_ss, inv_ds = inv_ds
        )
      },
      lower_qt = grids$lower, upper_qt = grids$upper, SIMPLIFY = F, mc.cores = ncores
    )
    res <- res[[which.max(rlist::list.mapv(res, cor_inv))]]
    if (res$cor_inv == -2) stop("Fail to determine the invariant set")
    inv <- res$inv
    info_inv <- c(
      max_cor = res$cor_inv,
      len_inv = length(inv), prop_inv = length(inv) / len,
      lower_qt = res$lower_qt, upper_qt = res$upper_qt,
      prop_inv_ss = round(length(res$inv_ss) / length(inv), 2),
      prop_inv_ds = round(length(res$inv_ds) / length(inv), 2)
    )
    if (verbose) print(info_inv)

    if (verbose) message("fit rlm")
    epsA <- if (fake0) {
      0
    } else {
      tmpA <- meanA[meanA != 0]
      id_out <- which(log(tmpA) %in% graphics::boxplot(log(tmpA), plot = F)$out)
      ifelse(length(id_out) > 0, min(tmpA[-id_out]), min(tmpA))
    }
    epsB <- if (fake0) {
      0
    } else {
      tmpB <- meanB[meanB != 0]
      id_out <- which(log(tmpB) %in% graphics::boxplot(log(tmpB), plot = F)$out)
      ifelse(length(id_out) > 0, min(tmpB[-id_out]), min(tmpB))
    }
    eps_anchor <- if (fake0) {
      NULL
    } else {
      which(meanA <= epsA & meanB <= epsB)
    }
    inv <- setdiff(inv, eps_anchor)
    if (length(inv) < min(len / 20, 10)) stop("Fail to determine the invariant set")
    log_learn <- data.frame(logA = log(meanA[inv] + epsA), logB = log(meanB[inv] + epsB))
    fit <- MASS::rlm(formula = logB ~ logA, data = log_learn, maxit = 200)
    b <- fit$coefficients
    if (any(is.na(b)) | (b[2] <= 0)) stop("missing/infeasible b")

    if (verbose) message("transform")
    if (b[1] >= 0) {
      normA <- exp(b[1]) * (rA + epsA)^(b[2]) - epsB
      normA[rA == 1] <- 1
      normA[normA > 1] <- 1
      normA[rA == 0] <- 0
      normA[normA < 0] <- 0
      normB <- rB
    } else {
      normB <- ((rB + epsB) / exp(b[1]))^(1 / b[2]) - epsA
      normB[rB == 1] <- 1
      normB[normB > 1] <- 1
      normB[rB == 0] <- 0
      normB[normB < 0] <- 0
      normA <- rA
    }

    if (!identical(FALSE, plotfile)) {
      p0 <- ggplot2::ggplot(
        data.frame(logA = log(meanA[plot_id] + 1), logB = log(meanB[plot_id] + 1)),
        ggplot2::aes(logA, logB)
      ) +
        ggplot2::geom_point() +
        ggplot2::geom_abline(intercept = 0, slope = 1, color = "red", linetype = 2) +
        ggplot2::labs(title = "Before fitting", x = "log(meanA + 1)", y = "log(meanB + 1)")
      print(p0)
      plot_id_inv <- sample(seq_along(inv), size = min(1e3, length(inv)))
      p1 <- ggplot2::ggplot(
        data.frame(x = log(epsA + meanA[plot_id]), y = log(epsB + meanB[plot_id])),
        ggplot2::aes(x, y)
      ) +
        ggplot2::geom_point(alpha = 0.1) +
        ggplot2::geom_point(
          data = log_learn[plot_id_inv, ],
          ggplot2::aes(x = logA, y = logB),
          color = "blue", alpha = 0.1
        ) +
        ggplot2::geom_line(
          data = dplyr::tibble(
            x = sort(fit$x[plot_id_inv, 2]),
            y = stats::predict(fit, newdata = data.frame(logA = x))
          ),
          ggplot2::aes(x, y),
          color = "red"
        ) +
        ggplot2::labs(
          title = paste0(
            "Fitting, b=(", format(b[1], digits = 2), ", ",
            format(b[2], digits = 2), ")"
          ),
          x = "logA", y = "logB"
        )
      print(p1)
      plot_id <- union(plot_id, sample(inv, size = min(length(inv), 100)))
      subr <- dplyr::tibble(
        A = rowMeans(as.matrix(normA[plot_id, ])),
        B = rowMeans(as.matrix(normB[plot_id, ])),
        inv = plot_id %in% inv
      )
      p2 <- ggplot2::ggplot(subr, ggplot2::aes(A, B, color = inv)) +
        ggplot2::geom_point(alpha = 0.3) +
        ggplot2::geom_abline(intercept = 0, slope = 1, color = "red", linetype = 2) +
        ggplot2::labs(title = "Normalized")
      print(p2)
      grDevices::dev.off()
    }

    return(list(
      r = dplyr::bind_cols(dplyr::as_tibble(normA), dplyr::as_tibble(normB)),
      b = b, info_inv = info_inv, eps = c(A = epsA, B = epsB)
    ))
  }
  res <- helper(
    rA = new_r %>>% dplyr::select(dplyr::starts_with("A")),
    rB = new_r %>>% dplyr::select(dplyr::starts_with("B"))
  )

  message("traceback")
  r <- parallel::mcmapply(
    FUN = function(r_rna, intvs) intvs %>>% rlist::list.map(r_rna[., ]),
    r_rna = intv_to_traceback_rna %>>% rlist::list.map(res$r[., ]),
    intvs = intv_to_traceback_seg, SIMPLIFY = F, mc.cores = ncores
  )

  return(list(r = r, fit = res[c("b", "info_inv", "eps")]))
}
