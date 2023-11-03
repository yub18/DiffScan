#' Scan to identify structurally variable regions (SVRs)
#'
#' This function detects structurally variable regions in nucleotide resolution and
#' controls family-wise error rate.
#'
#' @param r Initialized reactivity from \code{\link{init}} or normalized reactivity
#' from \code{\link{normalize}}.
#' @param seed Random seed for reproducible results.
#' @param alpha Significance level.
#' @param N Number of replications in Monte Carlo sampling.
#' @param ncores Number of cores to use for parallel computation. Automatically determined
#' by default.
#'
#' @return A list with each element corresponding to a transcript. \code{svr} is
#' SVRs with significance score \code{P}. \code{scan} includes the start and stop positions
#' of all scanned segments.
#' @export
#'
scan <- function(r, seed, alpha = 0.05, N = 1e3, ncores = "auto") {
  dir_pkg = find.package('DiffScan')
  if(!('Qmax_null.rds' %in% list.files(paste0(dir_pkg, '/R/')))) {
    message('First time to run scan. Downloading precalculated parameters ...')
    download.file(
      url = "https://public.dm.files.1drv.com/y4mWkxGL7J2ZVa7SOUZDCBulPrRpKSMLim7Z_eT_dBL8AZIQqPvZBRt9OR1APkA8heQCCBAoLNW8FfyFGkmTtxEyZXJgBQlJ0wobPqdtoSuSITxIXaxc5XY2LeBCnwiVlorIYjwi4ZyBToPn5NYrzpHISGkKNKC37TCpfintXrj-bHah23Nnk2jXXvdYfVB1qTFekOKsZZaJ-_3zutAdYVXVoOg89mBYPOiQJElW-udi70?AVOverride=1",
      destfile = paste0(dir_pkg, '/R/Qmax_null.rds'),
      mode = 'wb'
    )
  }
  if(!('Qmax_null' %in% objects(envir = .GlobalEnv))) {
    Qmax_null <<- readRDS(paste0(dir_pkg, '/R/Qmax_null.rds'))
  }
  if(.Platform$OS.type == 'windows') ncores <- 1
  else {
    if (ncores == "auto") ncores <- parallel::detectCores() - 1
  }
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  if(.Platform$OS.type != 'windows') parallel::mc.reset.stream()
  nseg <- 100
  radius <- 2
  gamma <- 0.5
  scan_ctrl <- list(
    radius = radius, N = N,
    Q = function(pv) sum(-2 * log(pv)) / (length(pv))^gamma,
    region_range = dplyr::tibble(length = 1:100, Lmin = 1, Lmax = pmin(length, 20))
  ) # fixed, nseg and region_range are concordant with pre-calculated Qmax_null

  message("Scan ====================")
  main <- function(r, radius, N, Q, region_range, Qmax_null, alpha, ncores) {
    message("sampling null ...")
    nseg_num <- r %>>% rlist::list.mapv(sapply(., nrow)) %>>% table() %>>% {
      dplyr::tibble(nseg = as.integer(names(.)), num = .)
    }
    QGlobalMax_null <- parallel::mcmapply(
      FUN = function(k) {
        max(apply(nseg_num, 1, function(z) {
          max(sample(x = Qmax_null[[paste0("nseg=", z[1])]], size = z[2]))
        }))
      },
      k = 1:N, mc.cores = ncores
    )
    Qthr <- ifelse(alpha >= 1, -Inf, quantile(QGlobalMax_null, 1 - alpha))

    message("calculate Q and P ...")
    scan_seg_r <- function(r) {
      scan_seg <- function(obs, Q, Lmin, Lmax) {
        n <- length(obs)
        regions <- lapply(1:(n - Lmax + 1), function(k) {
          cbind(Start = k, Stop = (k + Lmin - 1):(k + Lmax - 1))
        }) %>>% rlist::list.rbind()
        if (Lmax > Lmin) {
          regions <- rbind(
            regions,
            lapply((n - Lmax + 2):(n - Lmin + 1), function(k) {
              cbind(Start = k, Stop = (k + Lmin - 1):n)
            }) %>>% rlist::list.rbind()
          )
        }
        Q <- apply(regions, 1, function(z) Q(obs[z[1]:z[2]]))
        keep_id <- which(Q >= Qthr)
        if (length(keep_id) == 0) {
          return(NULL)
        }
        P <- sapply(Q[keep_id], function(z) mean(QGlobalMax_null >= z))
        cbind(
          if (length(keep_id) > 1) {
            regions[keep_id, ]
          } else {
            t(regions[keep_id, ])
          },
          Q = Q[keep_id],
          P = P
        )
      }

      rA <- r %>>% dplyr::select(dplyr::starts_with("A"))
      rB <- r %>>% dplyr::select(dplyr::starts_with("B"))
      if (radius == 0 & ncol(rA) == 1 & ncol(rB) == 1) {
        message("set radius > 0 when single replicate is supplied")
        stop()
      }
      len <- nrow(r)
      P <- sapply(1:len, function(j) {
        Cj <- max(1, j - radius):min(len, j + radius)
        Pj <- tryCatch(
          stats::wilcox.test(
            x = unlist(rA[Cj, ]), y = unlist(rB[Cj, ]), paired = F
          )$p.value,
          warning = function(w) {
            return("Fail")
          }
        )
        if (!identical("Fail", Pj)) {
          return(Pj)
        }
        coin_test <- tryCatch(
          coin::wilcox_test(
            formula = x ~ group,
            data = dplyr::tibble(
              x = c(unlist(rA[Cj, ]), unlist(rB[Cj, ])),
              group = factor(c(
                rep("a", length(Cj) * ncol(rA)),
                rep("b", length(Cj) * ncol(rB))
              ))
            ),
            distribution = "exact"
          ),
          warning = function(w) {
            return("Error")
          }
        )
        if (identical("Error", coin_test)) {
          return(1)
        } # coin test fails for all same values
        coin::pvalue(coin_test)
      })
      scan_seg(
        obs = P,
        Q = Q,
        Lmin = region_range$Lmin[region_range$length == len],
        Lmax = region_range$Lmax[region_range$length == len]
      )
    }
    pbmcapply::pbmcmapply(
      FUN = function(r_rna_cut) {
        r_rna_cut %>>% rlist::list.map(scan_seg_r(.))
      },
      r_rna_cut = r,
      SIMPLIFY = F,
      mc.cores = ncores
    )
  }
  res <- main(
    r = r,
    radius = scan_ctrl$radius,
    N = scan_ctrl$N,
    Q = scan_ctrl$Q,
    region_range = scan_ctrl$region_range,
    Qmax_null = Qmax_null,
    alpha = alpha,
    ncores = ncores
  )

  message('\t')
  message("Scan completed, pick regions ====================")
  pick <- function(res, ncores) {
    pick_seg <- function(cands) {
      if (identical(NULL, cands)) {
        return(
          dplyr::tibble(
            Start = integer(0), Stop = integer(0),
            Q = integer(0), P = integer(0)
          )
        )
      }
      cands <- cands %>>% dplyr::as_tibble() %>>% dplyr::arrange(dplyr::desc(Q))
      picks <- NULL
      while (nrow(cands) > 0) {
        new <- cands[1, ]
        picks <- rbind(picks, new)
        cands <- cands %>>% dplyr::filter(Stop < new$Start | Start > new$Stop)
      }
      picks
    }
    cuts_to_scanned_sites <- function(cuts) {
      raw_scanned_sites <- cuts %>>% rlist::list.map(
        as.integer(strsplit(., "_to_")[[1]])
      ) %>>% rlist::list.rbind()
      if (length(cuts) == 1) {
        return(raw_scanned_sites)
      }

      scanned_sites <- matrix(raw_scanned_sites[1, ], nrow = 1)
      for (k in 2:length(cuts)) {
        new <- raw_scanned_sites[k, ]
        if (new[1] == scanned_sites[1, 2] + 1) {
          scanned_sites[1, 2] <- new[2]
        } else {
          scanned_sites <- rbind(new, scanned_sites)
        }
      }
      rownames(scanned_sites) <- NULL
      colnames(scanned_sites) <- c("Start", "Stop")
      return(scanned_sites)
    }
    pick_rna <- function(res_rna) {
      scan <- cuts_to_scanned_sites(names(res_rna))
      svr <- res_rna %>>% rlist::list.map({
        abs_start <- strsplit(.name, "_to_")[[1]][1] %>>% as.integer()
        pick_seg(.) %>>% dplyr::mutate(
          Start = Start + abs_start - 1,
          Stop = Stop + abs_start - 1
        )
      }) %>>% rlist::list.rbind() %>>% dplyr::arrange(dplyr::desc(Q))
      list(scan = scan, SVR = svr)
    }

    parallel::mcmapply(
      FUN = function(res_rna) pick_rna(res_rna = res_rna),
      res_rna = res,
      SIMPLIFY = F,
      mc.cores = ncores
    )
  }
  pick(res = res, ncores = ncores)
}
