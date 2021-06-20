#' Simulate reactivity
#'
#' Given RNA sequences, this function simulates reactivity accounting for
#' multiple structural conformations and dependencies between adjacent nucleotide positions.
#'
#' Users need to have ViennaRNA software installed to simulate structural
#' conformations from sequences.
#'
#' @param seq List of RNA sequences.
#' @param dir_ViennaRNA Directory of ViennaRNA software. See Details.
#' @param types Vector of reactivity types ('icSHAPE', 'Cordero' and 'Sokusd' supported).
#' @param ncores Number of cores for parallel computation. Automatically determined by default.
#' @param plot_cor Logical value indicating whether to plot within/between group correlations.
#' @param clear Logical value indicating whether to clear intermediate files.
#'
#' @return A list containing the following components:
#' \item{\code{r}}{Simulated reactivity.}
#' \item{\code{stru}}{A list of simulated structural conformations.}
#' \item{\code{truth}}{A list of vectors indicating positions of structurally variable regions.}
#' \item{\code{w}}{Weights of dominant conformations.}
#' @export
#'
simulate <- function(seq, dir_ViennaRNA, types, ncores = "auto", plot_cor = FALSE, clear = TRUE) {
  if(.Platform$OS.type == 'windows') ncores <- 1
  else {
    if (ncores == "auto") ncores <- parallel::detectCores() - 1
  }

  # sample conformations
  system("mkdir seqs; mkdir strus")
  if (any(sapply(seq, length) > 5000)) message("RNA strands > 5000 nt detected. This may take a long time")
  seq %>>% rlist::list.map(seqinr::write.fasta(., names = .name, file.out = paste0("seqs/", .name, ".fasta")))
  system(paste0(
    "export PATH=", dir_ViennaRNA, ":$PATH \n ",
    "parallel --bar -j ", ncores,
    " 'RNAsubopt --infile {} --stochBT_en 20 --nonRedundant",
    " --outfile $(basename {})' ::: seqs/*.fasta \n ",
    "mv *.sub strus/"
  ))
  extract_sub_from_ViennaRNA <- function(sub_file) {
    tryCatch(
      sapply(read.table(sub_file, skip = 2)[, 1], function(x) strsplit(x, "")[[1]] == ".", USE.NAMES = F),
      error = function(e) message("Wrong format:", sub_file)
    )
  }
  stru <- sapply(paste0("strus/", list.files("strus/")), extract_sub_from_ViennaRNA)
  names(stru) <- strsplit(list.files("strus/"), "\\.sub") %>>% rlist::list.mapv(.[1])
  stru <- stru %>>% rlist::list.map(unique(., MARGIN = 2)[, 1:10])

  # simulate reactivity for each conformation
  simulate_r <- function(type = c("icSHAPE", "Cordero", "Sokusd"), from = c("ds", "ss"), size) {
    type <- match.arg(type)
    from <- match.arg(from)
    if (from == "ds") {
      sample_ds <- switch(type,
        icSHAPE = function(n) {
          is0 <- sample(x = c(T, F), size = n, replace = T, prob = c(0.613, 1 - 0.613))
          sapply(is0, function(x) ifelse(x, 0, exp(rnorm(n = 1, mean = -3.40, sd = 1.88))))
        },
        Cordero = function(n) SpatialExtremes::rgev(n, loc = 0.0947, scale = 0.0672, shape = 0.2352),
        Sokusd = function(n) SpatialExtremes::rgev(n, loc = 0.0523, scale = 0.0680, shape = 0.8681)
      )
      sample_ds(n = size)
    } else {
      sample_ss <- switch(type,
        icshape = function(n) {
          is0 <- sample(x = c(T, F), size = n, replace = T, prob = c(0, 1))
          sapply(is0, function(x) ifelse(x, 0, exp(rnorm(n = 1, mean = -2.69, sd = 0.932))))
        },
        Cordero = function(n) SpatialExtremes::rgev(n, loc = 0.2198, scale = 0.1852, shape = 0.5426),
        Sokusd = function(n) rexp(n, rate = 1.4638)
      )
      sample_ss(n = size)
    }
  }
  r_stru <- sapply(types, simplify = F, function(type) {
    parallel::mcmapply(function(confs) {
      apply(confs, 2, function(s) {
        r <- rep(NA, length(s))
        r[s] <- simulate_r(type = type, from = "ss", size = sum(s))
        r[!s] <- simulate_r(type = type, from = "ds", size = sum(!s))
        r
      })
    }, confs = stru, SIMPLIFY = F, mc.cores = ncores)
  })

  # assign weights
  ## first 2 structures are dominant
  ## minor structures are assigned weights from U(0,1) and then scaled to sum to 0.1
  w_domi <- 0.9
  w_raw <- list(
    Low = list(A = c(0.2, 0.7), B = c(0.7, 0.2)),
    Medium = list(A = c(0.1, 0.8), B = c(0.8, 0.1)),
    High = list(A = c(0, 0.9), B = c(0.9, 0))
  )
  sample_minor_w <- function(n) {
    runif(n) %>>% {
      . / sum(.) * (1 - w_domi)
    }
  }
  add_domi_noise <- function(w, sd = 0.1) {
    while (1) {
      w1 <- rnorm(n = 1, mean = w[1], sd = sd)
      if (w1 >= 0 & w1 <= w_domi) break
    }
    c(w1, w_domi - w1)
  }
  w <- w_raw %>>% list.map(
    sapply(names(stru), simplify = F, function(x) {
      dplyr::tibble(
        A1 = add_domi_noise(A), A2 = add_domi_noise(A),
        B1 = add_domi_noise(B), B2 = add_domi_noise(B)
      )
    })
  )

  # simulate reactivity
  r <- sapply(r_stru, simplify = F, function(rs) {
    sapply(w, simplify = F, function(ws) {
      mapply(function(.r, .w) {
        w_full <- rbind(as.matrix(.w), replicate(4, sample_minor_w(8)))
        dplyr::as_tibble(.r %*% w_full)
      }, .r = rs, .w = ws, SIMPLIFY = F)
    })
  })

  # summary
  truth <- stru %>>% rlist::list.map(.[, 1] != .[, 2])
  if (plot_cor) {
    pdf("Correlation_between_replicates.pdf")
    par(mfrow = c(length(types), length(w_raw)))
    sapply(types, function(type) {
      r[[type]] %>>% rlist::list.map({
        cors <- sapply(., simplify = F, cor, method = "spearman")
        cors <- data.frame(
          within = sapply(cors, function(x) mean(c(x[1, 2], x[3, 4]))),
          between = sapply(cors, function(x) mean(x[1:2, 3:4]))
        )
        boxplot(cors, ylim = c(0, 1), ylab = "cor", main = paste(type, ",", .name))
        NULL
      })
    })
    dev.off()
  }
  list(r = r, stru = stru, truth = truth, w = w)
}
