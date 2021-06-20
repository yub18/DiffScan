
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DESSERT

<!-- badges: start -->

<!-- badges: end -->

DESSERT is a statistical framework for DiffErential analySis of RNA
Structure probing Experiments in high ResoluTion. Specifically, it
identifies RNA structurally variable regions between two cellular
conditions.

DESSERT performs differential analysis with nucleotide resolution,
rigorous error control, and high statistical power. Moreover, it is
compatible with various RNA structure probing platforms.

## Installation

Install DESSERT in [R](https://cran.r-project.org/) or
[Rstudio](https://www.rstudio.com/products/rstudio/) with:

``` r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
    
devtools::install_github("yub18/DESSERT")
```

## Usage

Suppose there are two RNA transcripts with SP data in condition A and B,
with 2 replicates in each condition. To identify SVRs in the
transcripts, first format the data like the following.

``` r
library(DESSERT)
library(dplyr)
set.seed(123)
r = list(
  t1 = tibble(A1 = runif(50), A2 = A1 + runif(50, -0.1, 0.1),
              B1 = A1 + runif(50, -0.1, 0.1), B2 = B1 + runif(50, -0.1, 0.1)),
  # t1: 50 nt, without SVRs
  t2 = tibble(A1 = runif(150), A2 = A1 + runif(150, -0.1, 0.1),
              B1 = c(A1[1:90] + runif(90, -0.1, 0.1), rep(0, 10), 
                     A1[101:130] + runif(30, -0.1, 0.1), rep(1, 20)), 
              B2 = B1 + runif(150, -0.1, 0.1))
  # t2: 150 nt, 2 SVRs covering nucleotide positions 91~100 and 131~150
)
```

Second, initialize diffScan as follows. (See the document of  for
appropriate quality control.)

``` r
r_init = init(r, qc_ctrl = F, ncores = 1)
```

Third, perform normalization to remove systematic bias in reactivities.
(Users convinced of the comparability between reactivities can skip this
step.)

``` r
r_norm = normalize(r_init, ncores = 1)
```

Next, run the Scan algorithm to identify SVRs.

``` r
res = scan(r_norm$r, seed = 123, alpha = 0.05, ncores = 1)
```

## Citation
