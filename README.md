
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DiffScan

<!-- badges: start -->

<!-- badges: end -->

DiffScan is a computational framework for differential analysis of RNA
structure probing experiments at nucleotide resolution. Specifically, it
identifies RNA structurally variable regions (SVRs) between two cellular
conditions.

DiffScan is compatible with various RNA structure probing platforms.

## System Requirements

`DiffScan` is built and tested upon [R](https://cran.r-project.org/)
(version 4.0.2). It mainly depends on the following R packages.

`dplyr, rlist, pipeR, pbmcapply, MASS, coin, plyr, ggplot2, seqinr,
SpatialExtremes, preprocessCore`

## Installation

Install `DiffScan` in [R](https://cran.r-project.org/) or
[Rstudio](https://www.rstudio.com/products/rstudio/) with:

``` r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
if (!requireNamespace("preprocessCore", quietly = TRUE))
    install.packages("http://bioconductor.org/packages/3.11/bioc/src/contrib/preprocessCore_1.50.0.tar.gz")
    
devtools::install_github("yub18/DiffScan")
```

Typical install time: a few minutes or less.

## Usage

Suppose there are two RNA transcripts (`t1, t2`) with structure probing
data in condition A and B, with 2 replicates in each condition (`A1, A2,
B1, B2`). To identify SVRs in `t1, t2`, first format the data like the
following.

``` r
library(DiffScan)
library(dplyr)
set.seed(123)
r = list(
  t1 = tibble(A1 = runif(50), A2 = A1 + runif(50, -0.1, 0.1),
              B1 = A1 + runif(50, -0.1, 0.1), B2 = B1 + runif(50, -0.1, 0.1)),
  t2 = tibble(A1 = runif(150), A2 = A1 + runif(150, -0.1, 0.1),
              B1 = c(A1[1:90] + runif(90, -0.1, 0.1), rep(0, 10), 
                     A1[101:130] + runif(30, -0.1, 0.1), rep(1, 20)), 
              B2 = B1 + runif(150, -0.1, 0.1))
)
```

In this example, transcript `t1` has no SVRs, and transcript `t2` has 2
SVRs covering nucleotide positions 91 nt \~ 100 nt and 131 nt \~ 150 nt.

Second, initialize `DiffScan` as follows. (See the document of
`DiffScan::init` for appropriate quality control.)

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

`DiffScan` predicts no SVRs for transcript `t1`, and predicts 2 SVRs for
transcript `t2` covering nucleotide positions 132 nt \~ 150 nt and 92 nt
\~100 nt.

``` r
res
#> $t1
#> $t1$scan
#>         [,1] [,2]
#> 1_to_50    1   50
#> 
#> $t1$SVR
#> # A tibble: 0 x 4
#> # ... with 4 variables: Start <dbl>, Stop <dbl>, Q <int>, P <int>
#> 
#> 
#> $t2
#> $t2$scan
#>      Start Stop
#> [1,]     1  150
#> 
#> $t2$SVR
#> # A tibble: 2 x 4
#>   Start  Stop     Q     P
#> * <dbl> <dbl> <dbl> <dbl>
#> 1   132   150  79.0     0
#> 2    92   100  54.9     0
```

Expected run time: a few minutes or less.

## Flu dataset

The following code reproduces DiffScan’s results for the Flu dataset.

``` r
data(Flu)
r = list(Flu = Flu[,1:8])
r_init = init(r, ncores = 1)
r_norm = normalize(r_init, ncores = 1)
res = scan(r_norm$r, seed = 123, alpha = 0.05, ncores = 1)
```

## Citation

[Yu, et al. Differential Analysis of RNA Structure Probing Experiments
at Nucleotide Resolution: Uncovering Regulatory Functions of RNA
Structure.
bioRxiv, 2021.](https://www.biorxiv.org/content/10.1101/2021.08.24.457484v1)

## Scripts and data

The scripts and data for the reproduction of the analyses and results in
DiffScan paper are available
[here](https://cloud.tsinghua.edu.cn/published/diffscan_scripts_and_data).
