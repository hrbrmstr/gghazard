<!-- README.md is generated from README.Rmd. Please edit that file -->
gghazard : Improved Base and Grid Plots for Survival Hazard Cox Regression

The following functions are implemented:

The following data sets are included:

### News

-   Version 0.0.0.9000 released

### Installation

``` r
devtools::install_github("hrbrmstr/gghazard")
```

### Usage

``` r
library(survival)
library(ggplot2)
library(gghazard)
#> Loading required package: gridExtra

# current verison
packageVersion("gghazard")
#> [1] '0.0.0.9000'

vfit <- coxph(Surv(time, status) ~ trt + factor(celltype) +
              karno + age, data=veteran, x=TRUE)
```

``` r
gg_cox_zph(fortify(cox.zph(vfit)), pred_col="red")
```

<img src="README-unnamed-chunk-2-1.png" title="" alt="" width="864" />

``` r
plot(cox.zph(vfit), smoothcol="red", confcol="red")
```

<img src="README-unnamed-chunk-3-1.png" title="" alt="" width="864" /><img src="README-unnamed-chunk-3-2.png" title="" alt="" width="864" /><img src="README-unnamed-chunk-3-3.png" title="" alt="" width="864" /><img src="README-unnamed-chunk-3-4.png" title="" alt="" width="864" /><img src="README-unnamed-chunk-3-5.png" title="" alt="" width="864" /><img src="README-unnamed-chunk-3-6.png" title="" alt="" width="864" />

### Test Results

``` r
library(gghazard)
library(testthat)

date()
#> [1] "Sun Aug 30 23:13:45 2015"

test_dir("tests/")
#> testthat results ========================================================================================================
#> OK: 0 SKIPPED: 0 FAILED: 0
#> 
#> DONE
```

### Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.
