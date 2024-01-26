vse4ts
================

<!-- badges: start -->
<!-- badges: end -->

# vse4ts

The goal of vse4ts is to identify memory patterns in time series using
variance scaling exponents. It constructs strong and weak variance
scaling exponents, achieving the identification of memory patterns in
time series, especially long memory.

## Installation

You can install the development version of vse4ts from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
# devtools::install_github("your-github-username/vse4ts")
```

## Example

Here is a basic example of how to use the `wvse` function in the vse4ts
package:

``` r
library(vse4ts)
set.seed(123)
x <- rnorm(1000)
vse <- wvse(x)
print(vse)
```
