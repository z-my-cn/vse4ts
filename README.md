vse4ts
================

<!-- badges: start -->
<!-- badges: end -->

## Introduction

The goal of vse4ts is to identify memory patterns in time series using
variance scaling exponents. It constructs variance scaling exponents,
achieving the identification of memory patterns in time series,
especially long memory.

## Installation

You can install the development version of vse4ts from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
# devtools::install_github("z-my-cn/vse4ts")
```

## Example

Here is a basic example of how to use the `vse` function in the vse4ts
package:

``` r
library(vse4ts)
set.seed(123)
x <- rnorm(1024)
x.vse <- vse(x)
print(x.vse)
#> $vse
#> [1] 0.4987233
```

## License

MIT © 2024 vse4ts authors
