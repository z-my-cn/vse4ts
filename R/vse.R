#' @title Calculate the Variance Scale Exponent of a Time Series
#' @description
#' Calculate the variance scale exponent of a time series.
#' @param x A time series vector.
#' @param m A parameter to control the number of scales. Default is 0.5.
#' @param n The number of scales. If `NULL`, it will be calculated as `floor(N^m)`.
#' @param type The type of variance scale exponent. Default is "weak".
#' @return The variance scale exponent.
#' @references
#' Fu, H., Chen, W., & He, X.-J. (2018). On a class of estimation and test for long memory. In Physica A: Statistical Mechanics and its Applications (Vol. 509, pp. 906–920). Elsevier BV. https://doi.org/10.1016/j.physa.2018.06.092
#' @examples
#' ## Compute the variance scale exponent of a time series
#' # Generate a random time series
#' set.seed(123)
#' x <- rnorm(1024) # F = H = 0.5 also d = 0
#' vse(x)
#'
#' ## Compare the result with the Hurst exponent
#' library(pracma)
#'
#' # A time series with Hurst exponent 0.72
#' data("brown72")
#' x <- brown72     # F = H = 0.72 also d = 0.22
#' hurstexp(x)
#' vse(x)
#'
#' # A time series with Hurst exponent 0.43
#' xlm <- numeric(1024); xlm[1] <- 0.1
#' for (i in 2:1024) xlm[i] <- 4 * xlm[i-1] * (1 - xlm[i-1])
#' x <- xlm         # F = H = 0.43 also d = -0.07
#' hurstexp(x)
#' vse(x)
#'
#' @export
vse <- function(x, m=0.5, n = NULL, type=c("weak", "strong")){
    # 确保输入参数符合要求
    stopifnot(is.numeric(x), is.numeric(m))
    type <- match.arg(type)

    # 计算需要的迭代次数
    # calculate the number of iterations needed
    N <- length(x)
    # n <- floor(N^m)
    if (is.null(n)) {
        n <- floor(N^m)
    } else {
        stopifnot(is.numeric(n), n > 0, n == as.integer(n))
    }

    # 内部函数，计算时间序列的方差
    # inner function, calculate the variance of the time series
    vartss <- function(TS, j) {
        r <- length(TS) %/% j
        TS1 <- TS[1:(r * j)]
        TS2 <- matrix(TS1, j, r)
        TS3 <- colSums(TS2)
        return(list(varTS=var(TS3), Nadj=r * j))
    }

    # 计算方差并进行加权平均
    # calculate the variance and perform a weighted average
    # varTSjj <- numeric(n)
    varTS <- function(j){
        vartss.array <- sapply(1:j,function(i) {
            as.numeric(vartss(x[i:N], j))
        })
        return(sum(vartss.array[1,]*vartss.array[2,])/(sum(vartss.array[2,])))
    }
    # 应用varTS函数到1到n的序列上
    # apply the varTS function to the sequence from 1 to n
    varTSjj <- sapply(1:n, varTS)

    # 创建时间序列数据和其方差
    # create time series data and its variance
    ts <- list(tsn=1:n, tsvar=varTSjj)

    if(type=="weak"){
        # 估计弱方差标度指数
        # 优化函数，用于拟合参数
        # optimization function, used to fit parameters
        fn <- function(theta, x, y) {
            sum((y - theta[1]*x^(2*theta[2]))^2)
        }

        # 使用非线性最小二乘法进行参数估计
        # use nonlinear least squares to estimate parameters
        out <- nlminb(c(1, 0.5), fn, x = ts$tsn, y = ts$tsvar / ts$tsvar[1])

        # 返回估计得到的参数
        # return the estimated parameters
        return(out$par[2])
    }else if(type=="strong"){
        # 估计强方差标度指数
        # 优化函数，用于拟合参数
        # optimization function, used to fit parameters
        fn <- function(theta, x, y) {
            sum((y - x^(2*theta[1]))^2)
        }

        # 使用非线性最小二乘法进行参数估计
        # use nonlinear least squares to estimate parameters
        out <- nlminb(c(0.5), fn, x = ts$tsn, y = ts$tsvar / ts$tsvar[1])

        # 返回估计得到的参数
        # return the estimated parameters
        return(out$par[1])
    }else{
        stop("type must be 'weak' or 'strong'")
    }
}

#' @title Testing White Noise in Time Series
#' @description
#' The function Wnoise.test computes the test statistic for white noise in time series based on the variance scale exponent.
#' The null hypothesis is that the time series is independent white noise, while the alternative hypothesis is that the time series is a non-independent stochastic process.
#' @param x A time series vector.
#' @param m A parameter to control the number of scales. Default is 0.5.
#' @param n The number of scales. If `NULL`, it will be calculated as `floor(N^m)`.
#' @return A list with class "Wnoise.test" containing the following components:
#' \item{Wnoise}{the test statistic }
#' \item{df}{the degrees of freedom of the test.}
#' \item{p.value}{the p-value of the test.}
#' @references
#' Fu, H., Chen, W., & He, X.-J. (2018). On a class of estimation and test for long memory. In Physica A: Statistical Mechanics and its Applications (Vol. 509, pp. 906–920). Elsevier BV. https://doi.org/10.1016/j.physa.2018.06.092
#' @examples
#' ## Test white noise in time series
#' library(pracma)
#'
#' set.seed(123)
#' data("brown72")
#' x72 <- brown72                          #  H = 0.72
#' xgn <- rnorm(1024)                      #  H = 0.50
#' xlm <- numeric(1024); xlm[1] <- 0.1     #  H = 0.43
#' for (i in 2:1024) xlm[i] <- 4 * xlm[i-1] * (1 - xlm[i-1])
#'
#' Wnoise.test(x72)
#' Wnoise.test(xgn)
#' Wnoise.test(xlm)
#'
#' @export
Wnoise.test<-function(x, m=0.5, n = NULL){
    # 确保输入参数为数值类型
    # ensure the input parameters are numeric
    stopifnot(is.numeric(x), is.numeric(m))

    # 计算需要的迭代次数
    # calculate the number of iterations needed
    N <- length(x)
    # n <- floor(N^m)
    if (is.null(n)) {
        n <- floor(N^m)
    } else {
        stopifnot(is.numeric(n), n > 0, n == as.integer(n))
    }

    # 内部函数，计算时间序列的方差
    # inner function, calculate the variance of the time series
    vartss <- function(TS, j) {
        r <- length(TS) %/% j
        TS1 <- TS[1:(r * j)]
        TS2 <- matrix(TS1, j, r)
        TS3 <- colSums(TS2)
        return(list(varTS=var(TS3), Nadj=r * j))
    }

    # 计算方差并进行加权平均
    # calculate the variance and perform a weighted average
    varTS <- function(j){
        vartss.array <- sapply(1:j,function(i) {
            as.numeric(vartss(x[i:N], j))
        })
        return(sum(vartss.array[1,]*vartss.array[2,])/(sum(vartss.array[2,])))
    }

    # 计算检验统计量，自由度和p值
    # calculate the test statistic, degree of freedom, and p-value
    Wnoise <-(floor(N/n)-1)*varTS(n)/varTS(1)/n
    Degr<-floor(N/n)-1
    Pvalue<-min(pchisq(Wnoise, Degr),1-pchisq(Wnoise, Degr))

    # 返回检验统计量，自由度和p值
    # return the test statistic, degree of freedom, and p-value
    # class(result) <- "Wnoise.test"
    # # return(result)
    # print(result)
    # return(invisible(result))

    result <- structure(list(Wnoise=Wnoise,
                             df=Degr,
                             p.value=Pvalue),
                        class="Wnoise.test")
    return(result)
}


#' @exportS3Method base::print
print.Wnoise.test <- function(x, ...) {
    cat("Wnoise Test\n\n")
    cat("Wnoise statistic:", x$Wnoise, "\n")
    cat("degrees of freedom:", x$df, "\n")
    cat("p-value:", x$p.value, "\n")
    cat("\n")
    cat("alternative hypothesis: non-independent stochastic process\n")
}


#' @title Testing Long Memory in Time Series
#' @description
#' The function SLmemory.test computes the test statistic for long memory in time series based on the variance scale exponent.
#' The null hypothesis is that the time series is white noise or short memory, while the alternative hypothesis is that the time series has long memory.
#' @param x A time series vector.
#' @param m A parameter to control the number of scales. Default is 0.5.
#' @param n The number of scales. If `NULL`, it will be calculated as `floor(N^m)`.
#' @return A list with class "SLmemory.test" containing the following components:
#' \item{SLmemory}{the test statistic }
#' \item{df}{the degrees of freedom of the test.}
#' \item{p.value}{the p-value of the test.}
#' @references
#' Fu, H., Chen, W., & He, X.-J. (2018). On a class of estimation and test for long memory. In Physica A: Statistical Mechanics and its Applications (Vol. 509, pp. 906–920). Elsevier BV. https://doi.org/10.1016/j.physa.2018.06.092
#' @examples
#' ## Test long memory in time series
#' library(pracma)
#'
#' set.seed(123)
#' data("brown72")
#' x72 <- brown72                          #  H = 0.72
#' xgn <- rnorm(1024)                      #  H = 0.50
#' xlm <- numeric(1024); xlm[1] <- 0.1     #  H = 0.43
#' for (i in 2:1024) xlm[i] <- 4 * xlm[i-1] * (1 - xlm[i-1])
#'
#' SLmemory.test(x72)
#' SLmemory.test(xgn)
#' SLmemory.test(xlm)
#'
#' @export
SLmemory.test<-function(x, m=0.5, n = NULL){
    # 确保输入参数为数值类型
    # ensure the input parameters are numeric
    stopifnot(is.numeric(x), is.numeric(m))

    # 计算需要的迭代次数
    # calculate the number of iterations needed
    N <- length(x)
    # n <- floor(N^m)
    if (is.null(n)) {
        n <- floor(N^m)
    } else {
        stopifnot(is.numeric(n), n > 0, n == as.integer(n))
    }

    # 内部函数，计算时间序列的方差
    # inner function, calculate the variance of the time series
    vartss <- function(TS, j) {
        r <- length(TS) %/% j
        TS1 <- TS[1:(r * j)]
        TS2 <- matrix(TS1, j, r)
        TS3 <- colSums(TS2)
        return(list(varTS=var(TS3), Nadj=r * j))
    }

    # 计算方差并进行加权平均
    # calculate the variance and perform a weighted average
    varTS <- function(j){
        vartss.array <- sapply(1:j,function(i) {
            as.numeric(vartss(x[i:N], j))
        })
        return(sum(vartss.array[1,]*vartss.array[2,])/(sum(vartss.array[2,])))
    }

    # 计算检验统计量，自由度和p值
    # calculate the test statistic, degree of freedom, and p-value
    SLmemory <-(floor(N/n)-1)*varTS(n)/(varTS(n+1)-varTS(n))/n
    Degr<-floor(N/n)-1
    Pvalue<-min(pchisq(SLmemory, Degr),1-pchisq(SLmemory, Degr))

    # 返回检验统计量，自由度和p值
    # return the test statistic, degree of freedom, and p-value
    # result <- list(SLmemory=SLmemory, df=Degr, p.value=Pvalue)
    # class(result) <- "SLmemory.test"
    result <- structure(list(SLmemory=SLmemory,
                             df=Degr,
                             p.value=Pvalue),
                        class="SLmemory.test")
    return(result)
}

#' @exportS3Method base::print
print.SLmemory.test <- function(x, ...) {
    cat("SLmemory Test\n\n")
    cat("SLmemory statistic:", x$SLmemory, "\n")
    cat("degrees of freedom:", x$df, "\n")
    cat("p-value:", x$p.value, "\n")
    cat("\n")
    cat("alternative hypothesis: long memory\n")
}


