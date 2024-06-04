#' @title Variance Scale Exponent
#' @description
#' Calculate the variance scaling exponent of a time series.
#' @param x A time series vector.
#' @param m A parameter to control the number of scales. Default is 0.5.
#' @return A list containing the variance scaling exponent.
#' @examples
#' set.seed(123)
#' x <- rnorm(1024)
#' vse(x)
#' @importFrom stats nlminb var
#'
#' @export
vse <- function(x, m=0.5) {
    # 确保输入参数为数值类型
    # ensure the input parameters are numeric
    stopifnot(is.numeric(x), is.numeric(m))

    # 计算需要的迭代次数
    # calculate the number of iterations needed
    N <- length(x)
    n <- floor(N^m)

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
    return(list(vse=out$par[2]))
}



