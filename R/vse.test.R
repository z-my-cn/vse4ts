#' @title Variance Scale Exponent Test
#' @description
#' Use the variance scale exponent method to construct a hypothesis test about long memory in time series.
#' @param x A time series vector.
#' @param m A parameter to control the number of scales. Default is 0.5.
#' @return A list with class "vsetest" containing the following components:
#' \item{SLmemory}{the test statistic of the Short-Long memory test.}
#' \item{df}{the degrees of freedom of the test.}
#' \item{p.value}{the p-value of the test.}
#' @references
#' Fu, H., Chen, W., & He, X.-J. (2018). On a class of estimation and test for long memory. In Physica A: Statistical Mechanics and its Applications (Vol. 509, pp. 906–920). Elsevier BV. https://doi.org/10.1016/j.physa.2018.06.092
#' @examples
#' set.seed(123)
#' x <- rnorm(1024)
#' vse.test(x)
#' @export
vse.test<-function(x, m=0.5){
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
    Pvalue<-min(pchisq(Sta,Degr),1-pchisq(Sta,Degr))

    # 返回检验统计量，自由度和p值
    # return the test statistic, degree of freedom, and p-value
    # list(Sta=Sta,Degr=Degr,Pvalue=Pvalue)
    # return(list(Pvalue=Pvalue))
    result <- list(SLmemory=SLmemory, df=Degr, p.value=Pvalue)
    class(result) <- "vsetest"
    return(result)
}


