#' Calculate weak variance scaling exponents
#'
#' @param R A time series vector.
#'
#' @return The weak variance scaling exponent.
#' @export
#'
#' @examples
#' x <- rnorm(1024)
#' wvse(x)
wvse <- function(R){
    TS <- R
    vartss <- function(TS,i){
        NN=length(TS)
        r=floor(NN/i)
        cc=r*i
        TS1=TS[1:cc]
        TS2=matrix(TS1,i,r)
        TS3=apply(TS2,2,sum)
        Nadj=r*i
        varTS=stats::var(TS3)
        list(a=varTS,b=Nadj)
    }

    N = length(R)
    n = floor(N^0.5)
    varTSjj <- numeric(n)
    varTS <- function(j){
        vartss.array <- sapply(1:j,function(i) as.numeric(vartss(R[(1+i-1):N],j)))
        varTSjj[j] <- sum(vartss.array[1,]*vartss.array[2,])/(sum(vartss.array[2,]))
    }
    varTSjj <- sapply(1:n,function(j) varTS(j))
    ts <- list(tsn=c(1:n),tsvar=varTSjj)
    rhs <- function(tsn, b0, b1) {
        b0*tsn^(2*b1)
    }
    fn <- function(theta, x, y) {
        sum((y - theta[1]*x^(2*theta[2]))^2)
    }
    out <- stats::nlminb(c(1,0.5), fn, x = ts$tsn[1:n], y = ts$tsvar[1:n]/ts$tsvar[1])
    out$par[2]
}
