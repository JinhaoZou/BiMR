\name{Data}
\alias{Data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Data
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Data(n = 1000, b01 = 1, b02 = 1, p1 = 0.3, p2 = 0.3, seed = 321, g12 = 0.3, g21 = 0.3, b11w = 0.1, b22w = 0.1, nIVw = 0, b11s = 1, b22s = 1, nIVs = 0, con1_var = 1, Ve1 = 0.1, Ve2 = 0.1, bcy1 = 1, bcy2 = 1, Causal = c("uni", "bi_infi", "bi_fl"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{b01}{
%%     ~~Describe \code{b01} here~~
}
  \item{b02}{
%%     ~~Describe \code{b02} here~~
}
  \item{p1}{
%%     ~~Describe \code{p1} here~~
}
  \item{p2}{
%%     ~~Describe \code{p2} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
}
  \item{g12}{
%%     ~~Describe \code{g12} here~~
}
  \item{g21}{
%%     ~~Describe \code{g21} here~~
}
  \item{b11w}{
%%     ~~Describe \code{b11w} here~~
}
  \item{b22w}{
%%     ~~Describe \code{b22w} here~~
}
  \item{nIVw}{
%%     ~~Describe \code{nIVw} here~~
}
  \item{b11s}{
%%     ~~Describe \code{b11s} here~~
}
  \item{b22s}{
%%     ~~Describe \code{b22s} here~~
}
  \item{nIVs}{
%%     ~~Describe \code{nIVs} here~~
}
  \item{con1_var}{
%%     ~~Describe \code{con1_var} here~~
}
  \item{Ve1}{
%%     ~~Describe \code{Ve1} here~~
}
  \item{Ve2}{
%%     ~~Describe \code{Ve2} here~~
}
  \item{bcy1}{
%%     ~~Describe \code{bcy1} here~~
}
  \item{bcy2}{
%%     ~~Describe \code{bcy2} here~~
}
  \item{Causal}{
%%     ~~Describe \code{Causal} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (n = 1000, b01 = 1, b02 = 1, p1 = 0.3, p2 = 0.3, seed = 321,
    g12 = 0.3, g21 = 0.3, b11w = 0.1, b22w = 0.1, nIVw = 0, b11s = 1,
    b22s = 1, nIVs = 0, con1_var = 1, Ve1 = 0.1, Ve2 = 0.1, bcy1 = 1,
    bcy2 = 1, Causal = c("uni", "bi_infi", "bi_fl"), ...)
{
    set.seed(seed)
    e1 <- rnorm(n, 0, Ve1)
    e2 <- rnorm(n, 0, Ve2)
    nIV <- nIVw + nIVs
    B11 = c(abs(rnorm(nIVw, b11w, 0.01)), abs(rnorm(nIVs, b11s,
        0.01)))
    B22 = c(abs(rnorm(nIVw, b22w, 0.01)), abs(rnorm(nIVs, b22s,
        0.01)))
    con1 <- abs(rnorm(n, 1, con1_var))
    X1 <- sapply(abs(rnorm(nIV, p1, 0.01)), FUN = function(x) rbinom(n,
        2, x))
    X2 <- sapply(abs(rnorm(nIV, p2, 0.01)), FUN = function(x) rbinom(n,
        2, x))
    if (Causal == "uni") {
        Y1 <- b01 + X1 \%*\% B11 + bcy1 * con1 + e1
        Y2 <- b02 + g12 * Y1 + bcy2 * con1 + X2 \%*\% B22 + e2
    }
    else if (Causal == "bi_infi") {
        Y1 <- (b01 + g21 * b02 + g21 * bcy2 * con1 + X2 \%*\% B22 *
            g21 + g21 * e2 + bcy1 * con1 + X1 \%*\% B11 + e1)/(1 -
            g21 * g12)
        Y2 <- (b02 + g12 * b01 + g12 * bcy1 * con1 + X1 \%*\% B11 *
            g12 + g12 * e1 + bcy2 * con1 + X2 \%*\% B22 + e2)/(1 -
            g21 * g12)
    }
    else {
        iter_sta <- sample(2:20, n, replace = T)
        Y_chos <- rbinom(n, 1, 0.5)
        Gen_Y <- function(i) {
            iter_sta_one <- iter_sta[i]
            Y_chos_one <- Y_chos[i]
            X1_one <- X1[i, ]
            X2_one <- X2[i, ]
            con1_one <- con1[i]
            e1_one <- e1[i]
            e2_one <- e2[i]
            Y1_one <- rnorm(1, b01, 0.01)
            Y2_one <- rnorm(1, b02, 0.01)
            if (Y_chos_one == 0) {
                for (j in 1:iter_sta_one) {
                  Y1_one_pre <- Y1_one
                  Y2_one <- b02 + B22 \%*\% X2_one + bcy2 * con1_one +
                    g12 * Y1_one_pre + e2_one
                  Y2_one_pre <- Y2_one
                  Y1_one <- b01 + B11 \%*\% X1_one + bcy1 * con1_one +
                    g21 * Y2_one_pre + e1_one
                }
            }
            else {
                for (j in 1:iter_sta_one) {
                  Y2_one_pre <- Y2_one
                  Y1_one <- b01 + B11 \%*\% X1_one + bcy1 * con1_one +
                    g21 * Y2_one_pre + e1_one
                  Y1_one_pre <- Y1_one
                  Y2_one <- b02 + B22 \%*\% X2_one + bcy2 * con1_one +
                    g12 * Y1_one_pre + e2_one
                }
            }
            return(c(Y1_one = Y1_one, Y2_one = Y2_one, Y1_one_pre = Y1_one_pre,
                Y2_one_pre = Y2_one_pre))
        }
        Y <- sapply(seq(n), FUN = function(x) Gen_Y(x))
        Y1 <- Y[1, ]
        Y2 <- Y[2, ]
        Y1_pre <- Y[3, ]
        Y2_pre <- Y[4, ]
    }
    Fst_X1 <- summary(lm(Y1 ~ X1))$fstatistic[[1]]
    Fst_X2 <- summary(lm(Y2 ~ X2))$fstatistic[[1]]
    return(list(X1 = X1, X2 = X2, Y1 = Y1, Y2 = Y2, con = con1,
        Fst_X1 = Fst_X1, Fst_X2 = Fst_X2))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
