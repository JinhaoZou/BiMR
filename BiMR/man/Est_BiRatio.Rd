\name{Est_BiRatio}
\alias{Est_BiRatio}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Est_BiRatio
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Est_BiRatio(data = data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
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
function (data = data)
{
    nX1 <- dim(as.matrix(data$X1))[2]
    nX2 <- dim(as.matrix(data$X2))[2]
    nC <- dim(as.matrix(data$con))[2]
    data_variable <- cbind(data$X1, data$X2, data$con)
    colnames(data_variable) <- c(sapply(seq(nX1), FUN = function(x) paste("X1.",
        x, sep = "")), sapply(seq(nX2), FUN = function(x) paste("X2.",
        x, sep = "")), sapply(seq(nC), FUN = function(x) paste("Con.",
        x, sep = "")))
    eqY1Y2 <- Y2 ~ . - Y1
    eqY2Y1 <- Y1 ~ . - Y2
    system <- list(Y1Y2 = eqY1Y2, Y2Y1 = eqY2Y1)
    data_all <- as.data.frame(cbind(data$Y1, data$Y2, data_variable))
    colnames(data_all)[1:2] <- c("Y1", "Y2")
    fit1 <- summary(lm(eqY1Y2, data = data_all))
    fit2 <- summary(lm(eqY2Y1, data = data_all))
    coef1 <- as.data.frame(fit1$coefficients)
    est_Y2X1_g12b11 <- coef1$Estimate[2:(1 + nX1)]
    sd_Y2X1_g12b11 <- coef1$`Std. Error`[2:(1 + nX1)]
    est_Y2X2_b22 <- coef1$Estimate[(2 + nX1):(1 + nX1 + nX2)]
    sd_Y2X2_b22 <- coef1$`Std. Error`[(2 + nX1):(1 + nX1 + nX2)]
    coef2 <- as.data.frame(fit2$coefficients)
    est_Y1X1_b11 <- coef2$Estimate[2:(1 + nX1)]
    sd_Y1X1_b11 <- coef2$`Std. Error`[2:(1 + nX1)]
    est_Y1X2_g21b22 <- coef2$Estimate[(2 + nX1):(1 + nX1 + nX2)]
    sd_Y1X2_g21b22 <- coef2$`Std. Error`[(2 + nX1):(1 + nX1 +
        nX2)]
    g12_ivw <- mr_ivw(mr_input(bx = est_Y1X1_b11, bxse = sd_Y1X1_b11,
        by = est_Y2X1_g12b11, byse = sd_Y2X1_g12b11))
    g21_ivw <- mr_ivw(mr_input(bx = est_Y2X2_b22, bxse = sd_Y2X2_b22,
        by = est_Y1X2_g21b22, byse = sd_Y1X2_g21b22))
    est_g12 <- g12_ivw$Estimate
    est_g21 <- g21_ivw$Estimate
    est_g12_result <- c(g12 = est_g12, Sdg12 = g12_ivw$StdError,
        CIlowg12 = g12_ivw$CILower, CIUpg12 = g12_ivw$CIUpper)
    est_g21_result <- c(g21 = est_g21, Sdg21 = g21_ivw$StdError,
        CIlowg21 = g21_ivw$CILower, CIUpg21 = g21_ivw$CIUpper)
    return(c(est_g12_result, est_g21_result))
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
