\name{Est_Ratio}
\alias{Est_Ratio}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Est_Ratio
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Est_Ratio(data = data)
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
    fitX1toY2 <- summary(lm(Y2 ~ X1, data = data))
    fitX1toY1 <- summary(lm(Y1 ~ X1, data = data))
    fitX2toY1 <- summary(lm(Y1 ~ X2, data = data))
    fitX2toY2 <- summary(lm(Y2 ~ X2, data = data))
    est_g13 <- fitX1toY2$coefficients[-1, 1]
    est_b11 <- fitX1toY1$coefficients[-1, 1]
    sd_g13 <- fitX1toY2$coefficients[-1, 2]
    sd_b11 <- fitX1toY1$coefficients[-1, 2]
    g12_ivw <- mr_ivw(mr_input(bx = est_b11, bxse = sd_b11, by = est_g13,
        byse = sd_g13))
    est_g12 <- g12_ivw$Estimate
    est_g23 <- fitX2toY1$coefficients[-1, 1]
    est_b22 <- fitX2toY2$coefficients[-1, 1]
    sd_g23 <- fitX2toY1$coefficients[-1, 2]
    sd_b22 <- fitX2toY2$coefficients[-1, 2]
    g21_ivw <- mr_ivw(mr_input(bx = est_b22, bxse = sd_b22, by = est_g23,
        byse = sd_g23))
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
