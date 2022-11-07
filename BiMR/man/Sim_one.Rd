\name{Sim_one}
\alias{Sim_one}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sim_one
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Sim_one(method = "Ratio", GRSorIV = "IVs", seed_random = 1, g12 = 0.3, g21 = 0.3, nIVs = 0, nIVw = 0, con1_var = 1, Ve1 = 0.1, Ve2 = 0.1, bcy1 = 1, bcy2 = 1, datatype = "bi_infi", n_sim = 1000, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{GRSorIV}{
%%     ~~Describe \code{GRSorIV} here~~
}
  \item{seed_random}{
%%     ~~Describe \code{seed_random} here~~
}
  \item{g12}{
%%     ~~Describe \code{g12} here~~
}
  \item{g21}{
%%     ~~Describe \code{g21} here~~
}
  \item{nIVs}{
%%     ~~Describe \code{nIVs} here~~
}
  \item{nIVw}{
%%     ~~Describe \code{nIVw} here~~
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
  \item{datatype}{
%%     ~~Describe \code{datatype} here~~
}
  \item{n_sim}{
%%     ~~Describe \code{n_sim} here~~
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
function (method = "Ratio", GRSorIV = "IVs", seed_random = 1,
    g12 = 0.3, g21 = 0.3, nIVs = 0, nIVw = 0, con1_var = 1, Ve1 = 0.1,
    Ve2 = 0.1, bcy1 = 1, bcy2 = 1, datatype = "bi_infi", n_sim = 1000,
    ...)
{
    print(c(g12, g21))
    print(GRSorIV)
    sim <- seq(n_sim * 1.1)
    temp <- sapply(sim, FUN = function(x) Data(g12 = g12, g21 = g21,
        nIVs = nIVs, nIVw = nIVw, con1_var = con1_var, Ve1 = Ve1,
        Ve2 = Ve2, bcy1 = bcy1, bcy2 = bcy2, Causal = datatype,
        seed = x + seed_random * n_sim))
    if (method == "all") {
        Est <- sapply(sim, FUN = function(x) Est_all(data = temp[,
            x], GRSorIV = GRSorIV))
        Est_slt <- Est[, colSums(is.na(Est)) == 0][, 1:n_sim]
        abs_bias <- data.frame(matrix(NA, n_sim, 8))
        colnames(abs_bias) = c("Ratio.G12", "Ratio.G21", "BiRatio.G12",
            "BiRatio.G21", "Liml.G12", "Liml.G21", "BiLiml.G12",
            "BiLiml.G21")
        result1 <- matrix(NA, 8, 6)
        result2 <- matrix(NA, 8, 5)
        colnames(result1) <- c("Mean.abs_bias", "Median.abs_bias",
            "MSE", "RMSE", "Mean.SD", "Coverage")
        colnames(result2) <- c("Median.est", "Bias.median", "Mean.est",
            "Bias.mean", "SD.bootstrap")
        rownames(result1) <- rownames(result2) <- c("Ratio.g12",
            "Ratio.g21", "BiRatio.g12", "BiRatio.g21", "Liml.g12",
            "Liml.g21", "BiLiml.g12", "BiLiml.g21")
        CI <- function(vec) {
            return(stats::quantile(vec, probs = c(0.025, 0.975),
                na.rm = T))
        }
        for (i in seq(4)) {
            abs_bias.g12 <- abs(Est_slt[i * 8 - 7, ] - g12)
            abs_bias.g21 <- abs(Est_slt[i * 8 - 3, ] - g21)
            abs_bias[, i * 2 - 1] <- abs_bias.g12
            abs_bias[, i * 2] <- abs_bias.g21
            MSEg12 <- mean((Est_slt[i * 8 - 7, ] - g12)^2)
            MSEg21 <- mean((Est_slt[i * 8 - 3, ] - g21)^2)
            mean.abs.bias.g12 <- paste0(round(mean(abs_bias.g12),
                2), "(", 100 * round(mean(abs_bias.g12)/abs(g12),
                2), "\%)")
            median.abs.bias.g12 <- paste0(round(median(abs_bias.g12),
                2), "(", 100 * round(median(abs_bias.g12)/abs(g12),
                2), "\%)")
            mean.abs.bias.g21 <- paste0(round(mean(abs_bias.g21),
                2), "(", 100 * round(mean(abs_bias.g21)/abs(g21),
                2), "\%)")
            median.abs.bias.g21 <- paste0(round(median(abs_bias.g21),
                2), "(", 100 * round(median(abs_bias.g21)/abs(g21),
                2), "\%)")
            result1[(i * 2 - 1):(i * 2), ] <- t(matrix(c(mean.abs.bias.g12,
                median.abs.bias.g12, round(MSEg12, 2), round(sqrt(MSEg12),
                  2), round(mean(Est_slt[i * 8 - 6, ]), 2), round(mean(Est_slt[i *
                  8 - 5, ] <= g12 & Est_slt[i * 8 - 4, ] >= g12),
                  2), mean.abs.bias.g21, median.abs.bias.g21,
                round(MSEg21, 2), round(sqrt(MSEg21), 2), round(mean(Est_slt[i *
                  8 - 2, ]), 2), round(mean(Est_slt[i * 8 - 1,
                  ] <= g21 & Est_slt[i * 8, ] >= g21), 2)), 6))
            result2[(i * 2 - 1):(i * 2), ] <- t(matrix(c(median(Est_slt[i *
                8 - 7, ]), median(Est_slt[i * 8 - 7, ]) - g12,
                mean(Est_slt[i * 8 - 7, ]), mean(Est_slt[i *
                  8 - 7, ]) - g12, sd(Est_slt[i * 8 - 7, ]),
                median(Est_slt[i * 8 - 3, ]), median(Est_slt[i *
                  8 - 3, ]) - g21, mean(Est_slt[i * 8 - 3, ]),
                mean(Est_slt[i * 8 - 3, ]) - g21, sd(Est_slt[i *
                  8 - 3, ])), 5))
        }
        compare <- c(cmp.Ratio.G12 = mean(abs_bias$Ratio.G12 -
            abs_bias$BiRatio.G12 >= 0), cmp.Ratio.G21 = mean(abs_bias$Ratio.G21 -
            abs_bias$BiRatio.G21 >= 0), cmp.Liml.G12 = mean(abs_bias$Liml.G12 -
            abs_bias$BiLiml.G12 >= 0), cmp.Liml.G21 = mean(abs_bias$Liml.G21 -
            abs_bias$BiLiml.G21 >= 0))
        results <- list(result1 = result1, result2 = result2,
            Est_slt = Est_slt, compare = compare)
    }
    else {
        Est <- sapply(sim, FUN = function(x) Est_bi_IVs(data = temp[,
            x], method = method, GRSorIV = GRSorIV))
        Est_slt <- Est[, colSums(is.na(Est)) == 0][, 1:n_sim]
        CI <- function(vec) {
            return(stats::quantile(vec, probs = c(0.025, 0.975),
                na.rm = T))
        }
        i <- 1
        abs_bias.g12 <- abs(Est_slt[i * 8 - 7, ] - g12)
        abs_bias.g21 <- abs(Est_slt[i * 8 - 3, ] - g21)
        MSEg12 <- mean((Est_slt[i * 8 - 7, ] - g12)^2)
        MSEg21 <- mean((Est_slt[i * 8 - 3, ] - g21)^2)
        result1 <- t(matrix(c(mean(abs_bias.g12), median(abs_bias.g12),
            MSEg12, sqrt(MSEg12), mean(Est_slt[i * 8 - 6, ]),
            mean(Est_slt[i * 8 - 5, ] <= g12 & Est_slt[i * 8 -
                4, ] >= g12), mean(abs_bias.g21), median(abs_bias.g21),
            MSEg21, sqrt(MSEg21), mean(Est_slt[i * 8 - 2, ]),
            mean(Est_slt[i * 8 - 1, ] <= g21 & Est_slt[i * 8,
                ] >= g21)), 6))
        result2 <- t(matrix(c(median(Est_slt[i * 8 - 7, ]), median(Est_slt[i *
            8 - 7, ]) - g12, mean(Est_slt[i * 8 - 7, ]), mean(Est_slt[i *
            8 - 7, ]) - g12, sd(Est_slt[i * 8 - 7, ]), median(Est_slt[i *
            8 - 3, ]), median(Est_slt[i * 8 - 3, ]) - g21, mean(Est_slt[i *
            8 - 3, ]), mean(Est_slt[i * 8 - 3, ]) - g21, sd(Est_slt[i *
            8 - 3, ])), 5))
        colnames(result1) <- c("Mean.abs_bias", "Median.abs_bias",
            "MSE", "RMSE", "Mean.SD", "Coverage")
        colnames(result2) <- c("Median.est", "Bias.median", "Mean.est",
            "Bias.mean", "SD.bootstrap")
        rownames(result1) <- rownames(result2) <- c("g12", "g21")
        results <- list(result1 = result1, result2 = result2,
            Est_slt = Est_slt)
    }
    return(results)
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
