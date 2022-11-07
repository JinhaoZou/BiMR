\name{Est_Liml}
\alias{Est_Liml}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Est_Liml
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Est_Liml(data = data)
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
    model12 <- ivmodel(Y = data$Y2, D = data$Y1, Z = data$X1)
    est_g12_model <- model12$LIML
    est_g12 <- est_g12_model$point.est
    model21 <- ivmodel(Y = data$Y1, D = data$Y2, Z = data$X2)
    est_g21_model <- model21$LIML
    est_g21 <- est_g21_model$point.est
    est_g12_result <- c(g12 = est_g12_model$point.est, Sdg12 = est_g12_model$std.err,
        CIlowg12 = est_g12_model$ci[1], CIUpg12 = est_g12_model$ci[2])
    est_g21_result <- c(g21 = est_g21_model$point.est, Sdg21 = est_g21_model$std.err,
        CIlowg21 = est_g21_model$ci[1], CIUpg21 = est_g21_model$ci[2])
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
