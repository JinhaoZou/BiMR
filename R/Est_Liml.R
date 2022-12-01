#' Function to estimation causal effects using LIML method
#'
#' @param data data is a list inludes X1, X2, Y1, Y2, confounder, F statistics of X1 on Y1, and F statistics of X2 on Y2
#' The X1 and X2 is matrix with number of individuals * number of IVs, Y1, Y2, confounder are vector with length bumber of individuals
#'
#' @return a vector with length 8, includes: estimation of g12, sd of g12 estimation, CI of g12, estimation of g21, sd of g21 estimation, CI of g21
#' g12 is the causal effect of Y1 on Y2, g21 is the causal effect of Y2 on Y1
#' @export
#' @examples {
#' data(example_data)
#' attach(example_data)
#' Est_Liml(data = example_data)
#' }

library(ivmodel)

Est_Liml <- function(data = data){
  #print("liml")
  model12 <- ivmodel(Y = data$Y2, D = data$Y1, Z = data$X1)
  est_g12_model <- model12$LIML
  est_g12 <- est_g12_model$point.est

  model21 <- ivmodel(Y = data$Y1, D = data$Y2, Z = data$X2)
  est_g21_model <- model21$LIML
  est_g21 <- est_g21_model$point.est

  est_g12_result <- c(g12 = est_g12_model$point.est, Sdg12 = est_g12_model$std.err, CIlowg12 = est_g12_model$ci[1],  CIUpg12 = est_g12_model$ci[2])
  est_g21_result <- c(g21 = est_g21_model$point.est, Sdg21 = est_g21_model$std.err, CIlowg21 = est_g21_model$ci[1], CIUpg21 = est_g21_model$ci[2])

  return(c(est_g12_result, est_g21_result))
}
