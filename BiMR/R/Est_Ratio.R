#' Function to estimation causal effects using Ratio method
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
#' Est_Ratio(data = example_data)
#' }

library(MendelianRandomization)

Est_Ratio <- function(data = data){
  #res1 <- residuals(summary(lm(Y1 ~ con, data = data, na.action = "na.exclude")))
  #res2 <- residuals(summary(lm(Y2 ~ con, data = data, na.action = "na.exclude")))

  fitX1toY2 <- summary(lm(Y2 ~ X1 , data = data))
  fitX1toY1 <- summary(lm(Y1 ~ X1 , data = data))

  fitX2toY1 <- summary(lm(Y1 ~ X2 , data = data))
  fitX2toY2 <- summary(lm(Y2 ~ X2 , data = data))

  est_g13 <- fitX1toY2$coefficients[-1,1]
  est_b11 <- fitX1toY1$coefficients[-1,1]
  sd_g13 <- fitX1toY2$coefficients[-1, 2]
  sd_b11 <- fitX1toY1$coefficients[-1, 2]
  g12_ivw <- mr_ivw(mr_input(bx = est_b11, bxse = sd_b11, by = est_g13, byse = sd_g13))
  est_g12 <- g12_ivw$Estimate

  #est_g12 <- sum(fitX1toY1$coefficients[-1,1]*fitX1toY2$coefficients[-1,1]/sd_g13^2)/sum((fitX1toY1$coefficients[-1,1]/sd_g13)^2)

  est_g23 <- fitX2toY1$coefficients[-1,1]
  est_b22 <- fitX2toY2$coefficients[-1,1]
  sd_g23 <- fitX2toY1$coefficients[-1, 2]
  sd_b22 <- fitX2toY2$coefficients[-1, 2]

  g21_ivw <- mr_ivw(mr_input(bx = est_b22, bxse = sd_b22, by = est_g23, byse = sd_g23))
  est_g21 <- g21_ivw$Estimate
  #est_g21 <- sum(fitX2toY2$coefficients[-1,1]*fitX2toY1$coefficients[-1,1]/sd_g23^2)/sum((fitX2toY2$coefficients[-1,1]/sd_g23)^2)

  est_g12_result <- c(g12 = est_g12, Sdg12 = g12_ivw$StdError, CIlowg12 = g12_ivw$CILower,  CIUpg12 = g12_ivw$CIUpper)
  est_g21_result <- c(g21 = est_g21, Sdg21 = g21_ivw$StdError, CIlowg21 = g21_ivw$CILower,  CIUpg21 = g21_ivw$CIUpper)

  return(c(est_g12_result, est_g21_result))
}
