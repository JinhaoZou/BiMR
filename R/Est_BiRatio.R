#' Function to estimation causal effects using bidirectional Ratio method
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
#' Est_BiRatio(data = example_data)
#' }

library(MendelianRandomization)

Est_BiRatio <- function(data = data){
  nX1 <- dim(as.matrix(data$X1))[2]
  nX2 <- dim(as.matrix(data$X2))[2]
  #nC <- dim(as.matrix(data$con))[2]

  data_variable <-  cbind(data$X1, data$X2)
  colnames(data_variable) <- c(sapply(seq(nX1), FUN = function(x) paste("X1.", x, sep = "")),
                               sapply(seq( nX2), FUN = function(x) paste("X2.", x, sep = "")))

  eqY1Y2 <- Y2 ~. - Y1
  eqY2Y1 <- Y1 ~. - Y2
  system <- list(Y1Y2 = eqY1Y2, Y2Y1 = eqY2Y1)

  data_all <- as.data.frame(cbind(data$Y1, data$Y2, data_variable))
  colnames(data_all)[1:2]<- c("Y1", "Y2")

  fit1 <- summary(lm(eqY1Y2, data = data_all))
  fit2 <- summary(lm(eqY2Y1, data = data_all))

  coef1 <- as.data.frame(fit1$coefficients)
  est_Y2X1_g12b11 <- coef1$Estimate[2:(1+nX1)]
  sd_Y2X1_g12b11 <- coef1$`Std. Error`[2:(1+nX1)]
  est_Y2X2_b22 <- coef1$Estimate[(2+nX1):(1+nX1+nX2)]
  sd_Y2X2_b22 <- coef1$`Std. Error`[(2+nX1):(1+nX1+nX2)]

  coef2 <- as.data.frame(fit2$coefficients)
  est_Y1X1_b11 <- coef2$Estimate[2:(1+nX1)]
  sd_Y1X1_b11 <- coef2$`Std. Error`[2:(1+nX1)]
  est_Y1X2_g21b22 <- coef2$Estimate[(2+nX1):(1+nX1+nX2)]
  sd_Y1X2_g21b22 <- coef2$`Std. Error`[(2+nX1):(1+nX1+nX2)]


  g12_ivw <- mr_ivw(mr_input(bx = est_Y1X1_b11, bxse = sd_Y1X1_b11, by = est_Y2X1_g12b11, byse = sd_Y2X1_g12b11))
  g21_ivw <- mr_ivw(mr_input(bx = est_Y2X2_b22, bxse = sd_Y2X2_b22, by = est_Y1X2_g21b22, byse = sd_Y1X2_g21b22))

  est_g12 <- g12_ivw$Estimate
  est_g21 <- g21_ivw$Estimate

  est_g12_result <- c(g12 = est_g12, Sdg12 = g12_ivw$StdError, CIlowg12 = g12_ivw$CILower,  CIUpg12 = g12_ivw$CIUpper)
  est_g21_result <- c(g21 = est_g21, Sdg21 = g21_ivw$StdError, CIlowg21 = g21_ivw$CILower,  CIUpg21 = g21_ivw$CIUpper)

  return(c(est_g12_result, est_g21_result))
}
