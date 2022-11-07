#' Function to simulate multiple datasets and evaluate the estimation bias of different methods
#'  Model:
#'  UMR:  Y1 = b01 + X1B11 + bcy1con1 + e1; Y2 = b02 + X2B22 + g12Y1 + bcy2con1 + e2
#'  BMR:  Y1 = b01 + X1B11 + g21Y2 + bcy1con1 + e1; Y2 = b02 + X2B22 + g12Y1 + bcy2con1 + e2
#'  X1 and X2 are IVs for Y1 and Y2, Y1 and Y2 are interested phenotype
#'
#' @param seed_random seed of geenrating dataset; default: 1
#' @param n_sim number replicates of this simulation; default: 1000
#' @param g12 causal effect of Y1 on Y2; default: 0.3
#' @param g21 causal effect of Y2 on Y1; default: 0.3
#' @param b11 effect of IVs on Y1; default: 0.1
#' @param b22 effect of IVs on Y2; default: 0.1
#' @param nIV number of IVs used for Y1/Y2; default: 10
#' @param Ve1 the variance of random error for Y1, generating with normal distribution with mean 0; default: 0.1
#' @param Ve2 the variance of random error for Y2, generating with normal distribution with mean 0; default: 0.1
#' @param con1_var the variance of generating the confounder, generating with normal distribution with mean 1; default: 1
#' @param bcy1 effect of confounder on Y1; default: 0.3
#' @param bcy2 effect of coufounder on Y2; default: 0.3
#' @param causal the causal effects scenarios choose one from "uni", "bi_infi", "bi_fl";
#' @param method the method used for estimation, if it is all, four methods are all used for each simulated dataset,
#' choose one from "Ratio", "BiRatio", "Liml", "BiLiml", "all"
#'
#' @return returns three results (result1, result2, result3) including important matrics to evaluate the bias of estimations using different methods,
#' and all the estimations (Est_slt);If the estimation using multiple methods, results include the camparison of unidirectional method and bidirectional methods (compare).
#' @export
#' @examples {
#' Sim_one(n_sim = 5, causal = "bi_infi", method = "Ratio")
#' }

library(ivmodel)
library(MendelianRandomization)

Sim_one <- function(seed_start = 0, seed_random = 1, n_sim = 1000, g12 = 0.3, g21 = 0.3,
                    b11 = 0.1, b22 = 0.1, nIV = 10,
                    Ve1 = 0.1, Ve2= 0.1, con1_var = 1, bcy1 = 0.3, bcy2 = 0.3,
                    causal = c("uni","bi_infi","bi_fl"), method = c("Ratio", "BiRatio", "Liml", "BiLiml", "all"),...){
  # method = "all"
  # causal = "bi_infi"
  # seed_random = 364 #to make sure different secnarios(different g12 and g21) has different seeds for data
  # n_sim = 5
  # g12 = 1.7
  # g21 = -1.3
  #
  # b11 = 1
  # b22 = 1
  # nIV = 1
  #
  # Ve1 = 0.1
  # Ve2 = 0.1
  #
  # bcy1 = 0.3
  # bcy2 = 0.3
  # con1_var = 1
  # Method <- c("Ratio",  "BiRatio", "Liml", "BiLiml")
  # Causal <- c("uni", "bi_infi", "bi_fl")

  print(seed_random)
  print(c(g12, g21)) # print the causal effects value to make sure we have correct setting

  ############ Generate the datasets ##################################
  sim <- seq(n_sim*1.1)
  temp <- sapply(sim, FUN = function(x) Data(g12 = g12, g21 = g21, Ve1 = Ve1, Ve2= Ve2,
                                             b11 = b11, b22 = b22, nIV = nIV,
                                             con1_var = con1_var, bcy1 = bcy1, bcy2 = bcy2,
                                             causal = causal, seed = x + seed_random*n_sim*1.1 + seed_start))

  ######## Estimate the mean F-statistics ##############################################
  Fst.X1 <- mean(sapply(sim, FUN = function(x) temp[,x]$Fst[3]))
  Fst.X2 <- mean(sapply(sim, FUN = function(x) temp[,x]$Fst[4]))


  ######## Function to extract most interest metrics ###################################
  Get_metrics <- function(i, Est_slt, g12, g21){

    abs_bias.g12 <- abs(Est_slt[i*8-7,] - g12) #absolute bias
    abs_bias.g21 <- abs(Est_slt[i*8-3,] - g21)

    MSEg12 <- mean((Est_slt[i*8-7,] - g12)^2)
    MSEg21 <- mean((Est_slt[i*8-3,] - g21)^2)

    Mean.abs_bias.g12 <- round(mean(abs_bias.g12), 2)
    Pct.mean.abs_bias.g12 <- paste0(100*round(Mean.abs_bias.g12 / abs(g12), 2), "%")
    Mean.abs_bias.g21 <- round(mean(abs_bias.g21), 2)
    Pct.mean.abs_bias.g21 <- paste0(100*round(Mean.abs_bias.g21 / abs(g21), 2), "%")

    Median.abs_bias.g12 <- round(median(abs_bias.g12), 2)
    Pct.median.abs_bias.g12 <- paste0(100*round(Median.abs_bias.g12 / abs(g12), 2), "%")
    Median.abs_bias.g21 <- round(median(abs_bias.g21), 2)
    Pct.median.abs_bias.g21 <- paste0(100*round(Median.abs_bias.g21 / abs(g21), 2), "%")


    result1 <- t(matrix(c(Mean.abs_bias.g12, Pct.mean.abs_bias.g12, Median.abs_bias.g12, Pct.median.abs_bias.g12,
                          round(MSEg12,2), round(sqrt(MSEg12),2), round(mean(Est_slt[i*8-6, ]),2),
                          round(mean(Est_slt[i*8-5, ] <= g12 & Est_slt[i*8-4, ] >= g12),2),
                          Mean.abs_bias.g21, Pct.mean.abs_bias.g21, Median.abs_bias.g21, Pct.median.abs_bias.g21,
                          round(MSEg21,2), round(sqrt(MSEg21),2), round(mean(Est_slt[i*8-2, ]),2),
                          round(mean(Est_slt[i*8-1, ] <= g21 & Est_slt[i*8, ] >= g21),2)),
                        8))

    result2 <- t(matrix(c(median(Est_slt[i*8-7, ]), median(Est_slt[i*8-7, ]) - g12, mean(Est_slt[i*8-7, ]), mean(Est_slt[i*8-7, ])- g12, sd(Est_slt[i*8-7, ]),
                          median(Est_slt[i*8-3, ]), median(Est_slt[i*8-3, ]) - g21, mean(Est_slt[i*8-3, ]), mean(Est_slt[i*8-3, ])- g21, sd(Est_slt[i*8-3, ])),5))


    result3 <- t(matrix(c(round(median(Est_slt[i*8-7, ]), 2), Median.abs_bias.g12, Pct.median.abs_bias.g12,
                          round(median(Est_slt[i*8-3, ]), 2), Median.abs_bias.g21, Pct.median.abs_bias.g21),3))


    return(list(result1 = result1, result2 = result2, result3 = result3,
                abs_bias.g12 = abs_bias.g12, abs_bias.g21 = abs_bias.g21))

  }


  ###### Estimate the bidirectional causal effects based on method ###############################
  Est <- sapply(sim, FUN = function(x) Est_all(data = temp[,x]$data, method = method))

  # If there are to many NAs, which means the value g12 and g21 can not converge in bidirectional scenarios
  # We need to set all output as NA is this scenario happens
  tryCatch({
    Est_slt <- Est[, colSums(is.na(Est)) == 0][,1:n_sim]
    # if the method is not all, the dimention of Est_slt is 8*n_sim
    # if the method is all, the dimention of Est_slt is 32*n_sim.
    ## the 32 means Est_G12, Sd_G12, G12CIup, G12CIdown, EstG21, Sd_G21, G21CIup, G21down, for each method.

    # Simulation, estimation with all methods
    if(method == "all"){
      # collect information of absolute bias of each method
      abs_bias <- data.frame(matrix(NA, n_sim, 8))

      # Extract the information of all estimation bias, sd of bias, quantile of bias.
      ## All possibaly cared parameters
      result1 <- matrix(NA, 8, 8)
      result2 <- matrix(NA, 8, 5)
      ## The care most parameters
      result3 <- matrix(NA, 8, 3)

      for(i in seq(4)){
        #each i means one method
        metrics <- Get_metrics(i, Est_slt, g12, g21)
        result1[(i*2-1):(i*2),] <- metrics$result1
        result2[(i*2-1):(i*2),] <- metrics$result2
        result3[(i*2-1):(i*2),] <- metrics$result3


        abs_bias[,i*2-1] <- metrics$abs_bias.g12
        abs_bias[,i*2] <- metrics$abs_bias.g21
      }

      colnames(result1) <- c("Mean.abs_bias", "Pct.mean.abs_bias", "Median.abs_bias", "Pct.median.abs_bias", "MSE", "RMSE","Mean.SD", "Coverage")
      colnames(result2) <- c("Median.est", "Bias.median", "Mean.est", "Bias.mean", "SD.bootstrap")
      colnames(result3) <- c("Median.est", "Median.abs_bias", "Pct.median.abs_bias")
      rownames(result1) <- rownames(result2) <- rownames(result3) <- c("Ratio.g12", "Ratio.g21", "BiRatio.g12", "BiRatio.g21",
                                                                       "Liml.g12", "Liml.g21", "BiLiml.g12", "BiLiml.g21")

      # result of comparing Ratio and BiRatio methods, LIML and BiLIML methos
      colnames(abs_bias) = c("Ratio.G12", "Ratio.G21", "BiRatio.G12", "BiRatio.G21",
                             "Liml.G12", "Liml.G21", "BiLiml.G12", "BiLiml.G21")
      compare <- c(cmp.Ratio.G12 = mean(abs_bias$Ratio.G12 - abs_bias$BiRatio.G12 >= 0),
                   cmp.Ratio.G21 = mean(abs_bias$Ratio.G21 - abs_bias$BiRatio.G21 >= 0),
                   cmp.Liml.G12 = mean(abs_bias$Liml.G12 - abs_bias$BiLiml.G12 >= 0),
                   cmp.Liml.G21 = mean(abs_bias$Liml.G21 - abs_bias$BiLiml.G21 >= 0))

      results <- list(result1 = result1, result2 = result2, result3 = result3,
                      Est_slt = Est_slt, Fst = c(Fst.X1, Fst.X2), compare = compare)

      # simulate with one method
    }else{

      metrics <- Get_metrics(1, Est_slt, g12, g21)
      result1 <- metrics$result1
      result2 <- metrics$result2
      result3 <- metrics$result3

      colnames(result1) <- c("Mean.abs_bias", "Pct.mean.abs_bias", "Median.abs_bias", "Pct.median.abs_bias", "MSE", "RMSE","Mean.SD", "Coverage")
      colnames(result2) <- c("Median.est","Bias.median", "Mean.est", "Bias.mean", "SD.bootstrap")
      colnames(result3) <- c("Median.est", "Median.abs_bias", "Pct.median.abs_bias")

      rownames(result1) <- rownames(result2) <- rownames(result3) <- c("g12", "g21")

      results <- list(result1 = result1, result2 = result2, result3 = result3,
                      Est_slt = Est_slt, Fst = c(Fst.X1, Fst.X2), Compare = NA)
    }

    return(results)
  }, error = function(e)
    return(list(result1 = NA, result2 = NA, result3 = NA,
                    Est_slt = NA, Fst = c(Fst.X1, Fst.X2), Compare = NA ))
  )

}


##Add the compare = NA at the results for one method. otherwise the result for trycatch result and one method result will not be consistent


