#' Function to simulate all scenarios of g12 and g21 with a range from -2 to 2
#'  Model:
#'  UMR:  Y1 = b01 + X1B11 + bcy1con1 + e1; Y2 = b02 + X2B22 + g12Y1 + bcy2con1 + e2
#'  BMR:  Y1 = b01 + X1B11 + g21Y2 + bcy1con1 + e1; Y2 = b02 + X2B22 + g12Y1 + bcy2con1 + e2
#'  X1 and X2 are IVs for Y1 and Y2, Y1 and Y2 are interested phenotype
#'
#' @param n_sim number replicates of this simulation; default: 1000
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
#' sim_all(n_sim = 5, causal = "bi_infi", method = "Ratio")
#' }

sim_all <- function(n_sim = 1000, seed_start_i = 0,
                    b11 = 0.1, b22 = 0.1, nIV = 10, #Normally we will change value here
                    Ve1 = 0.1, Ve2= 0.1, con1_var = 1, bcy1 = 0.3, bcy2 = 0.3, #we perfer keep this consistant
                    causal = c("uni","bi_infi","bi_fl"), method = c("Ratio", "BiRatio", "Liml", "BiLiml", "all"),...){
  # n_sim = 5
  # b11 = 0.1
  # b22 = 0.1
  # nIV = 10
  # Ve1 = 0.1
  # Ve2 = 0.1
  # con1_var = 1
  # bcy1 = 0.3
  # bcy2 = 0.3
  # causal = "bi_infi"
  # method = "Ratio"

  set.seed(321)
  G12 <- G21 <- seq(1,20)*0.2 - 2.1
  G12G21 <- cbind(rep(G12, each = 20), rep(G21,20))

  #e.g. Separate the 1000 simulations to 10 parts, each have different value
  #seed_start

  if(causal == "uni"){
    seed_start <- seed_start_i*21*n_sim*1.1
    a <- sapply(seq(length(G12)), FUN = function(x)
      Sim_one(seed_start = seed_start, seed_random = x, n_sim = n_sim, g12 = G12[x], g21 = 0,
              b11 = b11, b22 = b22, nIV = nIV,
              Ve1 = Ve1, Ve2 = Ve2, con1_var = con1_var, bcy1 = bcy1, bcy2 = bcy2,
              causal= causal, method = method))

    }else{
    seed_start <- seed_start_i*401*n_sim*1.1
    a <- sapply(seq(dim(G12G21)[1]), FUN = function(x)
      Sim_one(seed_start = seed_start, seed_random = x, n_sim = n_sim, g12 = G12G21[x,1], g21 = G12G21[x,2],
              b11 = b11, b22 = b22, nIV = nIV,
              Ve1 = Ve1, Ve2 = Ve2, con1_var = con1_var, bcy1 = bcy1, bcy2 = bcy2,
              causal = causal, method = method))
    }

  return(a)
}
