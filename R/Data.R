#' Function to simulate one dataset based on the direction of causal effects
#'  UMR:  Y1 = b01 + X1B11 + bcy1con1 + e1; Y2 = b02 + X2B22 + g12Y1 + bcy2con1 + e2
#'  BMR:  Y1 = b01 + X1B11 + g21Y2 + bcy1con1 + e1; Y2 = b02 + X2B22 + g12Y1 + bcy2con1 + e2
#'  X1 and X2 are IVs for Y1 and Y2, Y1 and Y2 are interested phenotype
#' @param seed seed of geenrating dataset; default: 321
#' @param n number of individual in the generated dataset; default: 1000
#' @param b01 baseline Y1; default: 1
#' @param b02 baseline Y2; default: 1
#' @param p1  Minor allele frequency of X1; default: 0.3
#' @param p2  Minor allele frequency of X2; default: 0.3
#' @param b11 effect of IVs on Y1; default: 0.1
#' @param b22 effect of IVs on Y2; default: 0.1
#' @param nIV number of IVs used for Y1/Y2; default: 10
#' @param g12 causal effect of Y1 on Y2; default: 0.3
#' @param g21 causal effect of Y2 on Y1; default: 0.3
#' @param con1_var the variance of generating the confounder, generating with normal distribution with mean 1; default: 1
#' @param Ve1 the variance of random error for Y1, generating with normal distribution with mean 0; default: 0.1
#' @param Ve2 the variance of random error for Y2, generating with normal distribution with mean 0; default: 0.1
#' @param bcy1 effect of confounder on Y1; default: 0.3
#' @param bcy2 effect of coufounder on Y2; default: 0.3
#' @param causal the causal effects scenarios choose one from "uni", "bi_infi", "bi_fl";
#'
#' @return data is a list inludes X1 (X1), X2 (X2), Y1 (Y1), Y2 (Y2), confounder (con), F statistics of X1 on Y1 (Fst_X1), and F statistics of X2 on Y2 (Fst_X2)
#' The X1 and X2 is matrix with number of individuals * number of IVs, Y1, Y2, confounder are vector with length n
#' @export
#' @examples {
#' Data(causal = "bi_infi")
#' }

Data <- function(seed = 321, n = 1000, b01 = 1, b02 = 1, p1 = 0.3, p2 = 0.3,
                 b11 = 0.1, b22 = 0.1, nIV = 10,  g12 = 0.3, g21 = 0.3,
                 con1_var = 1, Ve1 = 0.1, Ve2= 0.1, bcy1 = 0.3, bcy2 = 0.3, causal= c("uni", "bi_infi", "bi_fl"),...){

  # if the causal effect is Unidirectional, the causal effect of Y2 on Y1 is 0
  if(causal == "uni"){
    g21 <- 0
  }

  #print(seed)
  #seed = 19
  set.seed(seed)
  e1 <- rnorm(n,0,Ve1)
  e2 <- rnorm(n,0,Ve2)

  B11 = rep(b11,nIV)
  B22 = rep(b22,nIV)

  con1 <- abs(rnorm(n, 1, con1_var))
  #con2 <- abs(rnorm(n, 0, 0.01))

  X1 <- sapply(abs(rnorm(nIV, p1, 0.01)), FUN = function(x) rbinom(n, 2, x))
  X2 <- sapply(abs(rnorm(nIV, p2, 0.01)), FUN = function(x) rbinom(n, 2, x))

  #X1 <- matrix(sample(seq(2)-1, nIV*n, replace = T),n) #Sampling method in paper Burgess.2010
  #X2 <- matrix(sample(seq(2)-1, nIV*n, replace = T),n)
  #X1 <- matrix(rnorm(n*nIV, 0, 5),n)
  #X2 <- matrix(rnorm(n*nIV, 0, 5),n)

  if(causal == "uni"){
    #print(X1)
    #print(B11)
    Y1 <- b01 + X1%*%B11 + bcy1*con1 + e1
    Y2 <- b02 + g12*Y1 + bcy2*con1 + X2%*%B22 + e2
    #Y2 <- b02 + g12*Y1 + bcy2*con1 + e2
  }else if(causal == "bi_infi"){
    Y1 <- (b01 + g21 * b02 + g21 * bcy2 * con1 + X2 %*% B22 * g21 + g21 * e2 + bcy1 * con1 +  X1 %*% B11  + e1) / (1 - g21 * g12)
    Y2 <- (b02 + g12 * b01 + g12 * bcy1 * con1 + X1 %*% B11 * g12 + g12 * e1 + bcy2 * con1 +  X2 %*% B22  + e2) / (1 - g21 * g12)
  }else{
    iter_sta <- sample(2:21, n, replace=T)
    #iter <- 15.5
    #iter_sta <- rep(iter, n)


    #make sure first half starts from Y1, second half starts from Y2. 0 means starting from Y2, 1 means starting from Y1
    #Y_chos <- rep(c(0,1), each = n/2)
    Y_chos <- rbinom(n,1,0.5)

    #Y1_all <- rnorm(n,22,1)
    #Y2_all <- rnorm(n,22,1)

    ## Y For every person
    Gen_Y <- function(i){
      #iterative status of different people
      iter_sta_one <- iter_sta[i]
      Y_chos_one <- Y_chos[i]

      #Y_start <- ifelse(Y_chos_one == 0, Y1_one, Y2_one)

      X1_one <- X1[i,]
      X2_one <- X2[i,]
      con1_one <- con1[i]
      #con2_one <- con2[i]
      e1_one <- e1[i]
      e2_one <- e2[i]

      #Y1_one <- Y1_all[i]
      #Y2_one <- Y2_all[i]
      Y1_one <- rnorm(1,b01,0.01)
      Y2_one <- rnorm(1,b02,0.01)

      if(Y_chos_one == 0){
        for(j in 1:iter_sta_one ){
          Y1_one_pre <- Y1_one
          Y2_one <- b02 + B22%*%X2_one + bcy2*con1_one + g12*Y1_one_pre + e2_one
          Y2_one_pre <- Y2_one
          Y1_one <- b01 + B11%*%X1_one + bcy1*con1_one + g21*Y2_one_pre + e1_one

        }
      }else{
        for(j in 1:iter_sta_one){
          Y2_one_pre <- Y2_one
          Y1_one <- b01 + B11%*%X1_one + bcy1*con1_one + g21*Y2_one_pre + e1_one
          Y1_one_pre <- Y1_one
          Y2_one <- b02 + B22%*%X2_one + bcy2*con1_one + g12*Y1_one_pre + e2_one
        }
      }

      return(c(Y1_one = Y1_one, Y2_one = Y2_one, Y1_one_pre = Y1_one_pre, Y2_one_pre = Y2_one_pre))
    }


    Y <- sapply(seq(n), FUN = function(x) Gen_Y(x))
    Y1 <- Y[1, ]
    Y2 <- Y[2, ]
    #con <- Y[3,]
    Y1_pre <- Y[3,]
    Y2_pre <- Y[4,]

    #VE1 <- (g21^2*Ve2 + Ve1)/(1-g21*g12)^2
    #VE2 <- (g12^2*Ve1 + Ve2)/(1-g21*g12)^2
  }

  Fst_X1 <- summary(lm(Y1 ~ X1))$fstatistic[[1]]
  Fst_X2 <- summary(lm(Y2 ~ X2))$fstatistic[[1]]

  Fst_X1.1 <- summary(lm(Y1 ~ X1[,1]))$fstatistic[[1]]
  Fst_X2.1 <- summary(lm(Y2 ~ X2[,1]))$fstatistic[[1]]

  Fst_X1.sum <- summary(lm(Y1 ~ rowSums(X1)))$fstatistic[[1]]
  Fst_X2.sum  <- summary(lm(Y2 ~ rowSums(X2)))$fstatistic[[1]]

  result <- list(data = list(X1 = X1, X2 = X2,  Y1 = Y1, Y2 = Y2, con = as.matrix(con1,n)),
                 Fst = c(Fst_X1 = Fst_X1, Fst_X2 = Fst_X2, Fst_X1.1 = Fst_X1.1, Fst_X2.1 = Fst_X2.1,
                                     Fst_X1.sum = Fst_X1.sum, Fst_X2.sum = Fst_X2.sum))

  return(result)
}

#need to output the variance it explained as well
