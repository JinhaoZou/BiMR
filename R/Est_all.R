#' Function to estimation causal effects using one of four methods: Ratio, BiRatio, LIML, BiLIML or all four method
#'
#' @param data data is a list inludes X1, X2, Y1, Y2, confounder, F statistics of X1 on Y1, and F statistics of X2 on Y2
#' The X1 and X2 is matrix with number of individuals * number of IVs, Y1, Y2, confounder are vector with length bumber of individuals
#' @param method one method choose from "Ratio", "BiRatio", "Liml", "BiLiml"
#'
#' @return a matrix with dimension 8*4, each column includes: estimation of g12, sd of g12 estimation, CI of g12, estimation of g21, sd of g21 estimation, CI of g21
#' g12 is the causal effect of Y1 on Y2, g21 is the causal effect of Y2 on Y1; each rows includes: results from four methods "Ratio", "BiRatio", "Liml", "BiLiml"
#' @export
#' @examples {
#' data(example_data)
#' attach(example_data)
#' Est_all(data = example_data, method = "all")
#' }

Est_all <- function(data = data, method = c("Ratio", "BiRatio", "Liml", "BiLiml", "all"),...){

  #Remove the F stat from dataset
  data <- list(X1 = data$X1, X2 = data$X2, Y1 = data$Y1, Y2 = data$Y2, con = data$con)

  if(method == "all"){
    Method <- c("Ratio",  "BiRatio", "Liml", "BiLiml")
    result <- sapply(seq(4), FUN = function(x) Est_oneMethod(data = data, method = Method[x]))
    colnames(result) <- Method
  }else{
    result <-  Est_oneMethod(data = data, method = method)
  }
  return(result)
}
