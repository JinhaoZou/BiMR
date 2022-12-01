#' Function to estimation causal effects using one method from four methods: Ratio, BiRatio, LIML, BiLIML
#'
#' @param data data is a list inludes X1, X2, Y1, Y2, confounder, F statistics of X1 on Y1, and F statistics of X2 on Y2
#' The X1 and X2 is matrix with number of individuals * number of IVs, Y1, Y2, confounder are vector with length bumber of individuals
#' @param method one method choose from "Ratio", "BiRatio", "Liml", "BiLiml"
#'
#' @return a vector with length 8, includes: estimation of g12, sd of g12 estimation, CI of g12, estimation of g21, sd of g21 estimation, CI of g21
#' g12 is the causal effect of Y1 on Y2, g21 is the causal effect of Y2 on Y1
#' @export
#' @examples {
#' data(example_data)
#' attach(example_data)
#' Est_oneMethod(data = example_data, method = "Ratio")
#' }

Est_oneMethod <- function(data = data,method = c("Ratio", "BiRatio", "Liml", "BiLiml"),...){
  tryCatch({
    if(method == "Ratio"){
      #print("ratio")
      return(Est_Ratio(data = data))
    }else if(method == "BiRatio"){
      return(Est_BiRatio(data = data))

    }else if(method == "Liml"){
      return(Est_Liml(data = data))

    }else{
      return(Est_BiLiml(data = data))
    }
  },
  error = function(e)
    return(rep(NA, 8))
  )
}
