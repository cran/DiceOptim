##' Check that the new point is not too close to already known observations to avoid numerical issues.
##' Closeness can be estimated with several distances.
##' @title Prevention of numerical instability for a new observation
##' @param x a vector representing the input to check,
##' @param model list of objects of class \code{\link[DiceKriging]{km}}, one for each objective functions,
##' @param threshold optional value for the minimal distance to an existing observation, default to \code{1e-4},
##' @param distance selection of the distance between new observations, between "\code{euclidean}",
##'  "\code{covdist}" (default) and "\code{covratio}", see details,
##' @param type "\code{SK}" or "\code{UK}" (default), depending whether uncertainty related to trend estimation has to be taken into account.
##' @details If the distance between \code{x} and the closest observations in \code{model} is below
##'  \code{threshold}, \code{x} should not be evaluated to avoid numerical instabilities.
##' The distance can simply be the Euclidean distance or the canonical distance associated with the kriging covariance k: 
##'   \deqn{d(x,y) = \sqrt{k(x,x) - 2k(x,y) + k(y,y)}.}{d(x,y) = \sqrt(k(x,x) - 2k(x,y) + k(y,y)).} 
##'   The last solution is the ratio between the prediction variance at \code{x} and the variance of the process.
##' @return \code{TRUE} if the point should not be tested.
##' @author Mickael Binois 
##' @export
checkPredict <- function(x, model, threshold = 1e-4, distance = "covdist", type = "UK"){
  
  if(is.null(dim(x))){
    x <- matrix(x, nrow = 1) 
  }
  
  if(is.null(distance))
    distance <- "covdist"
  if(is.null(threshold))
    threshold <- 1e-4

  if(distance == "none") return(rep(FALSE, nrow(x)))
    
  if(distance == "euclidean"){
    tp1 <- as.numeric(t(model[[1]]@X))
    tp2 <- matrix(tp1 - as.numeric(x), nrow = model[[1]]@n, byrow = TRUE)^2
    mindist <- min(sqrt(rowSums(tp2)))
  }else{
    # d <- length(which(sapply(model, class) == "km")) # to filter the fastobj (no troubles to update)
    if(distance == "covratio"){
      
      mindist <- Inf
      for(i in 1:length(model)){
        if(class(model[[i]]) != "fastfun"){
          pred.sd <- predict(object = model[[i]], newdata = x, type = type, checkNames = FALSE, light.return = TRUE)
          mindist <- pmin(mindist, pred.sd$sd/sqrt(model[[i]]@covariance@sd2))
        }
      }
    }else{
      
      # if one model is likely to be unable to take the new point without numerical instabilities,
      # it should not be evaluated
      mindist <- Inf
      for(i in 1:length(model)){
        if(class(model[[i]]) != "fastfun"){
          kn_xx <- (predict(model[[i]], newdata = x, type = type, checkNames = FALSE)$sd)^2 # k_n(x,x)
          kn_yy <- (predict(model[[i]], newdata = model[[i]]@X, type = type, checkNames = FALSE)$sd)^2 # k_n(y,y)
          
          kxy <- covMat1Mat2(model[[i]]@covariance, x, model[[i]]@X)
          kyy <- covMatrix(model[[i]]@covariance, model[[i]]@X)$C
          kn_xy <- kxy - kxy %*% chol2inv(model[[i]]@T) %*% kyy # k_n(x,y)

          mindist <- pmin(mindist,
                          sqrt(pmax(0, kn_xx - 2 * apply(kn_xy - matrix(kn_yy, nrow(x), model[[i]]@n, byrow = T), 1, max))/model[[i]]@covariance@sd2))
        }
      } 
    }
  }
  
  return(mindist < threshold)
}
