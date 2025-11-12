#' Strongly multimdoal constraint function from Parr et al. (standardized version)
#' @title 2D constraint function
#' @param x a 2-dimensional vector or a two-column matrix specifying the location(s) where the function
#' is to be evaluated.
#' @return A scalar
#' @export
#' @examples
#' n.grid <- 20
#' x.grid <- y.grid <- seq(0,1,length=n.grid)
#' design.grid <- expand.grid(x.grid, y.grid)
#' response.grid <- apply(design.grid, 1, ParrConstraint)
#' z.grid <- matrix(response.grid, n.grid, n.grid)
#' contour(x.grid,y.grid,z.grid,40)
#' title("Parr constraint function")
#' 
ParrConstraint <- function(x){
  # Constraint function from Parr et al. (standardized version)
  #--------------------------------------------
  # Dimension: n = 2
  
  if(is.null(dim(x))){
    x <- matrix(x, 1) 
  }
  
  x1 <- 2 * x[,1] - 1
  x2 <- 2 * x[,2] - 1
  
  g <- (4-2.1*x1^2+1/3*x1^4)*x1^2 + x1*x2 + (-4+4*x2^2)*x2^2+3*sin(6*(1-x1)) + 3*sin(6*(1-x2 ))
  
  return(-g)
}