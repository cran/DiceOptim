#' 4D test function
#' 
#' Rosenbrock 4-dimensional test function.
#' 
#' The rosenbrock4 (standardized version) function is defined over the domain
#' \code{[0,1]^4}. It has 1 global minimizer : x* = c(0.4,0.4,0.4,0.4), with
#' minimum f(x*) = -1.019701, and an additional local minimizer, x*,2 =
#' c(0.26667,0.4,0.4,0.4), with minimum f(x*,2) = -1.019691.
#' 
#' @param x a 4-dimensional vector specifying the location where the function
#' is to be evaluated.
#' @return A real number equal to the rosenbrock4 function values at \code{x}
#' @author Tobias Wagner  
#' 
#' Victor Picheny 
#' 
#' David Ginsbourger 
#' @keywords optimize internal
#' @examples
#'  
#' design <- matrix(runif(400), 100, 4)
#' response <- apply(design, 1, rosenbrock4)
#' 
#' @export rosenbrock4
rosenbrock4 <- function(x)
{
# 4D-Rosenbrock test function (standardized version)
#---------------------------------------------------
# Dimension: n = 4
# Number of local minima: 2
# The global minimum: 
# x* = c(0.4,0.4,0.4,0.4), f(x*) = -1.019701
# The local minimum:
# x*,2 = c(0.26667,0.4,0.4,0.4), f(x*,2) = -1.019691

m <- 382658.057227524
s <- 375264.858362295

x <- 15 * x - 5

x1 <- x[seq(1,3)]
x2 <- x[seq(2,4)]
f <- sum(100*(x2 - x1^2)^2 + (1 - x1)^2);
f <- (f-m)/s
return(f)
}
