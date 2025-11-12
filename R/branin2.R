#' 2D test function
#' 
#' Branin 2-dimensional test function (standardized version).
#' 
#' The branin2 (standardized version) function is defined over the domain
#' \code{[0,1]^2}. It has 3 global minimizers : x*,1 = c(0.1239, 0.8183), x*,2
#' = c(0.5428, 0.1517), x*.3 = c(0.9617, 0.1650), with minimum f(x*,i) =
#' -1.047410
#' 
#' @param x a 2-dimensional vector specifying the location where the function
#' is to be evaluated.
#' @return A real number equal to the branin2 function values at \code{x}
#' @author Tobias Wagner 
#' 
#' Victor Picheny 
#' 
#' David Ginsbourger 
#' @keywords optimize internal
#' @examples
#'  
#' design <- matrix(runif(200), 200, 2)
#' response <- apply(design, 1, branin2)
#' 
#' @export branin2
branin2 <- function(x)
{
# Branin test function (standardized version)
#--------------------------------------------
# Dimension: n = 2
# Number of local minima: 3 (the global ones)
# The global minima: 
# x*,1 = c(0.1239, 0.8183), f(x*,1) = -1.047410
# x*,2 = c(0.5428, 0.1517), f(x*,2) = -1.047410
# x*.3 = c(0.9617, 0.1650), f(x*,3) = -1.047410

m <- 54.8104
s <- 51.9496

xx <- 15 * x[1] - 5
y <- 15 * x[2]

f <- (y - 5.1*xx^2/(4*pi^2) + 5*xx/pi - 6)^2 + 10*(1 - 1/(8*pi))*cos(xx) + 10
f <- (f-m)/s

return(f)
}
