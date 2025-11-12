#' 6D sphere function
#' 
#' 6D Shifted and rotated weighted sphere test function.
#' 
#' The 6D Shifted and rotated weighted sphere (standardized version) function
#' is here defined over the domain \code{[0,1]^6}. It has 1 global minimizer :
#' x* = c(1,0.8,0.6,0.4,0.2,0), ), with minimum f(x*) = -1.941389. It has no
#' further local minima.
#' 
#' @param x a 6-dimensional vector specifying the location where the function
#' is to be evaluated.
#' @return A real number equal to the sphere6 function values at \code{x}
#' @author Tobias Wagner  
#' 
#' Victor Picheny 
#' 
#' David Ginsbourger 
#' @keywords optimize internal
#' @examples
#'  
#' design <- matrix(runif(400), 100, 4)
#' response <- apply(design, 1, sphere6)
#' 
#' @export sphere6
sphere6 <- function(x)
{
# 6D Shifted and rotated weighted sphere (standardized version)
# -------------------------------------------------------------
# Dimension: n = 6
# Number of local minima: 1
# The global minimum:
# x* = c(1,0.8,0.6,0.4,0.2,0), f(x*) = -1.941389

x <- x * 10 - 5
m <- 1745.3796
s <- 899.0367

# shift
shift <- seq(-5, 5, 2)
x <- x + shift

# rotate
rotationAngle <- 22.5
radial <- rotationAngle * pi / 180
Rot <- rbind(c(cos(radial), -sin(radial)), c(sin(radial),cos(radial)))
for (i in seq(1,5))
{
	for (j in seq(i+1,6))
	{
		x[c(i, j)] <- as.double(x[c(i, j)]) %*% Rot
	}
}

# evaluate
f <- sum(x^2*(2^seq(1,6)))
f <- (f-m)/s
return(f)
}
