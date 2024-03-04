##' 4D test function
##' 
##' Hartman 4-dimensional test function.
##' 
##' The hartman4 (standardized version) function is defined over the domain
##' \code{[0,1]^4}. It has 1 global minimizer : x* = c(0.1873, 0.1906, 0.5566,
##' 0.2647), with minimum f(x*) = -3.135474
##' 
##' @param x a 4-dimensional vector specifying the location where the function
##' is to be evaluated.
##' @return A real number equal to the hartman4 function values at \code{x}
##' @author Tobias Wagner  
##' 
##' Victor Picheny 
##' 
##' David Ginsbourger 
##' @keywords optimize internal
##' @examples
##'  
##' design <- matrix(runif(400), 100, 4)
##' response <- apply(design, 1, hartman4)
##' 
##' @export hartman4
hartman4 <- function(x)
{
# 4D-Hartman, q=4 test function (standardized version)
# ----------------------------------------------------
# Dimension: n = 4
# Number of local minima: 1
# The global optimum
# x* = c(0.1873, 0.1906, 0.5566, 0.2647), f(x*) = -3.135474

a <- rbind(c(10.00,  0.05,  3.00, 17.00),
          c(3.00, 10.00,  3.50,  8.00),
          c(17.00, 17.00,  1.70,  0.05),
          c(3.50,  0.10, 10.00, 10.00),
          c(1.70,  8.00, 17.00,  0.10),
          c(8.00, 14.00,  8.00, 14.00))

p <- rbind(c(0.1312, 0.2329, 0.2348, 0.4047),
           c(0.1696, 0.4135, 0.1451, 0.8828),
           c(0.5569, 0.8307, 0.3522, 0.8732),
           c(0.0124, 0.3736, 0.2883, 0.5743),
           c(0.8283, 0.1004, 0.3047, 0.1091),
           c(0.5886, 0.9991, 0.6650, 0.0381))

C <- c(1.0, 1.2, 3.0, 3.2)
m <- -1.1
s <- 0.8387

d <- rep(0,1,4)

for (i in 1:4)
    d[i] = sum(a[seq(1,4),i]*(x - p[seq(1,4),i])^2);
end

f <- (-sum(C*exp(-d)) - m)/s

return(f)
}
