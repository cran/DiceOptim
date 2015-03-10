\name{goldsteinprice}

\alias{goldsteinprice}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 2D test function}
\description{
  Goldstein-Price 2-dimensional test function.
}
\usage{goldsteinprice(x)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a 2-dimensional vector specifying the location where the function is to be evaluated.}
}

\value{
  A real number equal to the goldsteinprice function values at \code{x}
}

\details{
The goldsteinprice (standardized version) function is defined over the domain \code{[0,1]^2}. It has 1 global minimizer : 
x* = c(0.5, 0.25), with minimum f(x*) = -3.129172, and 
3 local minima 
x*,2 = c(0.35, 0.4),  
x*,3 = c(0.95, 0.55), 
x*,4 = c(0.8 , 0.7), with respective minima
f(x*,2) = -2.180396,
f(x*,3) = -1.756143,
f(x*,4) = -0.807367.
}

\author{ Tobias Wagner, T.U. Dortmund, Germany
Victor Picheny, CERFACS, Toulouse, France

D. Ginsbourger, University of Bern, Switzerland. }

\keyword{optimize}
\keyword{internal}


\examples{ 
design <- matrix(runif(200), 200, 2)
response <- apply(design, 1, goldsteinprice)
}


