\name{sphere6}
\alias{sphere6}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 6D sphere function}
\description{
  6D Shifted and rotated weighted sphere test function.
}
\usage{sphere6(x)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a 6-dimensional vector specifying the location where the function is to be evaluated.}
}

\value{
  A real number equal to the sphere6 function values at \code{x}
}

\details{
The 6D Shifted and rotated weighted sphere (standardized version)
function is here defined over the domain \code{[0,1]^6}. It has 1 global minimizer : 
x* = c(1,0.8,0.6,0.4,0.2,0), ), with minimum f(x*) = -1.941389.
It has no further local minima.
}

\author{ Tobias Wagner, T.U. Dortmund,

Victor Picheny, Ecole Centrale Paris

D. Ginsbourger, University of Bern. }

\keyword{optimize}
\keyword{internal}


\examples{ 
design <- matrix(runif(400), 100, 4)
response <- apply(design, 1, sphere6)
}


