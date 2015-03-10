\name{rosenbrock4}
\alias{rosenbrock4}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 4D test function}
\description{
  Rosenbrock 4-dimensional test function.
}
\usage{rosenbrock4(x)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a 4-dimensional vector specifying the location where the function is to be evaluated.}
}

\value{
  A real number equal to the rosenbrock4 function values at \code{x}
}

\details{
The rosenbrock4 (standardized version) function is defined over the domain \code{[0,1]^4}. It has 1 global minimizer : 
x* = c(0.4,0.4,0.4,0.4), with minimum f(x*) = -1.019701,
and an additional local minimizer,
x*,2 = c(0.26667,0.4,0.4,0.4), with minimum f(x*,2) = -1.019691.}

\author{ Tobias Wagner, T.U. Dortmund,

Victor Picheny, Ecole Centrale Paris

D. Ginsbourger, University of Bern. }

\keyword{optimize}
\keyword{internal}


\examples{ 
design <- matrix(runif(400), 100, 4)
response <- apply(design, 1, rosenbrock4)
}


