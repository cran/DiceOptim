\name{hartman4}
\alias{hartman4}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 4D test function}
\description{
  Hartman 4-dimensional test function.
}
\usage{hartman4(x)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a 4-dimensional vector specifying the location where the function is to be evaluated.}
}

\value{
  A real number equal to the hartman4 function values at \code{x}
}

\details{
The hartman4 (standardized version) function is defined over the domain \code{[0,1]^4}. It has 1 global minimizer : 
x* = c(0.1873, 0.1906, 0.5566, 0.2647), with minimum f(x*) = -3.135474
}

\author{ Tobias Wagner, T.U. Dortmund,
Victor Picheny, Ecole Centrale Paris
D. Ginsbourger, University of Bern. }

\keyword{optimize}
\keyword{internal}


\examples{ 
design <- matrix(runif(400), 100, 4)
response <- apply(design, 1, hartman4)
}


