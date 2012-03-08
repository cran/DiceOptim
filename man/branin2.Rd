\name{branin2}

\alias{branin2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 2D test function}
\description{
  Branin 2-dimensional test function (standardized version).
}
\usage{branin2(x)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a 2-dimensional vector specifying the location where the function is to be evaluated.}
}

\value{
  A real number equal to the branin2 function values at \code{x}
}

\details{
The branin2 (standardized version) function is defined over the domain \code{[0,1]^2}. It has 3 global minimizers :  
x*,1 = c(0.1239, 0.8183), 
x*,2 = c(0.5428, 0.1517), 
x*.3 = c(0.9617, 0.1650), with minimum f(x*,i) = -1.047410
}

\author{ Tobias Wagner, T.U. Dortmund,
Victor Picheny, Ecole Centrale Paris
D. Ginsbourger, University of Bern. }

\keyword{optimize}
\keyword{internal}


\examples{ 
design <- matrix(runif(200), 200, 2)
response <- apply(design, 1, branin2)
}


