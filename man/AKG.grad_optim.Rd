\name{AKG.grad_optim}
\Rdversion{1}
\alias{AKG.grad_optim}

\title{ AKG's Gradient, for usage in max_AKG }

\description{
Gradient of the Approximate Knowledge Gradient (AKG) criterion.
}

\usage{
AKG.grad_optim(x, model, new.noise.var=0, type = "UK", envir)
}

\arguments{
  \item{x}{ the input vector at which one wants to evaluate the criterion }
  \item{model}{ a Kriging model of "km" class }
  \item{new.noise.var}{ (scalar) noise variance of the future observation. Default value is 0 (noise-free observation).}  
  \item{type}{ Kriging type: "SK" or "UK" }
  \item{envir}{ environment for inheriting intermediate calculations from AKG }
}
%\details{}
\value{Gradient of the Approximate Knowledge Gradient}
%\references{}
\author{ 
Victor Picheny
}
%\note{}

%\seealso{}

%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory. \keyword{}
%\keyword{}
