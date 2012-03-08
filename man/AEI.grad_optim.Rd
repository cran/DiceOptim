\name{AEI.grad_optim}
\Rdversion{1}
\alias{AEI.grad_optim}

\title{ AEI's Gradient }

\description{
Gradient of the Augmented Expected Improvement (AEI) criterion.}

\usage{
AEI.grad_optim(x, model, new.noise.var=0, y.min, type = "UK", envir)
}

\arguments{
  \item{x}{ the input vector at which one wants to evaluate the criterion }
  \item{model}{ a Kriging model of "km" class }
  \item{new.noise.var}{ the (scalar) noise variance of the new observation.}
  \item{y.min}{ The kriging predictor at the current best point (point with smallest kriging quantile). If not provided, this quantity is evaluated. } 
  \item{type}{ Kriging type: "SK" or "UK" }
  \item{envir}{ environment for inheriting intermediate calculations from AEI }
}
%\details{}
\value{Gradient of the Augmented Expected Improvement}
%\references{}
\author{ 
Victor Picheny (CERFACS, Toulouse, France)
David Ginsbourger (IMSV, University of Bern)
}
%\note{}

%\seealso{}

%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory. \keyword{}
%\keyword{internal}
