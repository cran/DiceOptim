\name{EQI.grad_optim}
\Rdversion{1}
\alias{EQI.grad_optim}

\title{ EQI's Gradient, for usage in max_EQI }

\description{
Gradient of the Expected Quantile Improvement (EQI) criterion.
}

\usage{
EQI.grad_optim(x, model, new.noise.var=0, alpha=0.9, q.min=NULL, type = "UK", envir)
}

\arguments{
  \item{x}{ the input vector at which one wants to evaluate the criterion }
  \item{model}{ a Kriging model of "km" class }
  \item{new.noise.var}{ (scalar) noise variance of the future observation. 
Default value is 0 (noise-free observation).}  
  \item{alpha}{ Quantile level (default value is 0.9) }  
  \item{q.min}{ Best kriging quantile. If not provided, this quantity is evaluated. } 
  \item{type}{ Kriging type: "SK" or "UK" }
  \item{envir}{ environment for inheriting intermediate calculations from EQI }
}
%\details{}
\value{Gradient of the Expected Quantile Improvement}
%\references{}
\author{ 
Victor Picheny (CERFACS, Toulouse, France) 
David Ginsbourger (IMSV, University of Bern, Switzerland)
}
%\note{}

%\seealso{}

%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory. \keyword{}
%\keyword{}
