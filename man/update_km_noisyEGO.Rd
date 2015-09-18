\name{update_km_noisyEGO}
\alias{update_km_noisyEGO}
\title{ Update of one or two Kriging models when adding new observation}
\description{
  Update of a noisy Kriging model when adding new observation, with or without covariance parameter re-estimation.
  When the noise level is unkown, a twin model "estim.model" is also updated.
}
\usage{
update_km_noisyEGO(model, x.new, y.new, noise.var=0, type="UK", 
add.obs=TRUE, index.in.DOE=NULL, CovReEstimate=TRUE, 
NoiseReEstimate=FALSE, estim.model=NULL, nugget.LB=1e-5,
cluster=NULL)
}
\arguments{
  \item{model}{ a Kriging model of "km" class }
  \item{x.new}{ a matrix containing the new points of experiments }
  \item{y.new}{ a matrix containing the function values on the points NewX }
  \item{noise.var}{scalar: noise variance}
  \item{type}{kriging type: "SK" or "UK"}
  \item{add.obs}{boolean: if TRUE, the new point does not exist already in the design of experiment model@X}
  \item{index.in.DOE}{optional integer: if add.obs=TRUE, it specifies the index of the observation in model@X corresponding to x.new}
  \item{CovReEstimate}{ optional boolean specfiying if the covariance parameters should be re-estimated (default value = TRUE) }
  \item{NoiseReEstimate}{ optional boolean specfiying if the noise variance should be re-estimated (default value = TRUE) }
  \item{estim.model}{ optional input of "km" class. Required if NoiseReEstimate=TRUE, in order to deal with repetitions.}
  \item{nugget.LB}{optional scalar: is used to define a lower bound on the noise variance.}
  \item{cluster}{to invoke \code{\link{km}} and \code{\link{genoud}} in parallel, 
    an object of the 'cluster' class returned
    by one of the \code{\link{makeCluster}} commands in the
    \sQuote{parallel} 
    package, or a vector of machine names.  See
    \code{\link{genoud}} for further details.  Currently
    only supported for \dQuote{AKG} criterion.
  }
}

\value{
A list containing:
  \item{model}{ The updated Kriging model }
  \item{estim.model}{If NoiseReEstimate=TRUE, the updated estim.model}
  \item{noise.var}{If NoiseReEstimate=TRUE, the re-estimated noise variance}
}

% \details{}

\author{ 
Victor Picheny (INRA, Castanet-Tolosan, France)
}

%\keyword{optimize}
%\keyword{internal}

%\examples{}


