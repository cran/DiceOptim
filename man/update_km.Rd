\name{update_km}
\alias{update_km}
\title{ Update of a Kriging model when adding new observation }
\description{
  Update of a Kriging model when adding new observation, with or without covariance parameter re-estimation.
}
\usage{
update_km(model, NewX, NewY, NewX_AllreadyExist=FALSE,
CovReEstimate=TRUE,new.noise.var=NULL,kmcontrol=NULL)
}
\arguments{
  \item{model}{ a Kriging model of "km" class }
  \item{NewX}{ a matrix containing the new points of experiments }
  \item{NewY}{ a matrix containing the function values on the points NewX }
  \item{NewX_AllreadyExist}{a boolean: do the NewX points allready exist in the design of experiment model@X ? }
  \item{CovReEstimate}{ optional boolean specfiying if the covariance parameters should be re-estimated at every iteration (defualt value = TRUE) }
  \item{new.noise.var}{ a vector containing the noise variance at each new observations }
  \item{kmcontrol}{ kmcontrol (to develop) }
}

\value{
  Updated Kriging model
}

% \details{}

\author{ 
Clement Chevalier, University of Bern (IMSV) 
Victor Picheny, Ecole Centrale Paris
}

%\keyword{optimize}
%\keyword{internal}

%\examples{}


