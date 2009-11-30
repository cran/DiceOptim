\name{DiceOptim-package}
\alias{DiceOptim}
\docType{package}

\title{ Kriging-based optimization methods for computer experiments}

\description{
Sequential and parallel Kriging-based optimization methods relying on the 1-point and multipoints expected improvement criteria.
}

\details{
\tabular{ll}{
Package: \tab DiceOptim\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab November 2009\cr
License: \tab GPL-3\cr
}
}

\author{David Ginsbourger (Stochastic Characterization Hydrogeology Group, University of Neuchatel, Switzerland)  

Olivier Roustant (Ecole Nationale Superieure des Mines de Saint-Etienne, France).}

\note{
This work was conducted within the frame of the DICE (Deep Inside Computer Experiments) Consortium between ARMINES, Renault, EDF, IRSN, ONERA and TOTAL S.A. \url{http://www.dice-consortium.fr/}.

The authors would like to thank Yves Deville for his precious advices
in R programming and packaging, as well as the DICE members for useful feedbacks, and especially Yann Richet (IRSN) for numerous discussions concerning the user-friendliness of this package.

Package \code{rgenoud} >=5.3.3. is recommended.

Important functions or methods:
\tabular{ll}{
\code{EI} \tab One-point noise-free EI criterion\cr

\code{qEI} \tab q-points noise-free EI criterion\cr

\code{EGO.nsteps} \tab Sequential EI Algorithm ---model updates including re-estimation of covariance\cr

 \tab hyperparameters--- with a fixed number of iterations (nsteps)\cr

\code{max_EI} \tab One-point noise-free EI maximization. No need to specify any objective function\cr

\code{max_qEI.CL} \tab (sub-)maximization of the q-points EI, based on the Constant Liar heuristic.\cr

   \tab No need to specify any objective function\cr

\code{CL.nsteps} \tab Parallel EI Algorithm ---model updates including re-estimation of covariance\cr

 \tab hyperparameters--- with a fixed number of iterations (nsteps)
}

}



\references{
N.A.C. Cressie (1993), \emph{Statistics for spatial data}, Wiley
series in probability and mathematical statistics.

D. Ginsbourger (2009), \emph{Multiples metamodeles pour l'approximation et l'optimisation de fonctions numeriques multivariables}, Ph.D. thesis, Ecole Nationale Superieure des
Mines de Saint-Etienne, 2009. 
\url{http://members.unine.ch/david.ginsbourger/recherche/these.htm}

D. Ginsbourger, R. Le Riche, and L. Carraro (2009), chapter "Kriging is well-suited to parallelize optimization", to appear in \emph{Computational Intelligence in Expensive Optimization Problems}, Studies in Evolutionary Learning and Optimization, Springer.

D.R. Jones (2001), A taxonomy of global optimization methods based on response surfaces, \emph{Journal of Global Optimization}, 21, 345-383.

D.R. Jones, M. Schonlau, and W.J. Welch (1998), Efficient global optimization of expensive black-box functions, \emph{Journal of Global Optimization}, 13, 455-492.

W.R. Jr. Mebane and J.S. Sekhon (2009), in press, Genetic optimization using derivatives: The rgenoud package for R, \emph{Journal of Statistical Software}.

J. Mockus (1988), \emph{Bayesian Approach to Global Optimization}. Kluwer academic publishers.

C.E. Rasmussen and C.K.I. Williams (2006), \emph{Gaussian Processes for Machine Learning}, the MIT Press, \url{http://www.GaussianProcess.org/gpml}

B.D. Ripley (1987), \emph{Stochastic Simulation}, Wiley.

T.J. Santner, B.J. Williams, and W.J. Notz (2003), \emph{The design and analysis of computer experiments}, Springer.

M. Schonlau (1997), \emph{Computer experiments and global optimization}, Ph.D. thesis, University of Waterloo.

}



\keyword{package}