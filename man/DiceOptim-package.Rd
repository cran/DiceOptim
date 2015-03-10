\name{DiceOptim-package}
\alias{DiceOptim}
\docType{package}

\title{ Kriging-based optimization methods for computer experiments}

\description{
Sequential and parallel Kriging-based optimization methods relying on expected improvement criteria.
}

\details{
\tabular{ll}{
Package: \tab DiceOptim\cr
Type: \tab Package\cr
Version: \tab 1.5\cr
Date: \tab March 2015\cr
License: \tab GPL-2 | GPL-3\cr
}
}

\author{
David Ginsbourger (Departement of Mathematics and Statistics, University of Bern, Switzerland)  

Victor Picheny (INRA, Castanet-Tolosan, France)

Olivier Roustant (Ecole Nationale Superieure des Mines de Saint-Etienne, France).

with contributions by C. Chevalier, S. Marmin, and T. Wagner 
}

\note{
This work is a follow-up of DiceOptim 1.0, which was conducted within the frame of the DICE (Deep Inside Computer Experiments) Consortium 
between ARMINES, Renault, EDF, IRSN, ONERA and TOTAL S.A. (http://www.dice-consortium.fr/).

The authors would like to thank Yves Deville for his precious advice in R programming and packaging, as well as DICE and ReDICE members for useful feedbacks, and especially Yann Richet (IRSN) for numerous discussions concerning the user-friendliness of this package.

Package \code{rgenoud} >=5.3.3. is recommended.

Important functions or methods:
\tabular{ll}{
\code{EI} \tab One-point Expected Improvement criterion\cr

\code{qEI} \tab Analytical expression of the multipoint Expected Improvement (qEI) criterion \cr

\code{EGO.nsteps} \tab Sequential EI Algorithm ---model updates including re-estimation of covariance\cr

 \tab hyperparameters--- with a fixed number of iterations (nsteps)\cr

\code{max_EI} \tab One-point noise-free EI maximization. No need to specify any objective function\cr

\code{max_qEI} \tab Maximization of multipoint expected improvement criterion (qEI) \cr

\code{qEGO.nsteps} \tab Batch-sequential EI Algorithm ---model updates including re-estimation of covariance\cr

\code{EQI} \tab One-point noisy EI-like criterion (Expected Quantile Improvement)\cr

\code{AEI} \tab One-point noisy EI-like criterion (Augmented Expected Improvement)\cr

\code{AKG} \tab One-point noisy EI criterion (Approximate Knowledge Gradient)\cr

\code{noisy.optimizer} \tab Sequential EI-like Algorithms for Noisy Kriging-based Optimization
}
}



\references{
C. Chevalier and D. Ginsbourger (2014) Learning and Intelligent Optimization - 7th International Conference, Lion 7, Catania, Italy, January 7-11, 2013, 
Revised Selected Papers, chapter Fast computation of the multipoint Expected Improvement with applications in batch selection, pages 59-69, Springer.

N.A.C. Cressie (1993), \emph{Statistics for spatial data}, Wiley
series in probability and mathematical statistics.

D. Ginsbourger, R. Le Riche, L. Carraro (2007), A Multipoint Criterion for Deterministic Parallel Global Optimization based on Kriging. 
The International Conference on Non Convex Programming, 2007.

D. Ginsbourger (2009), \emph{Multiples metamodeles pour l'approximation et l'optimisation de fonctions numeriques multivariables}, Ph.D. thesis, Ecole Nationale Superieure des
Mines de Saint-Etienne, 2009. 
\url{http://www.ginsbourger.ch/recherche/these.php}

D. Ginsbourger, R. Le Riche, and L. Carraro (2010), chapter "Kriging is well-suited to parallelize optimization", in \emph{Computational Intelligence in Expensive Optimization Problems}, Studies in Evolutionary Learning and Optimization, Springer.

D.R. Jones (2001), A taxonomy of global optimization methods based on response surfaces, \emph{Journal of Global Optimization}, 21, 345-383.

D.R. Jones, M. Schonlau, and W.J. Welch (1998), Efficient global optimization of expensive black-box functions, \emph{Journal of Global Optimization}, 13, 455-492.

W.R. Jr. Mebane and J.S. Sekhon (2009), in press, Genetic optimization using derivatives: The rgenoud package for R, \emph{Journal of Statistical Software}.

J. Mockus (1988), \emph{Bayesian Approach to Global Optimization}. Kluwer academic publishers.

V. Picheny and D. Ginsbourger (2013), Noisy kriging-based optimization methods: A unified implementation within the DiceOptim package, 
\emph{Computational Statistics & Data Analysis}, \url{http://dx.doi.org/10.1016/j.csda.2013.03.018}

V. Picheny, D. Ginsbourger, Y. Richet and G. Caplin (2013), Quantile-Based Optimization of Noisy Computer Experiments With Tunable Precision (with discussion and rejoinder),
\emph{Technometrics}, \url{http://www.tandfonline.com/doi/abs/10.1080/00401706.2012.707580}

C.E. Rasmussen and C.K.I. Williams (2006), \emph{Gaussian Processes for Machine Learning}, the MIT Press, \url{http://www.GaussianProcess.org/gpml}

B.D. Ripley (1987), \emph{Stochastic Simulation}, Wiley.

O. Roustant, D. Ginsbourger and Yves Deville (2012), DiceKriging, DiceOptim: Two R Packages for the Analysis of Computer Experiments by Kriging-Based Metamodeling and Optimization, \emph{Journal of Statistical Software}, 
\bold{51}(1), 1-55, \url{http://www.jstatsoft.org/v51/i01/}.

T.J. Santner, B.J. Williams, and W.J. Notz (2003), \emph{The design and analysis of computer experiments}, Springer.

M. Schonlau (1997), \emph{Computer experiments and global optimization}, Ph.D. thesis, University of Waterloo.
}



\keyword{package}