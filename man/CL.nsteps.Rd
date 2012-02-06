\name{CL.nsteps}
\alias{CL.nsteps}

\title{ Parallelized version of EGO.nsteps, based on the CL strategy }

\description{ Adaptation of the EI Algorithm for synchronous parallel computing, based on the "Constant Liar" heuristic strategy.
 Indeed, \code{\link{max_qEI.CL}} is used at every iteration to provide a new q-points design. 
The Kriging model is updated after each batch of objective function evaluations, including re-estimation of the covariance hyperparameters.
 
}

\usage{
CL.nsteps(model, fun, npoints, nsteps, lower, upper, parinit = NULL, 
kmcontrol = NULL, control = NULL)
}


\arguments{

  \item{model}{  an object of class \code{\link[DiceKriging]{km}} , }

  \item{fun}{ the objective function to be minimized, }

  \item{npoints}{ an integer representing the desired number of parallel computations, }

  \item{nsteps}{ an integer representing the desired number of batches, }

  \item{lower}{ vector of lower bounds for the variables to be optimized over, }

  \item{upper}{ vector of upper bounds for the variables to be optimized over, }

  \item{parinit}{  optional vector of initial values for the variables to be optimized over, }

  \item{control}{  an optional list of control parameters for optimization. One can control  

		 \code{"pop.size"} (default : [4+3*log(nb of variables)]), 

	       \code{"max.generations"}  (default :5), 

		 \code{"wait.generations"} (default :2), 

		\code{"BFGSburnin"} (default :0),

		of the function \code{\link[rgenoud]{genoud}}. }

  \item{kmcontrol}{ an optional list representing the control variables for the re-estimation of the kriging model. 
		    The items are the same as in \code{\link[DiceKriging]{km}} : 

		\code{penalty}, \code{optim.method}, \code{parinit}, \code{control}. 

		The default values are those contained in \code{model}, typically corresponding to the variables used in \code{\link[DiceKriging]{km}} to estimate a kriging model from the initial design points.}
}


\value{ A list with components:

\item{par}{a data frame representing the additional points visited during the algorithm,}

\item{value}{a data frame representing the response values at the points given in \code{par},}

\item{npoints}{an integer representing the number of parallel computations (given in argument),}

\item{nsteps}{an integer representing the desired number of iterations (given in argument),}

\item{lastmodel}{an object of class \code{\link[DiceKriging]{km}} corresponding to the last kriging model fitted.}

}




\author{David Ginsbourger (Departement of Mathematics and Statistics, University of Berne, Switzerland)  

Olivier Roustant (Ecole Nationale Superieure des Mines de Saint-Etienne, France).}


\seealso{ \code{\link{max_qEI.CL}}, \code{\link{qEI}}, \code{\link{EI}}, \code{\link{max_EI}}, \code{\link{EGO.nsteps}} }

%\keyword{models}
\keyword{optimize}

\references{

D. Ginsbourger (2009), \emph{Multiples metamodeles pour l'approximation et l'optimisation de fonctions numeriques multivariables}, Ph.D. thesis, Ecole Nationale Superieure des
Mines de Saint-Etienne, 2009. 
\url{http://members.unine.ch/david.ginsbourger/recherche/these.htm}

D. Ginsbourger, R. Le Riche, and L. Carraro (2009), chapter "Kriging is well-suited to parallelize optimization", to appear in \emph{Computational Intelligence in Expensive Optimization Problems}, Studies in Evolutionary Learning and Optimization, Springer.

W.R. Jr. Mebane and J.S. Sekhon (2009), in press, Genetic optimization using derivatives: The rgenoud package for R, \emph{Journal of Statistical Software}.

B.D. Ripley (1987), \emph{Stochastic Simulation}, Wiley.

T.J. Santner, B.J. Williams, and W.J. Notz (2003), \emph{The design and analysis of computer experiments}, Springer.

M. Schonlau (1997), \emph{Computer experiments and global optimization}, Ph.D. thesis, University of Waterloo.

}

\examples{

set.seed(123)
# 3 Iterations of the Constant Liar are applied to Branin's function, 
# with 8 points in parallel at each iteration 

# a 9-points factorial design, and the corresponding response
d <- 2
n <- 9
design.fact <- expand.grid(seq(0,1,length=3), seq(0,1,length=3)) 
design.fact <- data.frame(design.fact)
names(design.fact)<-c("x1", "x2")
response.branin <- apply(design.fact, 1, branin)
response.branin <- data.frame(response.branin)
names(response.branin) <- "y" 

# model identification
fitted.model1 <- km(~1, design=design.fact, response=response.branin, 
covtype="gauss", control=list(pop.size=50,trace=FALSE), parinit=c(0.5, 0.5))

# parallel EGO 
library(rgenoud)
nsteps <- 3
npoints <- 8
lower <- rep(0,d) 
upper <- rep(1,d) 

oEGOparallel1 <- CL.nsteps(model=fitted.model1, fun=branin, npoints=npoints, 
nsteps=nsteps, lower, upper, control=list(pop.size=20, BFGSburnin=2)) 

print(oEGOparallel1)

# graphics
n.grid <- 20
x.grid <- y.grid <- seq(0,1,length=n.grid)
design.grid <- expand.grid(x.grid, y.grid)
response.grid <- apply(design.grid, 1, branin)
z.grid <- matrix(response.grid, n.grid, n.grid)
contour(x.grid,y.grid,z.grid,40)
title("3 iterations of Constant Liar with 8 parallel searches")
points(design.fact[,1], design.fact[,2], pch=17, col="blue")
points(oEGOparallel1$par[1:npoints,], pch=19, col="red")
points(oEGOparallel1$par[(npoints+1):(2*npoints),], pch=19, col="violet")
points(oEGOparallel1$par[(2*npoints+1):(3*npoints),], pch=19, col="green")
text(oEGOparallel1$par[,1], oEGOparallel1$par[,2], 
labels=1:(npoints*nsteps), pos=3)

}
