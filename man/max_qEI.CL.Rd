\name{max_qEI.CL}
\alias{max_qEI.CL}

\title{ One-shot pseudo-maximization of qEI using the Constant Liar strategy }

\description{
Approached maximization of the q-points Expected Improvement criterion \code{\link{qEI}} using the "constant liar" heuristic. 
First, the regular Expected Improvement \code{\link{EI}} is maximized. Then, for the next points, 
the Expected Improvement is maximized again, but with an artificially updated Kriging model. 
Since the response values corresponding to the last best point obtained are not availbale, 
the idea of CL is to replace them by an arbitrary constant value L (a "lie") set by the user.
}

\usage{
max_qEI.CL(model, npoints, L, lower, upper, parinit=NULL, control=NULL)
}

\arguments{

  \item{model}{  an object of class \code{\link[DiceKriging]{km}} , }

  \item{npoints}{ an integer representing the desired number of iterations, }

  \item{L}{ a real number equal to the "constant liar" value. Typically, L is fixed to the smaller observed function value.}

  \item{lower}{ vector of lower bounds for the variables to be optimized over, }

  \item{upper}{ vector of upper bounds for the variables to be optimized over, }

  \item{parinit}{  optional vector of initial values for the variables to be optimized over, }

  \item{control}{  an optional list of control parameters for optimization. One can control  

		 \code{"pop.size"} (default : [4+3*log(nb of variables)]), 

	       \code{"max.generations"}  (default :5), 

		 \code{"wait.generations"} (default :2), 

		\code{"BFGSburnin"} (default :0),

		of the function \code{\link[rgenoud]{genoud}}. }

}




\value{ A list with components:

\item{par}{A matrix containing the \code{npoints} input vectors found.}

\item{value}{A vector giving the \code{npoints} values of artificial responses (all equal to L here).}

}


\author{David Ginsbourger (Departement of Mathematics and Statistics, University of Berne, Switzerland)  

Olivier Roustant (Ecole Nationale Superieure des Mines de Saint-Etienne, France).}


\references{

D. Ginsbourger (2009), \emph{Multiples metamodeles pour l'approximation et l'optimisation de fonctions numeriques multivariables}, Ph.D. thesis, Ecole Nationale Superieure des
Mines de Saint-Etienne, 2009. 
\url{http://www.ginsbourger.ch/recherche/these.php}

D. Ginsbourger, R. Le Riche, and L. Carraro (2009), chapter "Kriging is well-suited to parallelize optimization", to appear in \emph{Computational Intelligence in Expensive Optimization Problems}, Studies in Evolutionary Learning and Optimization, Springer.

W.R. Jr. Mebane and J.S. Sekhon (2009), in press, Genetic optimization using derivatives: The rgenoud package for R, \emph{Journal of Statistical Software}.

T.J. Santner, B.J. Williams, and W.J. Notz (2003), \emph{The design and analysis of computer experiments}, Springer.

M. Schonlau (1997), \emph{Computer experiments and global optimization}, Ph.D. thesis, University of Waterloo.

}




\seealso{ \code{\link{CL.nsteps}}, \code{\link{EI}},  \code{\link{max_EI}} }

%\keyword{models}
\keyword{optimize}

\examples{
set.seed(123)
# a 9-points factorial design, and the corresponding response
d <- 2 
n <- 9
design.fact <- expand.grid(seq(0,1,length=3), seq(0,1,length=3))
names(design.fact)<-c("x1", "x2")
design.fact <- data.frame(design.fact) 
names(design.fact)<-c("x1", "x2")
response.branin <- apply(design.fact, 1, branin)
response.branin <- data.frame(response.branin) 
names(response.branin) <- "y" 

# model identification
fitted.model1 <- km(~1, design=design.fact, response=response.branin, 
covtype="gauss", control=list(pop.size=50, trace=FALSE), parinit=c(0.5, 0.5))

# EGO parallel
library(rgenoud)
npoints <- 3
lower <- rep(0,d) 
upper <- rep(1,d) 
oEGOparallel <- max_qEI.CL(model=fitted.model1, npoints=npoints, 
L=min(response.branin), lower, upper, 
control=list(pop.size=20, BFGSburnin=2)) 
print(oEGOparallel)

\dontrun{
# graphics
n.grid <- 20
x.grid <- y.grid <- seq(0,1,length=n.grid)
design.grid <- expand.grid(x.grid, y.grid)
response.grid <- apply(design.grid, 1, branin)
z.grid <- matrix(response.grid, n.grid, n.grid)
contour(x.grid,y.grid,z.grid,40)
title("Branin function")
points(design.fact[,1], design.fact[,2], pch=17, col="blue")
points(oEGOparallel$par, pch=19, col="red")
text(oEGOparallel$par[,1], oEGOparallel$par[,2], labels=1:npoints, pos=3)
}

}
