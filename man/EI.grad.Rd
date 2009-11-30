\name{EI.grad}

\alias{EI.grad}

\title{ Analytical gradient of the Expected Improvement criterion (noise-free version)}

\description{ Computes the gradient of the Expected Improvement at the current location. Only available for Universal Kriging with constant trend (Ordinary Kriging). }

\usage{
EI.grad(x, model, envir)}

%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{x}{ a vector representing the input for which one wishes to calculate \code{\link{EI}}. }

  \item{model}{ an object of class \code{\link[DiceKriging]{km}}. }

  \item{envir}{ an environment specifying where to get intermediate values calculated in \code{\link{EI}}. }

}

\value{
  The gradient of the expected improvement criterion with respect to x. 
Returns 0 at design points (where the gradient does not exist).
}

\author{David Ginsbourger (Departement of Mathematics and Statistics, University of Berne, Switzerland)  

Olivier Roustant (Ecole Nationale Superieure des Mines de Saint-Etienne, France).}


\references{

D. Ginsbourger (2009), \emph{Multiples metamodeles pour l'approximation et l'optimisation de fonctions numeriques multivariables}, Ph.D. thesis, Ecole Nationale Superieure des
Mines de Saint-Etienne, 2009. 
\url{http://members.unine.ch/david.ginsbourger/recherche/these.htm}

J. Mockus (1988), \emph{Bayesian Approach to Global Optimization}. Kluwer academic publishers.

T.J. Santner, B.J. Williams, and W.J. Notz (2003), \emph{The design and analysis of computer experiments}, Springer.

M. Schonlau (1997), \emph{Computer experiments and global optimization}, Ph.D. thesis, University of Waterloo.

}

\seealso{ \code{\link{EI}} }

\keyword{models}
\keyword{optimize}

\examples{
# a 9-points factorial design, and the corresponding response
d <- 2; n <- 9
design.fact <- expand.grid(seq(0,1,length=3), seq(0,1,length=3)) 
names(design.fact)<-c("x1", "x2")
design.fact <- data.frame(design.fact) 
names(design.fact)<-c("x1", "x2")
response.branin <- apply(design.fact, 1, branin)
response.branin <- data.frame(response.branin) 
names(response.branin) <- "y" 

# model identification
fitted.model1 <- km(~1, design=design.fact, response=response.branin, 
covtype="gauss", control=list(pop.size=50), parinit=c(0.5, 0.5))

# graphics
n.grid <- 50
x.grid <- y.grid <- seq(0,1,length=n.grid)
design.grid <- expand.grid(x.grid, y.grid)
#response.grid <- apply(design.grid, 1, branin)
EI.grid <- apply(design.grid, 1, EI,fitted.model1)
#EI.grid <- apply(design.grid, 1, EI.plot,fitted.model1, gr=TRUE)

z.grid <- matrix(EI.grid, n.grid, n.grid)

contour(x.grid,y.grid,z.grid,20)
title("Expected Improvement for the Branin function known at 9 points")
points(design.fact[,1], design.fact[,2], pch=17, col="blue")

# graphics
n.gridx <- 15
n.gridy <- 20
x.grid2 <- seq(0,1,length=n.gridx) 
y.grid2 <- seq(0,1,length=n.gridy) 
design.grid2 <- expand.grid(x.grid2, y.grid2)

EI.envir <- new.env()	
	environment(EI) <- environment(EI.grad) <- EI.envir 

for(i in seq(1, nrow(design.grid2)) )
{
	x <- design.grid2[i,]
	ei <- EI(x, model=fitted.model1, envir=EI.envir)
	eigrad <- EI.grad(x , model=fitted.model1, envir=EI.envir)
	if(!(is.null(ei)))
	{
	arrows(x$Var1,x$Var2,
	x$Var1 + eigrad[1]*2.2*10e-5, x$Var2 + eigrad[2]*2.2*10e-5, 
	length = 0.04, code=2, col="orange", lwd=2)
	}
}

}
