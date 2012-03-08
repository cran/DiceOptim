\name{EI}
\alias{EI}


\title{ Analytical expression of the Expected Improvement criterion (noise-free version)}

\description{ Computes the Expected Improvement at current location. }

\usage{
EI(x, model, type="UK", envir=NULL)
}

\arguments{

  \item{x}{ a vector representing the input for which one wishes to calculate EI. }

  \item{model}{ an object of class \code{\link[DiceKriging]{km}}. }

  \item{type}{ "SK" or "UK" (by default), depending whether uncertainty related to trend estimation has to be taken into account. }

  \item{envir}{ an optional environment specifying where to assign intermediate values for future gradient calculations. Default is NULL.}

}

\value{
  The expected improvement, defined as 
  \deqn{EI(x) := E[( min(Y(X)) - Y(x))^{+} | Y(X)=y(X)],} 
  where X is the current design of experiments and Y is the random process assumed to have generated the objective function y.
}

\author{
David Ginsbourger (Departement of Mathematics and Statistics, University of Bern, Switzerland)   

Olivier Roustant (Ecole Nationale Superieure des Mines de Saint-Etienne, France).}


\references{

D.R. Jones, M. Schonlau, and W.J. Welch (1998), Efficient global optimization of expensive black-box functions, \emph{Journal of Global Optimization}, 13, 455-492.

J. Mockus (1988), \emph{Bayesian Approach to Global Optimization}. Kluwer academic publishers.

T.J. Santner, B.J. Williams, and W.J. Notz (2003), \emph{The design and analysis of computer experiments}, Springer.

M. Schonlau (1997), \emph{Computer experiments and global optimization}, Ph.D. thesis, University of Waterloo.

}



\seealso{ \code{\link{max_EI}}, \code{\link{EGO.nsteps}}, \code{\link{qEI}} }

\keyword{models}
%\keyword{optimize}

\examples{
set.seed(123)
##########################################################################
### 	  EI SURFACE ASSOCIATED WITH AN ORDINARY KRIGING MODEL        ####
###    OF THE BRANIN FUNCTION KNOWN AT A 9-POINTS FACTORIAL DESIGN    ####
##########################################################################

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
covtype="gauss", control=list(pop.size=50,trace=FALSE), parinit=c(0.5, 0.5))

# graphics
n.grid <- 30
x.grid <- y.grid <- seq(0,1,length=n.grid)
design.grid <- expand.grid(x.grid, y.grid)
#response.grid <- apply(design.grid, 1, branin)
EI.grid <- apply(design.grid, 1, EI,fitted.model1)
z.grid <- matrix(EI.grid, n.grid, n.grid)
contour(x.grid,y.grid,z.grid,25)
title("Expected Improvement for the Branin function known at 9 points")
points(design.fact[,1], design.fact[,2], pch=17, col="blue")
}