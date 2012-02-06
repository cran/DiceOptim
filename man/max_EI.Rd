\name{max_EI}
\alias{max_EI}

\title{ One-shot maximization of the Expected Improvement criterion (noise-free version)}

\description{  
Given an object of class \code{\link[DiceKriging]{km}} and a set of tuning parameters (\code{lower},\code{upper},\code{parinit}, and \code{control}), 
\code{max_EI} performs the maximization of the Expected Improvement criterion
and delivers the next point to be visited in an EGO-like procedure. 

The latter maximization relies on a genetic algorithm using derivatives, \code{\link[rgenoud]{genoud}}.
This function plays a central role in the package since it is in constant use in the proposed algorithms.
It is important to remark that the information needed about the objective function reduces here
to the vector of response values embedded in \code{model} (no call to the objective function or simulator). 
}

\usage{
max_EI(model, lower, upper, parinit=NULL, control = NULL) }

\arguments{

  \item{model}{  an object of class \code{\link[DiceKriging]{km}} , }


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

\item{par}{The best set of parameters found.}

\item{value}{The value of expected improvement at par.}

}

\note{ Most EGO-like methods (EI algorithms) usually work with Ordinary Kriging (constant trend), by maximization of the expected improvement. 
Here, the EI maximization is also possible with any linear trend. However, note that the optimization may perform much faster and better
when the trend is a constant since it is the only case where the analytical gradient is available. }


\author{David Ginsbourger (Departement of Mathematics and Statistics, University of Berne, Switzerland)   

Olivier Roustant (Ecole Nationale Superieure des Mines de Saint-Etienne, France).}


\references{

D. Ginsbourger (2009), \emph{Multiples metamodeles pour l'approximation et l'optimisation de fonctions numeriques multivariables}, Ph.D. thesis, Ecole Nationale Superieure des
Mines de Saint-Etienne, 2009. 
\url{http://members.unine.ch/david.ginsbourger/recherche/these.htm}

D.R. Jones, M. Schonlau, and W.J. Welch (1998), Efficient global optimization of expensive black-box functions, \emph{Journal of Global Optimization}, 13, 455-492.

W.R. Jr. Mebane and J.S. Sekhon (2009), in press, Genetic optimization using derivatives: The rgenoud package for R, \emph{Journal of Statistical Software}.

}

%\keyword{models}
\keyword{optimize}

\examples{
set.seed(123)
##########################################################
### "ONE-SHOT" EI-MAXIMIZATION OF THE BRANIN FUNCTION ####
### 	KNOWN AT A 9-POINTS FACTORIAL DESIGN          ####
##########################################################

# a 9-points factorial design, and the corresponding response
d <- 2 
n <- 9
design.fact <- expand.grid(seq(0,1,length=3), seq(0,1,length=3)) 
names(design.fact) <- c("x1", "x2")
design.fact <- data.frame(design.fact) 
names(design.fact) <- c("x1", "x2")
response.branin <- apply(design.fact, 1, branin)
response.branin <- data.frame(response.branin) 
names(response.branin) <- "y" 

# model identification
fitted.model1 <- km(~1, design=design.fact, response=response.branin, 
covtype="gauss", control=list(pop.size=50,trace=FALSE), parinit=c(0.5, 0.5))

# EGO one step
library(rgenoud)
lower <- rep(0,d) 
upper <- rep(1,d)     # domain for Branin function
oEGO <- max_EI(fitted.model1, lower, upper, 
control=list(pop.size=20, BFGSburnin=2))
print(oEGO)

# graphics
n.grid <- 20
x.grid <- y.grid <- seq(0,1,length=n.grid)
design.grid <- expand.grid(x.grid, y.grid)
response.grid <- apply(design.grid, 1, branin)
z.grid <- matrix(response.grid, n.grid, n.grid)
contour(x.grid,y.grid,z.grid,40)
title("Branin Function")
points(design.fact[,1], design.fact[,2], pch=17, col="blue")
points(oEGO$par[1], oEGO$par[2], pch=19, col="red")


#############################################################
### "ONE-SHOT" EI-MAXIMIZATION OF THE CAMELBACK FUNCTION ####
###	KNOWN AT A 16-POINTS FACTORIAL DESIGN            ####
#############################################################

# a 16-points factorial design, and the corresponding response
d <- 2 
n <- 16
design.fact <- expand.grid(seq(0,1,length=4), seq(0,1,length=4)) 
names(design.fact)<-c("x1", "x2")
design.fact <- data.frame(design.fact) 
names(design.fact) <- c("x1", "x2")
response.camelback <- apply(design.fact, 1, camelback)
response.camelback <- data.frame(response.camelback) 
names(response.camelback) <- "y" 

# model identification
fitted.model1 <- km(~1, design=design.fact, response=response.camelback, 
covtype="gauss", control=list(pop.size=50,trace=FALSE), parinit=c(0.5, 0.5))

# EI maximization
library(rgenoud)
lower <- rep(0,d) 
upper <- rep(1,d)   
oEGO <- max_EI(fitted.model1, lower, upper, 
control=list(pop.size=20, BFGSburnin=2))
print(oEGO)

# graphics
n.grid <- 20
x.grid <- y.grid <- seq(0,1,length=n.grid)
design.grid <- expand.grid(x.grid, y.grid)
response.grid <- apply(design.grid, 1, camelback)
z.grid <- matrix(response.grid, n.grid, n.grid)
contour(x.grid,y.grid,z.grid,40)
title("Camelback Function")
points(design.fact[,1], design.fact[,2], pch=17, col="blue")
points(oEGO$par[1], oEGO$par[2], pch=19, col="red")


####################################################################
### "ONE-SHOT" EI-MAXIMIZATION OF THE GOLDSTEIN-PRICE FUNCTION #####
### 	     KNOWN AT A 9-POINTS FACTORIAL DESIGN              #####
####################################################################


# a 9-points factorial design, and the corresponding response
d <- 2 
n <- 9
design.fact <- expand.grid(seq(0,1,length=3), seq(0,1,length=3)) 
names(design.fact)<-c("x1", "x2")
design.fact <- data.frame(design.fact) 
names(design.fact)<-c("x1", "x2")
response.goldsteinPrice <- apply(design.fact, 1, goldsteinPrice)
response.goldsteinPrice <- data.frame(response.goldsteinPrice) 
names(response.goldsteinPrice) <- "y" 

# model identification
fitted.model1 <- km(~1, design=design.fact, response=response.goldsteinPrice, 
covtype="gauss", control=list(pop.size=50, max.generations=50, 
wait.generations=5, BFGSburnin=10, trace=FALSE), parinit=c(0.5, 0.5), optim.method="gen")

# EI maximization
library(rgenoud)
lower <- rep(0,d); upper <- rep(1,d);     # domain for Branin function
oEGO <- max_EI(fitted.model1, lower, upper, control
=list(pop.size=50, max.generations=50, wait.generations=5, BFGSburnin=10))
print(oEGO)

# graphics
n.grid <- 20
x.grid <- y.grid <- seq(0,1,length=n.grid)
design.grid <- expand.grid(x.grid, y.grid)
response.grid <- apply(design.grid, 1, goldsteinPrice)
z.grid <- matrix(response.grid, n.grid, n.grid)
contour(x.grid,y.grid,z.grid,40)
title("Goldstein-Price Function")
points(design.fact[,1], design.fact[,2], pch=17, col="blue")
points(oEGO$par[1], oEGO$par[2], pch=19, col="red")

}

