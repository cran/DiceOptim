\name{kriging.quantile.grad}

\alias{kriging.quantile.grad}

\title{ Analytical gradient of the Kriging quantile of level alpha}

\description{ Computes the gradient of the Kriging quantile of level alpha at the current location. Only available for Universal Kriging with constant trend (Ordinary Kriging). }

\usage{
kriging.quantile.grad(x, model, alpha=0.1, type="UK")}

%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{x}{ a vector representing the input for which one wishes to calculate kriging.quantile.grad. }

  \item{model}{ an object of class \code{\link[DiceKriging]{km}}. }

  \item{alpha}{ A quantile level (between 0 and 1) }

  \item{type}{ Kriging type: "SK" or "UK" }
}

\value{
  The gradient of the Kriging mean predictor with respect to x. 
% Returns 0 at design points (where the gradient does not exist).
}

\author{David Ginsbourger (IMSV, University of Berne, Switzerland)  
}


\references{

O. Roustant, D. Ginsbourger, Y. Deville,
\emph{DiceKriging, DiceOptim: Two R packages for the analysis of computer experiments by kriging-based metamodeling and optimization}, 
submitted to J. Stat. Soft., 2010.
\url{http://hal.archives-ouvertes.fr/hal-00495766_v2/}

D. Ginsbourger (2009), \emph{Multiples metamodeles pour l'approximation et l'optimisation de fonctions numeriques multivariables}, Ph.D. thesis, Ecole Nationale Superieure des
Mines de Saint-Etienne, 2009. 
\url{http://www.ginsbourger.ch/recherche/these.php}

}

\seealso{ \code{\link{EI.grad}} }

\keyword{models}
\keyword{optimize}

\examples{
##########################################################################
###    KRIGING QUANTILE SURFACE AND ITS GRADIENT FOR                  ####
###    THE BRANIN FUNCTION KNOWN AT A 12-POINT LATIN HYPERCUBE DESIGN ####
##########################################################################

set.seed(421)

# Set test problem parameters
doe.size <- 12
dim <- 2
test.function <- get("branin2")
lower <- rep(0,1,dim)
upper <- rep(1,1,dim)
noise.var <- 0.2

# Generate DOE and response
doe <- as.data.frame(optimumLHS(n=doe.size, k=dim))
y.tilde <- rep(0, 1, doe.size)
for (i in 1:doe.size)  {
y.tilde[i] <- test.function(doe[i,]) + sqrt(noise.var)*rnorm(n=1)
}
y.tilde <- as.numeric(y.tilde)

# Create kriging model
model <- km(y~1, design=doe, response=data.frame(y=y.tilde),
        covtype="gauss", noise.var=rep(noise.var,1,doe.size), 
	lower=rep(.1,dim), upper=rep(1,dim), control=list(trace=FALSE))

# Compute actual function and criterion on a grid
n.grid <- 21
x.grid <- y.grid <- seq(0,1,length=n.grid)
design.grid <- expand.grid(x.grid, y.grid)
#names(design.grid) <- c("V1", "V2")
nt <- nrow(design.grid)
crit.grid <- rep(0,1,nt)
crit.grad    <- matrix(0,nt,2)

for (i in 1:nt)
{ 
crit.grid[i] <- kriging.quantile(x=design.grid[i,], model=model, alpha=.1)
crit.grad[i,] <- kriging.quantile.grad(x=design.grid[i,], model=model, alpha=.1)
}

z.grid <- matrix(crit.grid, n.grid, n.grid)
contour(x.grid,y.grid, z.grid, 30)
title("kriging.quantile and its gradient")
points(model@X[,1],model@X[,2],pch=17,col="blue")

for (i in 1:nt)
{
 x <- design.grid[i,]
 arrows(x$Var1,x$Var2, x$Var1+crit.grad[i,1]*.01,x$Var2+crit.grad[i,2]*.01, 
length=0.04,code=2,col="orange",lwd=2)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory. \keyword{}
%\keyword{internal}
