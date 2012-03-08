\name{AKG.grad}
\Rdversion{1}
\alias{AKG.grad}

\title{ AKG's Gradient }

\description{
Gradient of the Approximate Knowledge Gradient (AKG) criterion.
}

\usage{
AKG.grad(x, model, new.noise.var=0, type = "UK")
}

\arguments{
  \item{x}{ the input vector at which one wants to evaluate the criterion }
  \item{model}{ a Kriging model of "km" class }
  \item{new.noise.var}{ (scalar) noise variance of the future observation. Default value is 0 (noise-free observation).}  
  \item{type}{ Kriging type: "SK" or "UK" }
}
%\details{}
\value{Gradient of the Approximate Knowledge Gradient}
%\references{}
\author{
Victor Picheny (CERFACS, Toulouse, France)
}
%\note{}

%\seealso{}

\examples{
##########################################################################
###    AKG SURFACE AND ITS GRADIENT ASSOCIATED WITH AN ORDINARY       ####
###                                            KRIGING MODEL          ####
### OF THE BRANIN FUNCTION KNOWN AT A 12-POINT LATIN HYPERCUBE DESIGN ####
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
crit.grid[i] <- AKG(x=design.grid[i,], model=model, new.noise.var=noise.var)
crit.grad[i,] <- AKG.grad(x=design.grid[i,], model=model, new.noise.var=noise.var)
}

z.grid <- matrix(crit.grid, n.grid, n.grid)
contour(x.grid,y.grid, z.grid, 30)
title("AKG and its gradient")
points(model@X[,1],model@X[,2],pch=17,col="blue")

for (i in 1:nt)
{
 x <- design.grid[i,]
 arrows(x$Var1,x$Var2, x$Var1+crit.grad[i,1]*.2,x$Var2+crit.grad[i,2]*.2, 
length=0.04,code=2,col="orange",lwd=2)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory. \keyword{}
%\keyword{internal}
