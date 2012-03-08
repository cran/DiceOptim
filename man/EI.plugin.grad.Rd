\name{EI.plugin.grad}

\alias{EI.plugin.grad}

\title{ Analytical gradient of the Expected Improvement criterion with plug-in}

\description{ Computes the gradient of the Expected Improvement at the current location. Only available for Universal Kriging with constant trend (Ordinary Kriging). }

\usage{
EI.plugin.grad(x, model, plugin, type="UK", envir=NULL)
}

\arguments{
  \item{x}{ the input vector at which one wants to evaluate the criterion }
  \item{model}{ a Kriging model of "km" class }
  \item{plugin}{ Plugin for the minimum of the observations. } 
  \item{type}{ Kriging type: "SK" or "UK" }
  \item{envir}{ an optional environment specifying where to assign intermediate values for future gradient calculations. Default is NULL.}
}

\value{
  The gradient of the expected improvement criterion with respect to x. 
Returns 0 at design points (where the gradient does not exist).
}
\author{
David Ginsbourger (University of Bern, Switzerland)
Victor Picheny (CERFACS, Toulouse, France)  
}


%\references{}

\seealso{ \code{\link{EI.plugin}} }

\keyword{models}
\keyword{optimize}

%
\examples{

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
for (i in 1:doe.size)  {y.tilde[i] <- test.function(doe[i,]) + sqrt(noise.var)*rnorm(n=1)}
y.tilde <- as.numeric(y.tilde)

# Create kriging model
model <- km(y~1, design=doe, response=data.frame(y=y.tilde),
            covtype="gauss", noise.var=rep(noise.var,1,doe.size), lower=rep(.1,dim), upper=rep(1,dim))

# Compute actual function and criterion on a grid
n.grid <- 21
x.grid <- y.grid <- seq(0,1,length=n.grid)
design.grid <- expand.grid(x.grid, y.grid)
nt <- nrow(design.grid)
crit.grid <- rep(0,1,nt)
crit.grad    <- matrix(0,nt,2)

EI.envir <- new.env()

for (i in 1:nt)
{ 
  crit.grid[i] <- EI.plugin(x=design.grid[i,], model=model, plugin=1, envir=EI.envir)
  crit.grad[i,] <- EI.plugin.grad(x=design.grid[i,], model=model, plugin=1, envir=EI.envir)
}

z.grid <- matrix(crit.grid, n.grid, n.grid)
contour(x.grid,y.grid, z.grid, 30)
title("EI.plugin and its gradient")
points(model@X[,1],model@X[,2],pch=17,col="blue")

for (i in 1:nt)
{
 x <- design.grid[i,]
 arrows(x$Var1,x$Var2, x$Var1+crit.grad[i,1]*.03,x$Var2+crit.grad[i,2]*.03, length=0.04,code=2,col="orange",lwd=2)
}
}
