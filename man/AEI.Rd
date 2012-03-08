\name{AEI}
\Rdversion{1}
\alias{AEI}

\title{ Augmented Expected Improvement}

\description{
Evaluation of the Augmented Expected Improvement (AEI) criterion, which is a modification of the classical EI criterion for noisy functions. 
The AEI consists of the regular EI multiplied by a penalization function that accounts for the disminishing payoff of observation replicates.
The current minimum y_min is chosen as the kriging predictor of the observation with smallest kriging quantile.}

\usage{
AEI(x, model, new.noise.var=0, y.min=NULL, type = "UK", envir=NULL)
}

\arguments{
  \item{x}{ the input vector at which one wants to evaluate the criterion }
  \item{model}{ a Kriging model of "km" class }
  \item{new.noise.var}{ the (scalar) noise variance of the future observation.}
  \item{y.min}{ The kriging predictor at the current best point (point with smallest kriging quantile). 
If not provided, this quantity is evaluated. } 
  \item{type}{ Kriging type: "SK" or "UK" }
  \item{envir}{ environment for saving intermediate calculations and reusing them within AEI.grad_optim }
}
%\details{}
\value{Augmented Expected Improvement}
\references{
D.R. Jones, M. Schonlau, and W.J. Welch (1998), Efficient global optimization of expensive black-box functions, \emph{Journal of Global Optimization}, 13, 455-492.

D. Huang, T.T. Allen, W.I. Notz, and N. Zeng (2006), Global Optimization of Stochastic Black-Box Systems via Sequential Kriging Meta-Models, \emph{Journal of Global Optimization}, 34, 441-466.
}
\author{ 
Victor Picheny (CERFACS, Toulouse, France)
David Ginsbourger (University of Bern, Switzerland)
}
%\note{}

%\seealso{}

\examples{

##########################################################################
###    AEI SURFACE ASSOCIATED WITH AN ORDINARY KRIGING MODEL        ####
### OF THE BRANIN FUNCTION KNOWN AT A 12-POINT LATIN HYPERCUBE DESIGN ####
##########################################################################

#library("lhs")
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
nt <- nrow(design.grid)
crit.grid <- rep(0,1,nt)
func.grid <- rep(0,1,nt)

for (i in 1:nt)
{ 
crit.grid[i] <- AEI(x=design.grid[i,], model=model, new.noise.var=noise.var)
func.grid[i] <- test.function(design.grid[i,])
}

# Compute kriging mean and variance on a grid
names(design.grid) <- c("V1","V2")
pred <- predict.km(model, newdata=design.grid, type="UK")
mk.grid <- pred$m
sk.grid <- pred$sd

# Plot actual function
z.grid <- matrix(func.grid, n.grid, n.grid)
filled.contour(x.grid,y.grid, z.grid, nlevels=50, color = rainbow,
plot.axes = {title("Actual function");
points(model@X[,1],model@X[,2],pch=17,col="blue"); 
axis(1); axis(2)})

# Plot Kriging mean
z.grid <- matrix(mk.grid, n.grid, n.grid)
filled.contour(x.grid,y.grid, z.grid, nlevels=50, color = rainbow,
plot.axes = {title("Kriging mean");
points(model@X[,1],model@X[,2],pch=17,col="blue"); 
axis(1); axis(2)})

# Plot Kriging variance
z.grid <- matrix(sk.grid^2, n.grid, n.grid)
filled.contour(x.grid,y.grid, z.grid, nlevels=50, color = rainbow,
plot.axes = {title("Kriging variance");
points(model@X[,1],model@X[,2],pch=17,col="blue"); 
axis(1); axis(2)})

# Plot AEI criterion
z.grid <- matrix(crit.grid, n.grid, n.grid)
filled.contour(x.grid,y.grid, z.grid, nlevels=50, color = rainbow,
plot.axes = {title("AEI");
points(model@X[,1],model@X[,2],pch=17,col="blue"); 
axis(1); axis(2)})
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory. 
\keyword{models}
