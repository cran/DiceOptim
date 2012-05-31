\name{noisy.optimizer}
\Rdversion{1}
\alias{noisy.optimizer}

\title{Optimization of homogenously noisy functions based on Kriging}

\description{
Sequential optimization of kriging-based criterion conditional on noisy observations, with model update after each evaluation. 
Eight criteria are proposed to choose the next observation: random search, sequential parameter optimization (SPO), reinterpolation, Expected Improvement (EI) with plugin, Expected Quantile Improvement (EQI), 
quantile minimization, Augmented Expected Improvement (AEI) and Approximate Knowledge Gradient (AKG).
The criterion optimization is based on the package rgenoud.
}

\usage{
noisy.optimizer(optim.crit, optim.param=NULL, model, n.ite, noise.var=NULL, funnoise, lower, upper, parinit=NULL, control=NULL,CovReEstimate=TRUE,
                             noiseReEstimate=FALSE, nugget.LB=1e-5, obs.n.rep=NULL, estim.model=NULL)
}

\arguments{

  \item{optim.crit}{String defining the criterion to be optimized at each iteration. Possible values are: "random.search", "SPO", "reinterpolation", "EI.plugin", "EQI", "min.quantile", "AEI", "AKG".}

  \item{optim.param}{List of parameters for the chosen criterion.
 
  For "SPO": optim.param$nrep is a vector of integers representing the number of repetitions for each observation; 
                 optim.param$history.best is a vector of the indices (integers) of the best observations. 
                 If not provided, nrep is set to a vector of ones and history.best is set to the index of the current best observation.

  For "EI.plugin": optim.param$plugin.type is a string defining which plugin is to be used. Possible values are "ytilde", "quantile" and "other".
                      If "quantile" is chosen, optim.param$quantile defines the quantile level. 
                      If "other" is chosen, optim.param$plugin directly sets the plugin value.

  For "EQI": optim.param$quantile defines the quantile level. If not provided, default value is 0.9.

  For "min.quantile": optim.param$quantile defines the quantile level. If not provided, default value is 0.1.

  For "AEI": optim.param$quantile defines the quantile level to choose the current best point. If not provided, default value is 0.75.
  }
    
  \item{model}{ a Kriging model of "km" class }

  \item{n.ite}{ Number of iterations}
  
  \item{noise.var}{ Noise variance (scalar). If noiseReEstimate=TRUE, it is an initial guess for the unknown variance 
(used in optimization).}
    
  \item{funnoise}{ objective (noisy) function }
    
  \item{lower}{ vector containing the lower bounds of the variables   to be optimized over }

  \item{upper}{ vector containing the upper bounds of the variables to be optimized over }

  \item{parinit}{ optional vector of initial values for the variables to be optimized over }

   \item{control}{  optional list of control parameters for optimization. 
One can control  \code{"pop.size"}  (default : [N=3*2^dim for dim<6 and N=32*dim otherwise]]),  \code{"max.generations"} (N),  
\code{"wait.generations"} (2) and  \code{"BFGSburnin"} (0) of function  \code{"genoud"} (see \code{\link[rgenoud]{genoud}}). 
Numbers into brackets are the default values }  
	
   \item{CovReEstimate}{ optional boolean specfiying if the covariance parameters should be re-estimated at every iteration (default value = TRUE)}

   \item{noiseReEstimate}{optional boolean specfiying if the noise variance should be re-estimated at every iteration (default value = FALSE)}

   \item{obs.n.rep}{optional vector containing the number of repetitions per observation point. 
Required if "model" has heterogeneous variances and noise is reestimated. In that case noise variances
must be of the form sigma^2/obs.n.rep[i].}

   \item{nugget.LB}{optional scalar of minimal value for the estimated noise variance. Default value is 1e-5.}

   \item{estim.model}{optional kriging model of "km" class with homogeneous nugget effect (no noise.var). 
Required if noise variance is reestimated and the initial "model" has heterogenenous noise variances.}

}


\value{A list with components:
\item{par}{ The added observations}
\item{value}{the added observation values}
\item{nsteps}{the number of added observations}
\item{lastmodel}{the current (last) kriging model of "km" class}
\item{estim.model}{If noiseReEstimate=TRUE, the current (last) kriging model of "km" class for estimating the noise variance.}
}
%\references{}
\author{ 
Victor Picheny (CERFACS, Toulouse, France)
}
%\note{}

%\seealso{}

\examples{
##########################################################################
### EXAMPLE 1: 5 OPTIMIZATION STEPS USING EQI WITH KNOWN NOISE         ###
### AND KNOWN COVARIANCE PARAMETERS FOR THE BRANIN FUNCTION            ###
##########################################################################

set.seed(10)

# Set test problem parameters
doe.size <- 9
dim <- 2
test.function <- get("branin2")
lower <- rep(0,1,dim)
upper <- rep(1,1,dim)
noise.var <- 0.1

# Build noisy simulator
funnoise <- function(x)
{     f.new <- test.function(x) + sqrt(noise.var)*rnorm(n=1)
      return(f.new)}

# Generate DOE and response
doe <- as.data.frame(optimumLHS(n=doe.size, k=dim))
y.tilde <- funnoise(doe)

# Create kriging model
model <- km(y~1, design=doe, response=data.frame(y=y.tilde),
     covtype="gauss", noise.var=rep(noise.var,1,doe.size), 
     lower=rep(.1,dim), upper=rep(1,dim), control=list(trace=FALSE))

# Optimisation with noisy.optimizer
n.ite <- 15
optim.param <- list()
optim.param$quantile <- .9
optim.result <- noisy.optimizer(optim.crit="EQI", optim.param=optim.param, model=model,
		n.ite=n.ite, noise.var=noise.var, funnoise=funnoise, lower=lower, upper=upper,
		noiseReEstimate=FALSE, CovReEstimate=FALSE)

new.model <- optim.result$model
new.doe   <- optim.result$X
new.y     <- optim.result$y

##### DRAW RESULTS #####
# Compute actual function on a grid
n.grid <- 21
x.grid <- y.grid <- seq(0,1,length=n.grid)
design.grid <- expand.grid(x.grid, y.grid)
names(design.grid) <- c("V1","V2")
nt <- nrow(design.grid)
func.grid <- rep(0,1,nt)

for (i in 1:nt)
{ func.grid[i] <- test.function(x=design.grid[i,])}

# Compute initial and final kriging on a grid
pred <- predict(object=model, newdata=design.grid, type="UK", checkNames = FALSE)
mk.grid1 <- pred$m
sk.grid1 <- pred$sd

pred <- predict(object=new.model, newdata=design.grid, type="UK", checkNames = FALSE)
mk.grid2 <- pred$m
sk.grid2 <- pred$sd

# Plot initial kriging mean
z.grid <- matrix(mk.grid1, n.grid, n.grid)
filled.contour(x.grid,y.grid, z.grid, nlevels=50, color = topo.colors,
plot.axes = {title("Initial kriging mean");
points(model@X[,1],model@X[,2],pch=17,col="black"); 
axis(1); axis(2)})

# Plot initial kriging variance
z.grid <- matrix(sk.grid1^2, n.grid, n.grid)
filled.contour(x.grid,y.grid, z.grid, nlevels=50, color = topo.colors,
plot.axes = {title("Initial kriging variance");
points(model@X[,1],model@X[,2],pch=17,col="black"); 
axis(1); axis(2)})

# Plot final kriging mean
z.grid <- matrix(mk.grid2, n.grid, n.grid)
filled.contour(x.grid,y.grid, z.grid, nlevels=50, color = topo.colors,
plot.axes = {title("Final kriging mean");
points(new.model@X[,1],new.model@X[,2],pch=17,col="black"); 
axis(1); axis(2)})

# Plot final kriging variance
z.grid <- matrix(sk.grid2^2, n.grid, n.grid)
filled.contour(x.grid,y.grid, z.grid, nlevels=50, color = topo.colors,
plot.axes = {title("Final kriging variance");
points(new.model@X[,1],new.model@X[,2],pch=17,col="black"); 
axis(1); axis(2)})

# Plot actual function and observations
z.grid <- matrix(func.grid, n.grid, n.grid)
tit <- "Actual function - Black: initial points; red: added points"
filled.contour(x.grid,y.grid, z.grid, nlevels=50, color = topo.colors,
plot.axes = {title(tit);points(model@X[,1],model@X[,2],pch=17,col="black"); 
points(new.doe[1,],new.doe[2,],pch=15,col="red");
axis(1); axis(2)})
} 


