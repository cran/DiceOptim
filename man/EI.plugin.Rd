\name{EI.plugin}
\Rdversion{1}
\alias{EI.plugin}

\title{ Expected Improvement with plugin}

\description{
Evaluation of the Expected Improvement criterion, where the current minimum of the observations (y.min) is replaced by an arbitrary value (plugin). Suitable for noisy frameworks.}

\usage{
EI.plugin(x, model, plugin=NULL, type="UK", envir=NULL)
}

\arguments{
  \item{x}{ the input vector at which one wants to evaluate the criterion }
  \item{model}{ a Kriging model of "km" class }
  \item{plugin}{ Plugin for the minimum of the observations. } 
  \item{type}{ Kriging type: "SK" or "UK" }
  \item{envir}{ an optional environment specifying where to assign intermediate values for future gradient calculations. Default is NULL.}
}
%\details{}
\value{
 The Expected Improvement with plug-in, defined as 
  \deqn{EI(x) := E[( m - Y(x))^{+} | Y(X)+\epsilon=\tilde{y}],} 
  where X is the current design of experiments, m is the plugin value for the current minimum (here called ymin.plugin), and Y is the random process assumed to have generated the objective function y.
}
%\references{}
\author{ 
Victor Picheny (CERFACS, Toulouse, France) 
David Ginsbourger (University of Bern, Switzerland)
}
%\note{}

%\seealso{}

\examples{

##########################################################################
###   EI.plugin SURFACE ASSOCIATED WITH AN ORDINARY KRIGING MODEL     ####
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
nt <- nrow(design.grid)
crit.grid <- rep(0,1,nt)
func.grid <- rep(0,1,nt)

for (i in 1:nt)
{ 
crit.grid[i] <- EI.plugin(x=design.grid[i,], model=model, plugin=1)
func.grid[i] <- test.function(design.grid[i,])[1,]
}

# Compute kriging mean and variance on a grid
names(design.grid) <- c("V1","V2")
pred <- predict(model, newdata=design.grid, type="UK", checkNames = FALSE)
mk.grid <- pred$m
sk.grid <- pred$sd

# Plot actual function
z.grid <- matrix(func.grid, n.grid, n.grid)
filled.contour(x.grid,y.grid, z.grid, nlevels=50, color =topo.colors,
plot.axes = {title("Actual function");
points(model@X[,1],model@X[,2],pch=17,col="blue"); 
axis(1); axis(2)})

# Plot Kriging mean
z.grid <- matrix(mk.grid, n.grid, n.grid)
filled.contour(x.grid,y.grid, z.grid, nlevels=50, color =topo.colors,
plot.axes = {title("Kriging mean");
points(model@X[,1],model@X[,2],pch=17,col="blue"); 
axis(1); axis(2)})

# Plot Kriging variance
z.grid <- matrix(sk.grid^2, n.grid, n.grid)
filled.contour(x.grid,y.grid, z.grid, nlevels=50, color =topo.colors,
plot.axes = {title("Kriging variance");
points(model@X[,1],model@X[,2],pch=17,col="blue"); 
axis(1); axis(2)})


# Plot EI.plugin criterion
z.grid <- matrix(crit.grid, n.grid, n.grid)
filled.contour(x.grid,y.grid, z.grid, nlevels=50, color =topo.colors,
plot.axes = {title("EI.plugin");
points(model@X[,1],model@X[,2],pch=17,col="blue"); 
axis(1); axis(2)})
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory. \keyword{}
%\keyword{internal}
