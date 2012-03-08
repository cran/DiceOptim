\name{max_EI.plugin}
\Rdversion{1}
\alias{max_EI.plugin}

\title{ Maximizer of the Expected Improvement criterion function with plugin of the minimum. }

\description{
Maximization, based on the package rgenoud of the Expected Improvement (EI) criterion with plugin of the minimum.}

\usage{
max_EI.plugin(model, plugin, type = "UK", lower, upper, parinit=NULL, control=NULL)
}

\arguments{
  \item{model}{ a Kriging model of "km" class }
  
  \item{plugin}{ Plugin for the current minimum of the observations. } 

  \item{type}{ Kriging type: "SK" or "UK" }
      
  \item{lower}{ vector containing the lower bounds of the variables to be optimized over }

  \item{upper}{ vector containing the upper bounds of the variables to be optimized over }

  \item{parinit}{ optional vector containing the initial values for the variables to be optimized over }

   \item{control}{  optional list of control parameters for optimization. 
One can control  \code{"pop.size"}  (default : [4+3*log(nb of variables)]),  \code{"max.generations"} (5),  
\code{"wait.generations"} (2) and  \code{"BFGSburnin"} (0) of function  \code{"genoud"} (see \code{\link[rgenoud]{genoud}}). 
Numbers into brackets are the default values }  

}

\value{A list with components:
\item{par}{the best set of parameters found.}
\item{value}{the value EI at par.}
}
%\references{}
\author{ 
Victor Picheny (CERFACS, Toulouse, France) 
David Ginsbourger (University of Bern, Switzerland)
}
%\note{}

%\seealso{}

%
\examples{
set.seed(10)

# Set test problem parameters
doe.size <- 10
dim <- 2
test.function <- get("branin2")
lower <- rep(0,1,dim)
upper <- rep(1,1,dim)
noise.var <- 0.2

# Generate DOE and response
doe <- as.data.frame(optimumLHS(n=doe.size, k=dim))
y.tilde <- rep(0, 1, doe.size)
for (i in 1:doe.size) {
y.tilde[i] <- test.function(doe[i,]) + sqrt(noise.var)*rnorm(n=1)}
y.tilde <- as.numeric(y.tilde)

# Create kriging model
model <- km(y~1, design=doe, response=data.frame(y=y.tilde),
     covtype="gauss", noise.var=rep(noise.var,1,doe.size), 
     lower=rep(.1,dim), upper=rep(1,dim), control=list(trace=FALSE))

# Optimisation using max_AEI
res <- max_AEI(model, new.noise.var=noise.var, type = "UK", 
lower=c(0,0), upper=c(1,1)) 
X.genoud <- res$par

# Compute actual function and criterion on a grid
n.grid <- 21
x.grid <- y.grid <- seq(0,1,length=n.grid)
design.grid <- expand.grid(x.grid, y.grid)
names(design.grid) <- c("V1","V2")
nt <- nrow(design.grid)
crit.grid <- rep(0,1,nt)

for (i in 1:nt)
{ 
  crit.grid[i] <- AEI(x=design.grid[i,], model=model, new.noise.var=noise.var)
}

# # 2D plots
z.grid <- matrix(crit.grid, n.grid, n.grid)
tit <- c("Green: best point found by optimizer")
tit <- paste(tit, collapse="")
# tit <- "; black: best on grid"
filled.contour(x.grid,y.grid, z.grid, nlevels=50, color = rainbow,
plot.axes = {title(tit);points(model@X[,1],model@X[,2],pch=17,col="blue"); 
points(X.genoud[1],X.genoud[2],pch=17,col="green");
axis(1); axis(2)})
}
