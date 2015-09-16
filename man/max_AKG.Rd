\name{max_AKG}
\Rdversion{1}
\alias{max_AKG}

\title{ Maximizer of the Expected Quantile Improvement criterion function }

\description{
Maximization, based on the package rgenoud of the Expected Quantile Improvement (AKG) criterion.}

\usage{
max_AKG(model, new.noise.var=0, type = "UK", lower, upper, parinit=NULL,
        control=NULL,
        cluster=FALSE
)
}

\arguments{
  \item{model}{ a Kriging model of "km" class }

  \item{new.noise.var}{ the (scalar) noise variance of an observation. Default value is 0 (noise-free observation).}  

  \item{type}{ Kriging type: "SK" or "UK" }
      
  \item{lower}{ vector containing the lower bounds of the variables to be optimized over }

  \item{upper}{ vector containing the upper bounds of the variables to be optimized over }

  \item{parinit}{ optional vector containing the initial values for the variables to be optimized over }

   \item{control}{  optional list of control parameters for optimization. 
One can control  \code{"pop.size"}  (default : [N=3*2^dim for dim<6 and N=32*dim otherwise]),  \code{"max.generations"} (12),  
\code{"wait.generations"} (2) and  \code{"BFGSburnin"} (2) of function  \code{"genoud"} (see \code{\link[rgenoud]{genoud}}). 
Numbers into brackets are the default values }  
  \item{cluster}{to invoke \code{\link{genoud}} in parallel, 
  an object of the 'cluster' class returned
  by one of the \code{\link{makeCluster}} commands in the
  \sQuote{parallel} 
  package or a vector of machine names.  See
  \code{\link{genoud}} for further details.
}
}

\value{A list with components:
\item{par}{the best set of parameters found.}
\item{value}{the value AKG at par.}
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
###    AKG SURFACE AND OPTIMIZATION PERFORMED BY GENOUD               ####
###    FOR AN ORDINARY KRIGING MODEL                                  ####
### OF THE BRANIN FUNCTION KNOWN AT A 12-POINT LATIN HYPERCUBE DESIGN ####
##########################################################################
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
for (i in 1:doe.size)  {y.tilde[i] <- test.function(doe[i,]) 
+ sqrt(noise.var)*rnorm(n=1)}
y.tilde <- as.numeric(y.tilde)

# Create kriging model
model <- km(y~1, design=doe, response=data.frame(y=y.tilde),
     covtype="gauss", noise.var=rep(noise.var,1,doe.size), 
     lower=rep(.1,dim), upper=rep(1,dim), control=list(trace=FALSE))

# Optimisation using max_AKG
res <- max_AKG(model, new.noise.var=noise.var, type = "UK", 
lower=c(0,0), upper=c(1,1)) 
X.genoud <- res$par

\dontrun{
# Compute actual function and criterion on a grid
n.grid <- 12 # Change to 21 for a nicer picture
x.grid <- y.grid <- seq(0,1,length=n.grid)
design.grid <- expand.grid(x.grid, y.grid)
names(design.grid) <- c("V1","V2")
nt <- nrow(design.grid)
crit.grid <- apply(design.grid, 1, AKG, model=model, new.noise.var=noise.var)

# # 2D plots
z.grid <- matrix(crit.grid, n.grid, n.grid)
tit <- "Green: best point found by optimizer"
filled.contour(x.grid,y.grid, z.grid, nlevels=50, color = topo.colors,
plot.axes = {title(tit);points(model@X[,1],model@X[,2],pch=17,col="blue"); 
points(X.genoud[1],X.genoud[2],pch=17,col="green");
axis(1); axis(2)})
}
}
