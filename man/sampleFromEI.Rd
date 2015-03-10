\name{sampleFromEI}
\alias{sampleFromEI}

\title{ Sampling points according to the expected improvement criterion }

\description{   
Samples \code{n} points from a distribution proportional to the 
expected improvement (EI) computed from a \code{km} object.
}

\usage{
sampleFromEI(model, minimization = TRUE, n = 1, initdistrib = NULL, 
lower = rep(0,model@d), upper = rep(1,model@d), T = NULL)
}



\arguments{

  \item{model}{  an object of class \code{\link[DiceKriging]{km}}, }

  \item{minimization}{ logical specifying if EI is used in minimiziation or in maximization, }

  \item{n}{ number of points to be sampled, }

  \item{initdistrib}{ matrix of candidate points. If \code{NULL}, (\code{model@d*1000}) points are obtained by latin hypercube sampling, }

  \item{lower}{ vector of lower bounds, }

  \item{upper}{ vector of upper bounds, }

  \item{T}{ optional scalar : if provided, it replaces the current minimum (or maximum) of observations. }
  
}	



\value{ A \code{n*model@d} matrix containing the sampled points. }




\author{Sebastien Marmin (Departement of Mathematics and Statistics, University of Berne, Switzerland ; Ecole Centrale de Marseille, France ; Institut de Radioprotection et de Surete Nucleaire, France),

Clement Chevalier (Institute of Mathematics, University of Zurich, Switzerland),

David Ginsbourger (Departement of Mathematics and Statistics, University of Berne, Switzerland).
}



\references{

D.R. Jones, M. Schonlau, and W.J. Welch (1998), Efficient global optimization of expensive black-box functions, \emph{Journal of Global Optimization}, 13, 455-492.

}


\seealso{ \code{\link{EI}}, \code{\link{km}}, \code{\link{qEI}} }

\keyword{optimization}

\examples{

\dontrun{

set.seed(004)

# a 9-points factorial design, and the corresponding responses
d <- 2
n <- 9
design.fact <- expand.grid(seq(0,1,length=3), seq(0,1,length=3)) 
names(design.fact)<-c("x1", "x2")
design.fact <- data.frame(design.fact) 
names(design.fact)<-c("x1", "x2")
response.branin <- apply(design.fact, 1, branin)
response.branin <- data.frame(response.branin) 
lower <- c(0,0)
upper <- c(1,1)
names(response.branin) <- "y" 


# model identification
fitted.model <- km(~1, design=design.fact, response=response.branin, 
                   covtype="gauss", control=list(pop.size=50,trace=FALSE), parinit=c(0.5, 0.5))

# sample a 30 point batch
batchSize <- 30
x <- sampleFromEI(model = fitted.model, n = batchSize, lower = lower, upper = upper)

# graphics 
# displays the EI criterion, the design points in black and the EI-sampled points in red.
nGrid <- 15
gridAxe1 <- seq(lower[1],upper[1],length=nGrid)
gridAxe2 <- seq(lower[2],upper[2],length=nGrid)
grid <- expand.grid(gridAxe1,gridAxe2)
aa <- apply(grid,1,EI,model=fitted.model)
myMat <- matrix(aa,nrow=nGrid)
image(x = gridAxe1, y = gridAxe2, z = myMat, 
      col = colorRampPalette(c("darkgray","white"))(5*10), 
      ylab = names(design.fact)[1], xlab=names(design.fact)[2], 
      main = "Sampling from the expected improvement criterion", 
      axes = TRUE, zlim = c(min(myMat), max(myMat)))
contour(x = gridAxe1, y = gridAxe2, z = myMat, 
        add = TRUE, nlevels = 10)
points(x[,1],x[,2],pch=19,col='red')
points(fitted.model@X[,1],fitted.model@X[,2],pch=19)
}
}
