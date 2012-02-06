\name{qEI}
\alias{qEI}

\title{ Monte-Carlo estimation of the multipoints Expected Improvement criterion (noise-free version)}

\description{ Computes the Expected Improvement at any vector of new locations. }

\usage{
qEI(newdata, model, type="UK", MC.samples=10000, return.I=FALSE) 
}

\arguments{

  \item{newdata}{ a matrix representing the set of input vectors (columns) where to estimate qEI. }

  \item{model}{ an object of class \code{km}. }

  \item{type}{ "SK" or "UK", distinguishing the cases of known or unknown trend parameters. }

  \item{MC.samples}{ integer, size of the Monte-Carlo sample to be used within the estimation of qEI. }

  \item{return.I}{ logical, standing for returning improvements.matrix in the list of outputs or not }
}

\value{
  The multipoints Expected Improvement, defined as 

 \deqn{EI(X_{new}) := E[ ( min(Y(X)) - min(Y( X_{new} )) )^{+} | Y(X)=y(X)],} 

where \eqn{X} is the current design of experiments, \eqn{ X_{new} } is a new candidate design, and \eqn{Y} is a random process assumed to have generated the objective function \eqn{y}.
}



\author{David Ginsbourger (Departement of Mathematics and Statistics, University of Berne, Switzerland)  

Olivier Roustant (Ecole Nationale Superieure des Mines de Saint-Etienne, France).}


\references{

D. Ginsbourger (2009), \emph{Multiples metamodeles pour l'approximation et l'optimisation de fonctions numeriques multivariables}, Ph.D. thesis, Ecole Nationale Superieure des
Mines de Saint-Etienne, 2009. 
\url{http://members.unine.ch/david.ginsbourger/recherche/these.htm}

D. Ginsbourger, R. Le Riche, and L. Carraro (2009), chapter "Kriging is well-suited to parallelize optimization", to appear in \emph{Computational Intelligence in Expensive Optimization Problems}, Studies in Evolutionary Learning and Optimization, Springer.

J. Mockus (1988), \emph{Bayesian Approach to Global Optimization}. Kluwer academic publishers.

B.D. Ripley (1987), \emph{Stochastic Simulation}, Wiley.

T.J. Santner, B.J. Williams, and W.J. Notz (2003), \emph{The design and analysis of computer experiments}, Springer.

M. Schonlau (1997), \emph{Computer experiments and global optimization}, Ph.D. thesis, University of Waterloo.

}

\seealso{ \code{\link{max_qEI.CL}}, \code{\link{CL.nsteps}} }

\keyword{models}
%\keyword{optimize}

\examples{
set.seed(123)
# ----------------------------------
# A 1D example with known parameters
# ----------------------------------

x <- c(0, 0.4, 0.6, 0.8, 1)
y <- 10*c(-0.6, 0, -2, 0.5, 0.9)

theta <- 0.1 
sigma <- 10 
trend <- 5*c(-2,1)
model <- km(~x,design=data.frame(x=x), response=data.frame(y=y), 
coef.trend=trend, covtype="gauss", coef.cov=theta, coef.var=sigma^2, control <- list(trace=FALSE))

## Kriging with gaussian covariance, and linear trend t(x) = -10 + 5x
type="SK"
t <- seq(from=0, to=1, by=0.005)
p <- predict(model, newdata=t, type=type)
plot(t, p$mean, type="l", ylim=c(min(p$lower95),max(p$upper95)), 
xlab="x", ylab="y")
lines(t, p$lower95, col="black", lty=2)
lines(t, p$upper95, col="black", lty=2)
points(x, y, col="red", pch=19)
abline(h=0)

## A first try of qEI
candidate.design <- seq(0,1,,101)
res <- qEI(newdata=candidate.design, model=model, type=type, 
MC.samples=10000, return.I=TRUE)

## One_point EI can easly be estimated on the basis of the matrix of improvements
EI_simples <- colMeans(res$I)
plot(candidate.design,EI_simples, type="l")

## Let us check that the analytical formula of EI and the MC estimate match
EI_analytical <- apply(as.data.frame(candidate.design), 1, EI, model, type=type)
plot(EI_simples, EI_analytical)

# Computation of the 2-points EI grid associated to candidate.design
# (Straightforward computation method: faster not available in this version yet) 

two_points_EI <- matrix(0,ncol=length(candidate.design),
nrow=length(candidate.design))

for(i in seq(1,length(candidate.design)) ){
for(j in seq(i,length(candidate.design)) ){
qI <- pmax(res$I[,i],res$I[,j])
two_points_EI[i,j] <- mean(qI)
two_points_EI[j,i] <- two_points_EI[i,j]
}}

# Plotting the 2_points EI as a function of both points, 
# based on the latter qEI estimation

contour(two_points_EI)
}
