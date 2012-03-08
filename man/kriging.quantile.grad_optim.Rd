\name{kriging.quantile.grad_optim}

\alias{kriging.quantile.grad_optim}

\title{ Analytical gradient of the Kriging quantile of level alpha}

\description{ Computes the gradient of the Kriging quantile of level alpha at the current location. Only available for Universal Kriging with constant trend (Ordinary Kriging). }

\usage{
kriging.quantile.grad_optim(x, model, alpha=0.1, type="UK", envir)}

%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{x}{ a vector representing the input for which one wishes to calculate kriging.quantile.grad_optim. }

  \item{model}{ an object of class \code{\link[DiceKriging]{km}}. }

  \item{alpha}{ A quantile level (between 0 and 1) }

  \item{type}{ Kriging type: "SK" or "UK" }

  \item{envir}{ an environment specifying where to get intermediate values calculated in \code{\link{kriging.quantile}}. }
}

\value{
  The gradient of the Kriging mean predictor with respect to x. 
% Returns 0 at design points (where the gradient does not exist).
}

\author{David Ginsbourger (IMSV, University of Bern, Switzerland)  
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

############# 1-D example with noisy observations #############

fundet <- function (x) {
return((sin(10 * x)/(1 + x) + 2 * cos(5 * x) * x^3 + 0.841)/1.6)}
theta <- 1/sqrt(30)
p <- 2
x <- c(0,0.15,0.40,0.9)
n <- length(x)
t <- seq(0,1,by=0.01)
covtype="gauss"#"matern5_2"

###########

noise.var <- c(0.1,0.01,0.3,0.1)^2
z <- fundet(x); 
y <- z + sqrt(noise.var)*rnorm(length(x))
model <- km(y~1, design=data.frame(x=x), response=
data.frame(y=y), coef.trend=0, covtype=covtype, 
coef.cov=theta, coef.var=1, noise.var=noise.var, 
control=list(trace=FALSE)) 
p <- predict.km(model, newdata=t, type="SK")

par(mfrow=c(2,1))
lower <- p$lower95; upper <- p$upper95
plot(t, p$mean, type="l", ylim=c(1.1*min(c(lower,y)) , 
1.5*max(c(upper,y))), xlab="x", ylab="y",col="blue", lwd=1.5)
polygon(c(t,rev(t)), c(lower, rev(upper)), col=gray(0.9), 
border = gray(0.9))
lines(t, p$mean, type="l", ylim=c(min(lower) ,max(upper)), 
xlab="x", ylab="y",col="blue", lwd=1)
lines(t, lower , col="blue", lty=4, lwd=1.7); 
lines(t, upper , col="blue", lty=4, lwd=1.7)
lines(t,fundet(t),col="black",lwd=2);
points(x, y, pch=8,col="blue", cex=2, lwd=2)  
segments(x, y-1.96*sqrt(noise.var), x, 
y+1.96*sqrt(noise.var), col= 'red', lwd=2)
segments(x-0.2/n, y-1.96*sqrt(noise.var), x+0.2/n, 
y-1.96*sqrt(noise.var), col= 'red', lwd=2)
segments(x-0.2/n, y+1.96*sqrt(noise.var), x+0.2/n, 
y+1.96*sqrt(noise.var), col= 'red', lwd=2)

###########


plot(t, p$upper, type="l", xlab="x", ylab="y",col="blue", lwd=1.5)

quantile.envir <- new.env()	
	environment(kriging.quantile) <- environment(kriging.quantile.grad_optim) <- quantile.envir 

q <- matrix("NA", ncol=1, nrow=length(t))
gq <- matrix("NA", ncol=1, nrow=length(t))
for(i in seq(1,length(t))){
q[i] <- kriging.quantile(t[i], model, alpha=0.95, envir=quantile.envir)
gq[i] <- kriging.quantile.grad_optim(t[i], model, alpha=0.95,
envir=quantile.envir)
}

plot(t, gq, type="l", xlab="x", ylab="y",col="blue", lwd=2)
abline(h=0, col="red")

###########
}
