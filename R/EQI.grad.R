##' EQI's Gradient
##' 
##' Analytical gradient of the Expected Quantile Improvement (EQI) criterion.
##' 
##' 
##' @param x the input vector at which one wants to evaluate the criterion
##' @param model a Kriging model of "km" class
##' @param new.noise.var (scalar) noise variance of the future observation.
##' Default value is 0 (noise-free observation).
##' @param beta Quantile level (default value is 0.9)
##' @param q.min Best kriging quantile. If not provided, this quantity is
##' evaluated.
##' @param type Kriging type: "SK" or "UK"
##' @param envir environment for inheriting intermediate calculations from EQI
##' @return Gradient of the Expected Quantile Improvement
##' @author Victor Picheny  
##' 
##' David Ginsbourger 
##' @examples
##' 
##' \dontrun{
##' set.seed(421)
##' 
##' # Set test problem parameters
##' doe.size <- 12
##' dim <- 2
##' test.function <- get("branin2")
##' lower <- rep(0,1,dim)
##' upper <- rep(1,1,dim)
##' noise.var <- 0.2
##' 
##' # Generate DOE and response
##' doe <- as.data.frame(matrix(runif(doe.size*dim),doe.size))
##' y.tilde <- rep(0, 1, doe.size)
##' for (i in 1:doe.size)  {
##' y.tilde[i] <- test.function(doe[i,]) + sqrt(noise.var)*rnorm(n=1)
##' }
##' y.tilde <- as.numeric(y.tilde)
##' 
##' # Create kriging model
##' model <- km(y~1, design=doe, response=data.frame(y=y.tilde),
##'         covtype="gauss", noise.var=rep(noise.var,1,doe.size), 
##' 	lower=rep(.1,dim), upper=rep(1,dim), control=list(trace=FALSE))
##' 
##' # Compute actual function and criterion on a grid
##' n.grid <- 12
##' x.grid <- y.grid <- seq(0,1,length=n.grid)
##' design.grid <- expand.grid(x.grid, y.grid)
##' nt <- nrow(design.grid)
##' 
##' crit.grid <- apply(design.grid, 1, EQI, model=model, new.noise.var=noise.var, beta=.9)
##' crit.grad <- t(apply(design.grid, 1, EQI.grad, model=model, new.noise.var=noise.var, beta=.9))
##' 
##' z.grid <- matrix(crit.grid, n.grid, n.grid)
##' contour(x.grid,y.grid, z.grid, 30)
##' title("EQI and its gradient")
##' points(model@@X[,1],model@@X[,2],pch=17,col="blue")
##' 
##' options(warn=-1)
##' for (i in 1:nt)
##' {
##'  x <- design.grid[i,]
##'  arrows(x$Var1,x$Var2, x$Var1+crit.grad[i,1]*.2,x$Var2+crit.grad[i,2]*.2, 
##' length=0.04,code=2,col="orange",lwd=2)
##' }
##' }
##' 
##' @export
EQI.grad <- function(x, model, new.noise.var=0, beta=0.9, q.min=NULL, type = "UK", envir=NULL){

########## Convert x in proper format(s) ###
d <- length(x)
if (d != model@d){ stop("x does not have the right size") }
newdata.num <- as.numeric(x)
newdata <- data.frame(t(newdata.num))
colnames(newdata) = colnames(model@X)

######### Compute q.min if missing #########
if (is.null(q.min))
{ pred <- predict.km(model, newdata=model@X, type=type, checkNames = FALSE)
  q.min <- min(pred$mean + qnorm(beta)*pred$sd)  }

######### Intermediate values ##########
T <- model@T
z <- model@z
u <- model@M

if (is.null(envir))
{
  ######### Compute prediction at x #########
  predx <- predict(model, newdata=newdata, type=type, checkNames = FALSE)
  mk.old <- predx$mean
  sk.old <- predx$sd
  v <- predx$Tinv.c
  c <- predx$c
  
  ######### mq and sq ##########
  if (sk.old < 1e-16)
  { mq <- mk.old
    sq <- 0
  } else
  { mq <- mk.old + qnorm(beta) * sqrt((new.noise.var * sk.old^2)/(new.noise.var + sk.old^2))
    sq <- sk.old^2/sqrt(new.noise.var + sk.old^2)
  } 
  
  xcr <- (q.min - mq)/sq
  xcr.prob <- as.numeric(pnorm(xcr))
  xcr.dens <- as.numeric(dnorm(xcr))
  
} else
{ toget <- matrix(c("xcr.prob", "xcr.dens", "c", "v", "mq", "sq"),1,12)
  apply(toget, 2, get, envir=envir)
  xcr.prob <- envir$xcr.prob
  xcr.dens <- envir$xcr.dens
  c        <- envir$c
  v        <- envir$v
  mq       <- envir$mq
  sq       <- envir$sq
  sk.old   <- envir$sk.old
}

if (sk.old < sqrt(model@covariance@sd2)/1e6)
{ eqi.grad.val <- rep(0,d)
} else
{
  ######### Intermediate values ##########
  F.newdata <- model.matrix(model@trend.formula, data=newdata)
  dc <- covVector.dx(x=newdata.num, X=model@X, object=model@covariance, c=c)
  f.deltax <- trend.deltax(x=newdata.num, model=model)
  W <- backsolve(t(T), dc, upper.tri=FALSE)
  
  ######### Gradient of mk and sk2 of the old model ##########
  mk.old.grad <- t(W)%*%z + t(model@trend.coef%*%f.deltax)
  if (type=="UK")
  { tuuinv <- solve(t(u)%*%u)
    sk2.old.grad <-  t( -2*t(v)%*%W + 2*(F.newdata - t(v)%*%u )%*% tuuinv %*% (f.deltax - t(t(W)%*%u) ) )
  } else
  { sk2.old.grad <-  t( -2*t(v)%*%W) }
  
  ######### Gradients of mq and sq ##########
  mq.grad  <- mk.old.grad + as.numeric(qnorm(beta)*(new.noise.var)^(3/2)/(2*sk.old*(sk.old^2 + new.noise.var)^(3/2)))*sk2.old.grad
  sq2.grad <- as.numeric((sk.old^2)*(2*new.noise.var + sk.old^2)/(new.noise.var + sk.old^2)^2)*sk2.old.grad
  sq.grad <- sq2.grad/as.numeric(2*sq)
  
  ######### Gradient of EQI ##########
  eqi.grad.val <- - xcr.prob * mq.grad + xcr.dens * sq.grad
}
return(eqi.grad.val)
}
