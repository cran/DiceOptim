##' Analytical gradient of the Expected Improvement criterion
##' 
##' Computes the gradient of the Expected Improvement at the current location.
##' The current minimum of the observations can be replaced by an arbitrary
##' value (plugin), which is usefull in particular in noisy frameworks.
##' 
##' 
##' @param x a vector representing the input for which one wishes to calculate
##' \code{\link{EI}}.
##' @param model an object of class \code{\link[DiceKriging]{km}}.
##' @param plugin optional scalar: if provided, it replaces the minimum of the
##' current observations,
##' @param type Kriging type: "SK" or "UK"
##' @param minimization logical specifying if EI is used in minimiziation or in
##' maximization,
##' @param envir an optional environment specifying where to get intermediate
##' values calculated in \code{\link{EI}}.
##' @return The gradient of the expected improvement criterion with respect to
##' x.  Returns 0 at design points (where the gradient does not exist).
##' @author David Ginsbourger 
##' 
##' Olivier Roustant 
##' 
##' Victor Picheny 
##' @seealso \code{\link{EI}}
##' @references
##' 
##' D. Ginsbourger (2009), \emph{Multiples metamodeles pour l'approximation et
##' l'optimisation de fonctions numeriques multivariables}, Ph.D. thesis, Ecole
##' Nationale Superieure des Mines de Saint-Etienne, 2009.
##' \url{http://www.ginsbourger.ch/recherche/these.php}
##' 
##' J. Mockus (1988), \emph{Bayesian Approach to Global Optimization}. Kluwer
##' academic publishers.
##' 
##' T.J. Santner, B.J. Williams, and W.J. Notz (2003), \emph{The design and
##' analysis of computer experiments}, Springer.
##' 
##' M. Schonlau (1997), \emph{Computer experiments and global optimization},
##' Ph.D. thesis, University of Waterloo.
##' @keywords models optimize
##' @examples
##' 
##' \dontrun{
##' set.seed(123)
##' # a 9-points factorial design, and the corresponding response
##' d <- 2; n <- 9
##' design.fact <- expand.grid(seq(0,1,length=3), seq(0,1,length=3)) 
##' names(design.fact)<-c("x1", "x2")
##' design.fact <- data.frame(design.fact) 
##' names(design.fact)<-c("x1", "x2")
##' response.branin <- apply(design.fact, 1, branin)
##' response.branin <- data.frame(response.branin) 
##' names(response.branin) <- "y" 
##' 
##' # model identification
##' fitted.model1 <- km(~1, design=design.fact, response=response.branin, 
##' covtype="gauss", control=list(pop.size=50,trace=FALSE), parinit=c(0.5, 0.5))
##' 
##' # graphics
##' n.grid <- 50
##' x.grid <- y.grid <- seq(0,1,length=n.grid)
##' design.grid <- expand.grid(x.grid, y.grid)
##' #response.grid <- apply(design.grid, 1, branin)
##' EI.grid <- apply(design.grid, 1, EI,fitted.model1)
##' #EI.grid <- apply(design.grid, 1, EI.plot,fitted.model1, gr=TRUE)
##' 
##' z.grid <- matrix(EI.grid, n.grid, n.grid)
##' 
##' contour(x.grid,y.grid,z.grid,20)
##' title("Expected Improvement for the Branin function known at 9 points")
##' points(design.fact[,1], design.fact[,2], pch=17, col="blue")
##' 
##' # graphics
##' n.gridx <- 15
##' n.gridy <- 20
##' x.grid2 <- seq(0,1,length=n.gridx) 
##' y.grid2 <- seq(0,1,length=n.gridy) 
##' design.grid2 <- expand.grid(x.grid2, y.grid2)
##' 
##' EI.envir <- new.env()	
##' 	environment(EI) <- environment(EI.grad) <- EI.envir 
##' 
##' options(warn=-1)
##' for(i in seq(1, nrow(design.grid2)) )
##' {
##' 	x <- design.grid2[i,]
##' 	ei <- EI(x, model=fitted.model1, envir=EI.envir)
##' 	eigrad <- EI.grad(x , model=fitted.model1, envir=EI.envir)
##' 	if(!(is.null(ei)))
##' 	{
##' 	arrows(x$Var1,x$Var2,
##' 	x$Var1 + eigrad[1]*2.2*10e-5, x$Var2 + eigrad[2]*2.2*10e-5, 
##' 	length = 0.04, code=2, col="orange", lwd=2)
##' 	}
##' }
##' }
##' 
##' @export EI.grad
EI.grad <- function(x, model, plugin=NULL, type="UK", minimization = TRUE, envir=NULL){ 
	
  ########################################################################################
  if (is.null(plugin)){ 
    if (minimization) {
      plugin <- min(model@y)
    } else {
      plugin <- -max(model@y)
    }
  }
  m <- plugin
  
  ########################################################################################
  # Convert x in proper format(s)
  d <- length(x)
  if (d != model@d){ stop("x does not have the right size") }
  newdata.num <- as.numeric(x)
  newdata <- data.frame(t(newdata.num))
  colnames(newdata) = colnames(model@X)
  ########################################################################################
  # Get quantities related to the model
  T <- model@T
  X <- model@X
  z <- model@z
  u <- model@M
  covStruct <- model@covariance
  
  # Get quantities related to the prediction
  if (is.null(envir))
  {  
    predx <- predict(object=model, newdata=newdata, type=type, checkNames = FALSE,se.compute=TRUE,cov.compute=FALSE)
     kriging.mean <- predx$mean
     if(!minimization) kriging.mean <- -kriging.mean
     kriging.sd <- predx$sd
     v <- predx$Tinv.c
     c <- predx$c
  
     xcr <- (m - kriging.mean)/kriging.sd
     xcr.prob <- pnorm(xcr)
     xcr.dens <- dnorm(xcr)    
  } else
  {  # If uploaded through "envir", no prediction computation is necessary 
     toget <- matrix(c("xcr", "xcr.prob", "xcr.dens", "kriging.sd", "c", "Tinv.c"), 1, 6)
     apply(toget, 2, get, envir=envir)
     xcr        <- envir$xcr
     xcr.prob   <- envir$xcr.prob
     xcr.dens   <- envir$xcr.dens
     kriging.sd <- envir$kriging.sd
     c          <- envir$c
     v          <- envir$Tinv.c
  }

  F.newdata <- model.matrix(model@trend.formula, data=newdata)

  ########################################################################################
  # Pursue calculation only if standard deviation is non-zero
  if ( kriging.sd/sqrt(model@covariance@sd2) < 1e-06) 
  { ei.grad <- rep(0,d)
  } else 
  { # Compute derivatives of the covariance and trend functions
    dc <- covVector.dx(x=newdata.num, X=X, object=covStruct, c=c)  
    f.deltax <- trend.deltax(x=newdata.num, model=model)
  
    # Compute gradients of the kriging mean and variance
    W <- backsolve(t(T), dc, upper.tri=FALSE)
    kriging.mean.grad <- t(W)%*%z + t(model@trend.coef%*%f.deltax)
    if (!minimization) kriging.mean.grad <- -kriging.mean.grad
    if (type=="UK")
    { tuuinv <- solve(t(u)%*%u)
      kriging.sd2.grad <-  t( -2*t(v)%*%W +
                            2*(F.newdata - t(v)%*%u )%*% tuuinv %*%
                            (f.deltax - t(t(W)%*%u) ))
    } else
    { kriging.sd2.grad <-  t( -2*t(v)%*%W) }
    
    kriging.sd.grad <- kriging.sd2.grad / (2*kriging.sd)

    # Compute gradient of EI
    ei.grad <- - kriging.mean.grad * xcr.prob + kriging.sd.grad * xcr.dens
  }
  ########################################################################################  
  return(ei.grad)
}