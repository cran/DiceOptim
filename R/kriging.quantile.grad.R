##' Analytical gradient of the Kriging quantile of level beta
##' 
##' Computes the gradient of the Kriging quantile of level beta at the current
##' location. Only available for Universal Kriging with constant trend
##' (Ordinary Kriging).
##' 
##' 
##' @param x a vector representing the input for which one wishes to calculate
##' kriging.quantile.grad.
##' @param model an object of class \code{\link[DiceKriging]{km}}.
##' @param beta A quantile level (between 0 and 1)
##' @param type Kriging type: "SK" or "UK"
##' @param envir environment for inheriting intermediate calculations from
##' \code{"kriging.quantile"}
##' @return The gradient of the Kriging mean predictor with respect to x.  %
##' Returns 0 at design points (where the gradient does not exist).
##' @author Victor Picheny 
##' 
##' David Ginsbourger 
##' @seealso \code{\link{EI.grad}}
##' @references
##' 
##' O. Roustant, D. Ginsbourger, Y. Deville, \emph{DiceKriging, DiceOptim: Two
##' R packages for the analysis of computer experiments by kriging-based
##' metamodeling and optimization}, submitted to J. Stat. Soft., 2010.
##' \url{http://hal.archives-ouvertes.fr/hal-00495766_v2/}
##' 
##' D. Ginsbourger (2009), \emph{Multiples metamodeles pour l'approximation et
##' l'optimisation de fonctions numeriques multivariables}, Ph.D. thesis, Ecole
##' Nationale Superieure des Mines de Saint-Etienne, 2009.
##' \url{http://www.ginsbourger.ch/recherche/these.php}
##' @keywords models optimize
##' @examples
##' 
##' \dontrun{
##' ##########################################################################
##' ###    KRIGING QUANTILE SURFACE AND ITS GRADIENT FOR                  ####
##' ###    THE BRANIN FUNCTION KNOWN AT A 12-POINT LATIN HYPERCUBE DESIGN ####
##' ##########################################################################
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
##' n.grid <- 12 # Change to 21 for a nicer picture
##' x.grid <- y.grid <- seq(0,1,length=n.grid)
##' design.grid <- expand.grid(x.grid, y.grid)
##' nt <- nrow(design.grid)
##' 
##' crit.grid <- apply(design.grid, 1, kriging.quantile, model=model, beta=.1)
##' crit.grad <- t(apply(design.grid, 1, kriging.quantile.grad, model=model, beta=.1))
##' 
##' z.grid <- matrix(crit.grid, n.grid, n.grid)
##' contour(x.grid,y.grid, z.grid, 30)
##' title("kriging.quantile and its gradient")
##' points(model@@X[,1],model@@X[,2],pch=17,col="blue")
##' 
##' for (i in 1:nt)
##' {
##'  x <- design.grid[i,]
##'  arrows(x$Var1,x$Var2, x$Var1+crit.grad[i,1]*.01,x$Var2+crit.grad[i,2]*.01, 
##' length=0.04,code=2,col="orange",lwd=2)
##' }
##' }
##' 
##' @export
kriging.quantile.grad <- function(x, model, beta=0.1, type="UK", envir=NULL)
{
  ########## Convert x in proper format(s) ###
  d <- length(x)
  if (d != model@d){ stop("x does not have the right size") }
  newdata.num <- as.numeric(x)
  newdata <- data.frame(t(newdata.num))
  colnames(newdata) = colnames(model@X)

	T <- model@T
	X <- model@X
	z <- model@z
  u <- model@M
  
	covStruct <- model@covariance
  
	if (is.null(envir))
	{ predx <- predict(object=model, newdata=newdata, type=type, checkNames = FALSE)
	  mk <- predx$mean
	  sk <- predx$sd
	  c <- predx$c 
	  v <- predx$Tinv.c
	} else
	{ toget <- matrix(c("mk", "sk",  "c", "Tinv.c"),1,4)
	  apply(toget, 2, get, envir=envir)
	  c   <- envir$c
	  v   <- envir$Tinv.c
	  mk  <- envir$mk
	  sk  <- envir$sk
	}
  F.newdata <- model.matrix(model@trend.formula, data=newdata)
  
	# Compute derivatives of the covariance and trend functions
	dc <- covVector.dx(x=newdata.num, X=X, object=covStruct, c=c)  
	f.deltax <- trend.deltax(x=newdata.num, model=model)
	
	# Compute gradients of the kriging mean and variance
	W <- backsolve(t(T), dc, upper.tri=FALSE)
	mk.grad <- t(W)%*%z + t(model@trend.coef%*%f.deltax)
	
  if (sk < sqrt(model@covariance@sd2)/1e6)
  { sk.grad <- 0
  } else
  { 
    if (type=="UK")
    { tuuinv <- solve(t(u)%*%u)
	    sk2.grad <-  t( -2*t(v)%*%W +
	             2*(F.newdata - t(v)%*%u )%*% tuuinv %*%
	                          (f.deltax - t(t(W)%*%u) ))
    } else { sk2.grad <-  t( -2*t(v)%*%W)}
	  
	  sk.grad <- sk2.grad / (2*sk)
  }

  return(quantile.grad <- mk.grad + qnorm(beta)*sk.grad)
}
