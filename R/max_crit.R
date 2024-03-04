##' Maximization of the Expected Improvement criterion
##' 
##' For a number of \code{control$restarts}, generates a large number of random samples,
##' then picks the one with best EI value to start L-BFGS.
##' 
##' @param model an object of class \code{\link[DiceKriging]{km}} ,
##' @param type Kriging type: "SK" or "UK"
##' @param lower,upper vectors of lower and upper bounds for the variables to be optimized over,
##' @param minimization logical specifying if EI is used in minimiziation or in
##' maximization,
##' @param control optional list of control parameters for optimization. For now only the number of \code{restarts} can be set.
##' @param proxy Boolean, if TRUE, then EI maximization is replaced by the minimization of the kriging mean.
##' @param trcontrol an optional list to activate the Trust-region management (see \code{\link[DiceOptim]{TREGO.nsteps}})
##' @param n.cores Number of cores if parallel computation is used
##' @return A list with components:
##' 
##' \item{par}{The best set of parameters found.}
##' 
##' \item{value}{The value of expected improvement at par.}
##' @author Victor Picheny 
##' @keywords optimize
##' @importFrom parallel mclapply
##' @examples
##' 
##' set.seed(123)
##' library(parallel)
##' ##########################################################
##' ### "ONE-SHOT" EI-MAXIMIZATION OF THE BRANIN FUNCTION ####
##' ### 	KNOWN AT A 9-POINTS FACTORIAL DESIGN          ####
##' ##########################################################
##' 
##' # a 9-points factorial design, and the corresponding response
##' d <- 2 
##' n <- 9
##' design.fact <- expand.grid(seq(0,1,length=3), seq(0,1,length=3)) 
##' names(design.fact) <- c("x1", "x2")
##' design.fact <- data.frame(design.fact) 
##' names(design.fact) <- c("x1", "x2")
##' response.branin <- apply(design.fact, 1, branin)
##' response.branin <- data.frame(response.branin) 
##' names(response.branin) <- "y" 
##' 
##' # model identification
##' fitted.model1 <- km(~1, design=design.fact, response=response.branin, 
##' covtype="gauss", control=list(pop.size=50,trace=FALSE), parinit=c(0.5, 0.5))
##' 
##' # EGO one step
##' lower <- rep(0,d) 
##' upper <- rep(1,d)     # domain for Branin function
##' oEGO <- max_crit(fitted.model1, lower=lower, upper=upper)
##' print(oEGO)
##' 
##' # graphics
##' n.grid <- 20
##' x.grid <- y.grid <- seq(0,1,length=n.grid)
##' design.grid <- expand.grid(x.grid, y.grid)
##' response.grid <- apply(design.grid, 1, branin)
##' z.grid <- matrix(response.grid, n.grid, n.grid)
##' contour(x.grid,y.grid,z.grid,40)
##' title("Branin Function")
##' points(design.fact[,1], design.fact[,2], pch=17, col="blue")
##' points(oEGO$par[1], oEGO$par[2], pch=19, col="red")
##' @export max_crit
max_crit <-function(model, type = "UK", lower, upper, minimization = TRUE, control=NULL, proxy=FALSE, trcontrol=NULL, n.cores=1) {
  
  # Check trcontrol
  if (is.null(trcontrol$sigma)) trcontrol$sigma <- (upper - lower) / 2
  if (is.null(trcontrol$x0)) trcontrol$x0 <- model@X[which.min(model@y),]
  if (is.null(trcontrol$beta)) trcontrol$beta <- 0.7
  if (is.null(trcontrol$alpha)) trcontrol$alpha <-  1/trcontrol$beta
  if (is.null(trcontrol$kappa)) trcontrol$kappa <- 1e-4
  if (is.null(trcontrol$crit)) trcontrol$crit <- "EI"
  if (is.null(trcontrol$inTR)) trcontrol$inTR <- FALSE
  
  # Simple change of bounds if inside Trust-region
  if (trcontrol$inTR) {
    TRlower <- trcontrol$x0 - trcontrol$sigma
    TRupper <- trcontrol$x0 + trcontrol$sigma
    lower <- pmax(TRlower, lower)
    upper <- pmin(TRupper, upper)
    if (trcontrol$crit!="EI") proxy <- TRUE
  }
  d <- ncol(model@X)
  
  # Set the number of parallel searches (i.e. random search + BFGS)
  if (is.null(control$restarts))  control$restarts <- min(10, d)
  
  # To use BFGS after random search 
  if (is.null(control$BFGS)) control$BFGS <- TRUE
  
  # Sample size for the random search
  if (!is.null(control$n.ei)) N <- control$n.ei else N <- min(2e3, 500*d)
  
  if (minimization) {
    plugin <- min(model@y)
  } else {
    plugin <- -max(model@y)
  }
  
  bestval <- Inf
  bestpar <- c()
  
  my.optim <- function(seed) {
    x <- matrix(rep(lower,N), byrow=TRUE,ncol=d) + matrix(runif(d*N),ncol=d)*matrix(rep((upper-lower),N), byrow=TRUE,ncol=d)
    if(d==1) integration.points <- matrix(integration.points,ncol=1)
    EI.vect <- EI(x, model, plugin=NULL, type=type, minimization = minimization, envir=NULL, proxy=proxy)
    I <- which.max(EI.vect)
    val <- EI.vect[I]
    
    par <- x[I,,drop=FALSE]
    
    o <- list(val=val, par=par)
    
    if (control$BFGS) {
        EI.envir <- new.env()	
        environment(EI) <- environment(EI.grad) <- EI.envir 
        gr <- EI.grad
        o <- optim(par, fn=EI, gr = EI.grad, plugin=plugin, model=model,
                   type=type, minimization = TRUE, envir=EI.envir, proxy=proxy,
                   method = "L-BFGS-B", lower = lower, upper = upper,
                   control = list(maxit=100, fnscale=-1, trace=0), hessian = FALSE)
    }
    return(o)
  }
  oo <- mclapply(as.list(1:control$restarts), my.optim, mc.cores=n.cores)
  allval <- Reduce(cbind, lapply(oo, function(alist) alist$val))
  I <- which.max(allval)
  par <- as.numeric(oo[[I]]$par)
  return(list(par=par, value= oo[[I]]$val))
}
