#' Maximizer of the Augmented Expected Improvement criterion function
#' 
#' Maximization, based on the package rgenoud of the Augmented Expected
#' Improvement (AEI) criterion.
#' 
#' 
#' @param model a Kriging model of "km" class
#' @param new.noise.var the (scalar) noise variance of the new observation.
#' @param y.min The kriging mean prediction at the current best point (point
#' with smallest kriging quantile).  If not provided, this quantity is
#' evaluated inside the AEI function (may increase computational time).
#' @param type Kriging type: "SK" or "UK"
#' @param lower vector containing the lower bounds of the variables to be
#' optimized over
#' @param upper optional vector containing the upper bounds of the variables
#' to be optimized over
#' @param parinit optional vector containing the initial values for the
#' variables to be optimized over
#' @param control optional list of control parameters for optimization.  One
#' can control \code{"pop.size"} (default : [N=3*2^dim for dim<6 and N=32*dim
#' otherwise]), \code{"max.generations"} (12), \code{"wait.generations"} (2)
#' and \code{"BFGSburnin"} (2) of function \code{"genoud"} (see
#' \code{\link[rgenoud]{genoud}}).  Numbers into brackets are the default
#' values
#' @return A list with components: \item{par}{the best set of parameters
#' found.} \item{value}{the value AEI at par.}
#' @author Victor Picheny 
#' 
#' David Ginsbourger 
#' @examples
#' 
#' library(DiceDesign)
#' set.seed(100)
#' 
#' # Set test problem parameters
#' doe.size <- 10
#' dim <- 2
#' test.function <- get("branin2")
#' lower <- rep(0,1,dim)
#' upper <- rep(1,1,dim)
#' noise.var <- 0.2
#' 
#' # Generate DOE and response
#' doe <- as.data.frame(lhsDesign(doe.size, dim)$design)
#' y.tilde <- rep(0, 1, doe.size)
#' for (i in 1:doe.size)  {y.tilde[i] <- test.function(doe[i,]) 
#' + sqrt(noise.var)*rnorm(n=1)}
#' y.tilde <- as.numeric(y.tilde)
#' 
#' # Create kriging model
#' model <- km(y~1, design=doe, response=data.frame(y=y.tilde),
#'      covtype="gauss", noise.var=rep(noise.var,1,doe.size), 
#'      lower=rep(.1,dim), upper=rep(1,dim), control=list(trace=FALSE))
#' 
#' # Optimisation using max_AEI
#' res <- max_AEI(model, new.noise.var=noise.var, type = "UK", 
#' lower=c(0,0), upper=c(1,1)) 
#' X.genoud <- res$par
#' 
#' # Compute actual function and criterion on a grid
#' n.grid <- 12 # Change to 21 for a nicer picture
#' x.grid <- y.grid <- seq(0,1,length=n.grid)
#' design.grid <- expand.grid(x.grid, y.grid)
#' names(design.grid) <- c("V1","V2")
#' nt <- nrow(design.grid)
#' crit.grid <- apply(design.grid, 1, AEI, model=model, new.noise.var=noise.var)
#' 
#' \dontrun{
#' # # 2D plots
#' z.grid <- matrix(crit.grid, n.grid, n.grid)
#' tit <- "Green: best point found by optimizer"
#' filled.contour(x.grid,y.grid, z.grid, nlevels=50, color = rainbow,
#' plot.axes = {title(tit);points(model@@X[,1],model@@X[,2],pch=17,col="blue"); 
#' points(X.genoud[1],X.genoud[2],pch=17,col="green");
#' axis(1); axis(2)})
#' }
#' 
#' @export
max_AEI <-function(model, new.noise.var=0, y.min=NULL, type = "UK", lower, upper, parinit=NULL, control=NULL) {

  # Compute y.min if missing
  if (is.null(y.min))
  { pred <- predict(object=model, newdata=model@X, type=type, checkNames = FALSE)
    y.min <- pred$mean[which.min(pred$mean + qnorm(0.75)*pred$sd)] }
  
  AEI.envir <- new.env()
	environment(AEI) <- environment(AEI.grad) <- AEI.envir 
	gr = AEI.grad

	d <- ncol(model@X)
	if (is.null(control$print.level)) control$print.level <- 1
	if(d<=6) N <- 3*2^d else N <- 32*d 
  if (is.null(control$pop.size))  control$pop.size <- N
  if (is.null(control$BFGSmaxit)) control$BFGSmaxit <- N
	if (is.null(control$solution.tolerance))  control$solution.tolerance <- 1e-21
  if (is.null(control$max.generations))  control$max.generations <- 12
  if (is.null(control$wait.generations))  control$wait.generations <- 2
  if (is.null(control$BFGSburnin)) control$BFGSburnin <- 2
	if (is.null(parinit))  parinit <- lower + runif(d) * (upper - lower)

	domaine <- cbind(lower, upper)

	o <- genoud(AEI, nvars=d, max=TRUE, 
	            pop.size=control$pop.size, max.generations=control$max.generations, wait.generations=control$wait.generations,
	            hard.generation.limit=TRUE, starting.values=parinit, MemoryMatrix=TRUE, 
	            Domains=domaine, default.domains=10, solution.tolerance=control$solution.tolerance,
	            gr=gr, boundary.enforcement=2, lexical=FALSE, gradient.check=FALSE, BFGS=TRUE,
	            data.type.int=FALSE, hessian=FALSE, unif.seed=floor(runif(1,max=10000)), int.seed=floor(runif(1,max=10000)), 
	            print.level=control$print.level, 
	            share.type=0, instance.number=0, output.path="stdout", output.append=FALSE, project.path=NULL,
	            P1=50, P2=50, P3=50, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0, P9mix=NULL, 
	            BFGSburnin=control$BFGSburnin, BFGSfn=NULL, BFGShelp=NULL, control=list("maxit"=control$BFGSmaxit), 
	            cluster=FALSE, balance=FALSE, debug=FALSE, 
              model=model, new.noise.var=new.noise.var, y.min=y.min, type=type, envir=AEI.envir 
		)
                            
    o$par <- t(as.matrix(o$par))
	colnames(o$par) <- colnames(model@X)
	o$value <- as.matrix(o$value)
	colnames(o$value) <- "AEI"   
	return(list(par=o$par, value=o$value)) 
}
