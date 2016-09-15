##' Maximizer of the Expected Quantile Improvement criterion function
##' 
##' Maximization, based on the package rgenoud of the Expected Quantile
##' Improvement (AKG) criterion.
##' 
##' 
##' @param model a Kriging model of "km" class
##' @param new.noise.var the (scalar) noise variance of an observation. Default
##' value is 0 (noise-free observation).
##' @param type Kriging type: "SK" or "UK"
##' @param lower vector containing the lower bounds of the variables to be
##' optimized over
##' @param upper vector containing the upper bounds of the variables to be
##' optimized over
##' @param parinit optional vector containing the initial values for the
##' variables to be optimized over
##' @param control optional list of control parameters for optimization.  One
##' can control \code{"pop.size"} (default : [N=3*2^dim for dim<6 and N=32*dim
##' otherwise]), \code{"max.generations"} (12), \code{"wait.generations"} (2)
##' and \code{"BFGSburnin"} (2) of function \code{"genoud"} (see
##' \code{\link[rgenoud]{genoud}}).  Numbers into brackets are the default
##' values
##' @return A list with components: \item{par}{the best set of parameters
##' found.} \item{value}{the value AKG at par.}
##' @author Victor Picheny 
##' 
##' David Ginsbourger 
##' @examples
##' 
##' ##########################################################################
##' ###    AKG SURFACE AND OPTIMIZATION PERFORMED BY GENOUD               ####
##' ###    FOR AN ORDINARY KRIGING MODEL                                  ####
##' ### OF THE BRANIN FUNCTION KNOWN AT A 12-POINT LATIN HYPERCUBE DESIGN ####
##' ##########################################################################
##' set.seed(10)
##' 
##' # Set test problem parameters
##' doe.size <- 10
##' dim <- 2
##' test.function <- get("branin2")
##' lower <- rep(0,1,dim)
##' upper <- rep(1,1,dim)
##' noise.var <- 0.2
##' 
##' # Generate DOE and response
##' library(DiceDesign)
##' doe <- as.data.frame(lhsDesign(doe.size, dim)$design)
##' y.tilde <- rep(0, 1, doe.size)
##' for (i in 1:doe.size)  {y.tilde[i] <- test.function(doe[i,]) 
##' + sqrt(noise.var)*rnorm(n=1)}
##' y.tilde <- as.numeric(y.tilde)
##' 
##' # Create kriging model
##' model <- km(y~1, design=doe, response=data.frame(y=y.tilde),
##'      covtype="gauss", noise.var=rep(noise.var,1,doe.size), 
##'      lower=rep(.1,dim), upper=rep(1,dim), control=list(trace=FALSE))
##' 
##' # Optimisation using max_AKG
##' res <- max_AKG(model, new.noise.var=noise.var, type = "UK", 
##' lower=c(0,0), upper=c(1,1)) 
##' X.genoud <- res$par
##' 
##' \dontrun{
##' # Compute actual function and criterion on a grid
##' n.grid <- 12 # Change to 21 for a nicer picture
##' x.grid <- y.grid <- seq(0,1,length=n.grid)
##' design.grid <- expand.grid(x.grid, y.grid)
##' names(design.grid) <- c("V1","V2")
##' nt <- nrow(design.grid)
##' crit.grid <- apply(design.grid, 1, AKG, model=model, new.noise.var=noise.var)
##' 
##' # # 2D plots
##' z.grid <- matrix(crit.grid, n.grid, n.grid)
##' tit <- "Green: best point found by optimizer"
##' filled.contour(x.grid,y.grid, z.grid, nlevels=50, color = topo.colors,
##' plot.axes = {title(tit);points(model@@X[,1],model@@X[,2],pch=17,col="blue"); 
##' points(X.genoud[1],X.genoud[2],pch=17,col="green");
##' axis(1); axis(2)})
##' }
##' 
##' @export
max_AKG <-function(model, new.noise.var=0, type = "UK", lower, upper, parinit=NULL, control=NULL) {

	AKG.envir <- new.env()	
	environment(AKG) <- environment(AKG.grad) <- AKG.envir 
	gr = AKG.grad

	d <- ncol(model@X)
	if (is.null(control$print.level)) control$print.level <- 1
	if(d<=6) N <- 3*2^d else N <- 32*d 
	if (is.null(control$BFGSmaxit)) control$BFGSmaxit <- N
	if (is.null(control$pop.size))  control$pop.size <- N
	if (is.null(control$solution.tolerance))  control$solution.tolerance <- 1e-21
	if (is.null(control$max.generations))  control$max.generations <- 12
	if (is.null(control$wait.generations))  control$wait.generations <- 2
	if (is.null(control$BFGSburnin)) control$BFGSburnin <- 2
	if (is.null(parinit))  parinit <- lower + runif(d) * (upper - lower)
     
	domaine <- cbind(lower, upper)

	o <- genoud(AKG, nvars=d, max=TRUE,
	            pop.size=control$pop.size, max.generations=control$max.generations, wait.generations=control$wait.generations,
	            hard.generation.limit=TRUE, starting.values=parinit, MemoryMatrix=TRUE, 
	            Domains=domaine, default.domains=10, solution.tolerance=control$solution.tolerance,
              gr=gr, boundary.enforcement=2, lexical=FALSE, gradient.check=FALSE, BFGS=TRUE,
	            data.type.int=FALSE, hessian=FALSE, unif.seed=floor(runif(1,max=10000)), int.seed=floor(runif(1,max=10000)),
	            print.level=control$print.level, 
              share.type=0, instance.number=0, output.path="stdout", output.append=FALSE, project.path=NULL,
	            P1=50, P2=50, P3=50, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0,	P9mix=NULL, 
              BFGSburnin=control$BFGSburnin, BFGSfn=NULL, BFGShelp=NULL, control=list("maxit"=control$BFGSmaxit),
	            cluster=FALSE, balance=FALSE, debug=TRUE, 
	            model=model, new.noise.var=new.noise.var, type=type, envir=AKG.envir
		)
                            
    o$par <- t(as.matrix(o$par))
	colnames(o$par) <- colnames(model@X)
	o$value <- as.matrix(o$value)
	colnames(o$value) <- "AKG"   
	return(list(par=o$par, value=o$value)) 
}
