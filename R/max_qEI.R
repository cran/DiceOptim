##' Maximization of multipoint expected improvement criterion (qEI)
##' 
##' Maximization of the \code{\link{qEI}} criterion. Two options are available
##' : Constant Liar (CL), and brute force qEI maximization with
##' Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm, or GENetic Optimization
##' Using Derivative (genoud) algorithm.
##' 
##' - CL is a heuristic method. First, the regular Expected Improvement EI is
##' maximized (\code{\link{max_EI}}). Then, for the next points, the Expected
##' Improvement is maximized again, but with an artificially updated Kriging
##' model. Since the response values corresponding to the last best point
##' obtained are not available, the idea of CL is to replace them by an
##' arbitrary constant value L (a "lie") set by the user (default is the
##' minimum of all currently available observations).
##' 
##' - The BFGS algorithm is implemented in the standard function
##' \code{\link{optim}}. Analytical formulae of \code{\link{qEI}} and its
##' gradient \code{\link{qEI.grad}} are used. The \code{nStarts} starting
##' points are by default sampled with respect to the regular EI
##' (\code{\link{sampleFromEI}}) criterion.
##' 
##' - The "genoud" method calls the function \code{\link{genoud}} using
##' analytical formulae of \code{\link{qEI}} and its gradient
##' \code{\link{qEI.grad}}.
##' 
##' The parameters of list \code{optimcontrol} are :
##' 
##' - \code{optimcontrol$method} : "BFGS" (default), "genoud" ; a string
##' specifying the method used to maximize the criterion (irrelevant when
##' \code{crit} is "CL" because this method always uses genoud),
##' 
##' - when \code{crit="CL"} :
##' 
##' + \code{optimcontrol$parinit} : optional matrix of initial values (must
##' have model@@d columns, the number of rows is not constrained),
##' 
##' + \code{optimcontrol$L} : "max", "min", "mean" or a scalar value specifying
##' the liar ; "min" takes \code{model@@min}, "max" takes \code{model@@max},
##' "mean" takes the prediction of the model ; When L is \code{NULL}, "min" is
##' taken if \code{minimization==TRUE}, else it is "max".
##' 
##' + The parameters of function \code{\link{genoud}}. Main parameters are :
##' \code{"pop.size"} (default : [N=3*2^model@@d for dim<6 and N=32*model@@d
##' otherwise]), \code{"max.generations"} (default : 12),
##' \code{"wait.generations"} (default : 2) and \code{"BFGSburnin"} (default :
##' 2).
##' 
##' - when \code{optimcontrol$method = "BFGS"} :
##' 
##' + \code{optimcontrol$nStarts} (default : 4),
##' 
##' + \code{optimcontrol$fastCompute} : if TRUE (default), a fast approximation
##' method based on a semi-analytic formula is used, see [Marmin 2014] for
##' details,
##' 
##' + \code{optimcontrol$samplingFun} : a function which sample a batch of
##' starting point (default : \code{\link{sampleFromEI}}),
##' 
##' + \code{optimcontrol$parinit} : optional 3d-array of initial (or candidate)
##' batches (for all \code{k}, parinit[,,k] is a matrix of size
##' \code{npoints*model@@d} representing one batch). The number of initial
##' batches (length(parinit[1,1,])) is not contrained and does not have to be
##' equal to \code{nStarts}. If there is too few initial batches for
##' \code{nStarts}, missing batches are drawn with \code{samplingFun} (default
##' : \code{NULL}),
##' 
##' - when \code{optimcontrol$method = "genoud"} :
##' 
##' + \code{optimcontrol$fastCompute} : if TRUE (default), a fast approximation
##' method based on a semi-analytic formula is used, see [Marmin 2014] for
##' details,
##' 
##' + \code{optimcontrol$parinit} : optional matrix of candidate starting
##' points (one row corresponds to one point),
##' 
##' + The parameters of the \code{\link{genoud}} function. Main parameters are
##' \code{"pop.size"} (default : \code{[50*(model@@d)*(npoints)]}),
##' \code{"max.generations"} (default : 5), \code{"wait.generations"} (default
##' : 2), \code{"BFGSburnin"} (default : 2).
##' 
##' 
##' @param model an object of class \code{\link[DiceKriging]{km}} ,
##' @param npoints an integer representing the desired number of iterations,
##' @param lower vector of lower bounds,
##' @param upper vector of upper bounds,
##' @param crit "exact", "CL" : a string specifying the criterion used. "exact"
##' triggers the maximization of the multipoint expected improvement at each
##' iteration (see \code{\link{max_qEI}}), "CL" applies the Constant Liar
##' heuristic,
##' @param minimization logical specifying if the qEI to be maximized is used
##' in minimiziation or in maximization,
##' @param optimcontrol an optional list of control parameters for
##' optimization. See details.
##' @return A list with components:
##' 
##' \item{par}{A matrix containing the \code{npoints} input vectors found.}
##' 
##' \item{value}{A value giving the qEI computed in \code{par}.}
##' @author Sebastien Marmin 
##' 
##' Clement Chevalier 
##' 
##' David Ginsbourger 
##' @seealso \code{\link{qEI}}, \code{\link{qEI.grad}}
##' @references
##' 
##' C. Chevalier and D. Ginsbourger (2014) Learning and Intelligent
##' Optimization - 7th International Conference, Lion 7, Catania, Italy,
##' January 7-11, 2013, Revised Selected Papers, chapter Fast computation of
##' the multipoint Expected Improvement with applications in batch selection,
##' pages 59-69, Springer.
##' 
##' D. Ginsbourger, R. Le Riche, L. Carraro (2007), A Multipoint Criterion for
##' Deterministic Parallel Global Optimization based on Kriging. The
##' International Conference on Non Convex Programming, 2007.
##' 
##' D. Ginsbourger, R. Le Riche, and L. Carraro. Kriging is well-suited to
##' parallelize optimization (2010), In Lim Meng Hiot, Yew Soon Ong, Yoel
##' Tenne, and Chi-Keong Goh, editors, \emph{Computational Intelligence in
##' Expensive Optimization Problems}, Adaptation Learning and Optimization,
##' pages 131-162. Springer Berlin Heidelberg.
##' 
##' J. Mockus (1988), \emph{Bayesian Approach to Global Optimization}. Kluwer
##' academic publishers.
##' 
##' M. Schonlau (1997), \emph{Computer experiments and global optimization},
##' Ph.D. thesis, University of Waterloo.
##' @keywords optimize
##' @examples
##' 
##' 
##' \dontrun{
##' 
##' set.seed(000)
##' # 3-points EI maximization.
##' # 9-points factorial design, and the corresponding response
##' d <- 2
##' n <- 9
##' design.fact <- expand.grid(seq(0,1,length=3), seq(0,1,length=3)) 
##' names(design.fact)<-c("x1", "x2")
##' design.fact <- data.frame(design.fact) 
##' names(design.fact)<-c("x1", "x2")
##' response.branin <- apply(design.fact, 1, branin)
##' response.branin <- data.frame(response.branin) 
##' names(response.branin) <- "y" 
##' lower <- c(0,0)
##' upper <- c(1,1)
##' 
##' # number of point in the bacth
##' batchSize <- 3
##' 
##' # model identification
##' fitted.model <- km(~1, design=design.fact, response=response.branin, 
##'                    covtype="gauss", control=list(pop.size=50,trace=FALSE), parinit=c(0.5, 0.5))
##' 
##' # maximization of qEI
##' 
##' # With a multistarted BFGS algorithm
##' maxBFGS <- max_qEI(model = fitted.model, npoints = batchSize, lower = lower, upper = upper, 
##' crit = "exact",optimcontrol=list(nStarts=3,method = "BFGS"))
##' # With a genetic algorithme using derivatives
##' maxGen  <- max_qEI(model = fitted.model, npoints = batchSize, lower = lower, upper = upper, 
##' crit = "exact", optimcontrol=list(nStarts=3,method = "genoud",pop.size=100,max.generations = 15))
##' # With the constant liar heuristic
##' maxCL   <- max_qEI(model = fitted.model, npoints = batchSize, lower = lower, upper = upper, 
##' crit = "CL",optimcontrol=list(pop.size=20))
##' 
##' # comparison
##' print(maxBFGS$value)
##' print(maxGen$value)
##' print(maxCL$value)
##' 
##' }
##' 
##' 
##' @export max_qEI
max_qEI <- function(model, npoints, lower, upper, crit="exact", minimization = TRUE, optimcontrol=NULL) {  
  if (is.null(optimcontrol$method)) optimcontrol$method <- "BFGS"
  optim.method <- optimcontrol$method
  d <- model@d 
  parinit <- optimcontrol$parinit
   if (crit == "CL") {
    res <- max_qEI.CL(model, npoints, optimcontrol$L, lower, upper, parinit=parinit,minimization = minimization, control=optimcontrol)
    res <- list(par = res$par, value = matrix(qEI(x = res$par, model = model)))
   } else if(crit=="exact") {
     EI.envir <- new.env()
     environment(qEI) <- environment(qEI.grad) <- EI.envir
     
     LOWER <- c(apply(matrix(lower,d,1),1,rep,npoints))
     UPPER <- c(apply(matrix(upper,d,1),1,rep,npoints))
      if (optim.method == "BFGS") {
        if (is.null(parinit))
          parinit <- array(NaN,dim=c(npoints,model@d,0))
        if (is.null(optimcontrol$nStarts))
          optimcontrol$nStarts <- 4
        if (is.null(optimcontrol$fastCompute))
          optimcontrol$fastCompute <- TRUE
        if (is.null(optimcontrol$sampleFun))
          optimcontrol$sampleFun <- sampleFromEI
        if (is.null(optimcontrol$gradNum))
          optimcontrol$gradNum <- FALSE
        if (is.null(optimcontrol$maxit))
          optimcontrol$maxit <- 100
        maxi <- 0
        startPoints <- optimcontrol$sampleFun(model=model,minimization=TRUE,n = (npoints*optimcontrol$nStarts),lower=lower,upper=upper)
        for (i in 1:(optimcontrol$nStarts)) {
          if (i > length(parinit[1,1,])){
            x <- startPoints[((i-1)*npoints+1):(i*(npoints)),]
          } else {
            x <- parinit[,,i]
          }
          if (!(optimcontrol$gradNum)) {
            o <- optim(c(x), qEI ,qEI.grad, control=list(trace = 0,REPORT=1,fnscale = -1,maxit=optimcontrol$maxit),method="L-BFGS-B",lower=LOWER,upper=UPPER,model = model,envir=EI.envir,fastCompute=optimcontrol$fastCompute,minimization = minimization)
          } else {
            o <- optim(par = c(x), fn = qEI ,gr=NULL, control=list(trace = 0,REPORT=1,fnscale = -1,maxit=optimcontrol$maxit),method="L-BFGS-B",lower=LOWER,upper=UPPER,model = model,envir=EI.envir,fastCompute=optimcontrol$fastCompute,minimization=minimization)
          }
          par <- matrix(o$par,ncol=d)
          value <- as.matrix(o$value)
          if (value>=maxi) {
            maxi <- value
            parmax <- par
          }
        }
        colnames(parmax) <- colnames(model@X)
        colnames(maxi) <- "EI"
        res <- list(par = parmax, value = maxi)
        #res <- max_qEI.BFGS(model,npoints,lower,upper, parinit = NULL,minimization = minimization, control = control)
      } else if (optim.method =="genoud") {
        if (is.null(optimcontrol$pop.size)) optimcontrol$pop.size <- 50*d
        if (is.null(optimcontrol$max.generations)) optimcontrol$max.generations <- 5
        if (is.null(optimcontrol$wait.generations)) optimcontrol$wait.generations <- 2
        if (is.null(optimcontrol$BFGSburnin)) optimcontrol$BFGSburnin <- 2
        if (is.null(optimcontrol$fastCompute)) optimcontrol$fastCompute <- TRUE
        if (is.null(optimcontrol$print.level)) optimcontrol$print.level <- 0
        domaine <- cbind(LOWER,UPPER)
        
        o <- genoud(qEI, nvars = (d*npoints), max = TRUE,
                    gr = qEI.grad,
                    pop.size = optimcontrol$pop.size, 
                    max.generations = optimcontrol$max.generations, 
                    wait.generations = optimcontrol$wait.generations, 
                    hard.generation.limit = TRUE, starting.values = c(optimcontrol$parinit), 
                    MemoryMatrix = TRUE, Domains = domaine, default.domains = 10, 
                    solution.tolerance = 1e-09, boundary.enforcement = 2, 
                    lexical = FALSE, gradient.check = FALSE, BFGS = TRUE, 
                    data.type.int = FALSE, hessian = FALSE, 
                    unif.seed = floor(runif(1, max = 10000)), 
                    int.seed = floor(runif(1, max = 10000)),
                    print.level = 0, share.type = 0, instance.number = 0, 
                    output.path = "stdout", output.append = FALSE, project.path = NULL, 
                    P1 = 50, P2 = 50, P3 = 50, P4 = 50, P5 = 50, P6 = 50, 
                    P7 = 50, P8 = 50, P9 = 0, P9mix = NULL, BFGSburnin = optimcontrol$BFGSburnin, 
                    BFGSfn = NULL, BFGShelp = NULL, cluster = FALSE, balance = FALSE, 
                    debug = FALSE, model = model,fastCompute = optimcontrol$fastCompute,minimization=minimization,envir=EI.envir)
        
        o$par <- matrix(o$par,ncol=d)
        colnames(o$par) <- colnames(model@X)
        o$value <- as.matrix(o$value)
        colnames(o$value) <- "EI"
        res <- list(par = o$par, value = o$value)
        #res <- max_qEI.gen(model,npoints,lower,upper, parinit = NULL,minimization = minimization, control = control)
      } else {
        stop(paste(paste("\'",optim.method,"\'",sep=""),"is unknown."))
      }
   } else {
     stop(paste(paste("\'",crit,"\'",sep=""),"is unknown."))
   }
   return(res)
}



max_qEI.CL <- function(model, npoints, L, lower, upper, parinit=NULL, minimization = TRUE, control=NULL) {
  KB <- FALSE
  n1 <- nrow(model@X)
  lengthParinit <- length(parinit[1,])
  if (is.null(L)) {
    liar <- minimization*min(model@y) + (!minimization)*max(model@y)
  } else if (L == "max") {
    liar <- max(model@y)
  } else if (L == "min") {
    liar <- min(model@y)
  } else if (L == "mean") {
    KB <- TRUE
  } else {
    if ((!is.numeric(L))||(length(L)!=1)) stop("control$L must be NULL, \"max\", \"min\", \"mean\" or a scalar specifying the plugin.")
    liar <- L
  }
  for (s in 1:npoints) {
    if (s > lengthParinit){
      startPoint <- NULL
    } else {
      startPoint <- parinit[,s]
    }
    oEGO <- max_EI(model=model, lower=lower, upper=upper, parinit=startPoint, minimization = minimization, control=c(control,list(print.level=0)))
  		model@X <- rbind(model@X, oEGO$par)
    if (KB) liar <- predict(object=model,newdata=oEGO$par,se.compute=FALSE,cov.compute=FALSE,light.return=TRUE,checkNames=FALSE)$mean
    model@y <- rbind(model@y, liar, deparse.level=0)   		
    model@F <- trendMatrix.update(model, Xnew=data.frame(oEGO$par))	
    if (model@noise.flag) {
	# heterogenous case : use 0 nugget for new points
	model@noise.var = c(model@covariance@nugget, 0)
    }	
    model <- computeAuxVariables(model)
  }
  return(list(par = model@X[(n1+1):(n1+npoints),, drop=FALSE], value = model@y[(n1+1):(n1+npoints),, drop=FALSE]))
}
