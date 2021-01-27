##' User-friendly wrapper of the functions \code{\link[DiceOptim]{fastEGO.nsteps}} and \code{\link[DiceOptim]{TREGO.nsteps}}.
##' Generates initial DOEs and kriging models (objects of class \code{\link[DiceKriging]{km}}), 
##' and executes \code{nsteps} iterations of either EGO or TREGO.
##' 
##' @details Does not require knowledge on kriging models (objects of class \code{\link[DiceKriging]{km}})  \cr
##' The \code{control} argument is a list that can supply any of the following components: \cr
##' \itemize{
##' \item \code{trace}{: between -1 and 3}
##' \item \code{seed}{: to fix the seed of the run}
##' \item \code{cov.reestim}{: Boolean, if TRUE (default) the covariance parameters are re-estimated at each iteration}
##' \item \code{model.trend}{: trend for the GP model}
##' \item \code{lb, ub}{: lower and upper bounds for the GP covariance ranges}
##' \item \code{nugget}{: optional nugget effect}
##' \item \code{covtype}{: covariance of the GP model (default "matern5_2")}
##' \item \code{optim.method}{: optimisation of the GP hyperparameters (default "BFGS")}
##' \item \code{multistart}{: number of restarts of BFGS}
##' \item \code{gpmean.trick, gpmean.freq}{: Boolean and integer, resp., for the gpmean trick}
##' \item \code{scaling}{: Boolean, activates input scaling}
##' \item \code{warping}{: Boolean, activates output warping}
##' \item \code{TR}{: Boolean, activates TREGO instead of EGO}
##' \item \code{trcontrol}{: list of parameters of the trust region, see \code{\link[DiceOptim]{TREGO.nsteps}}}
##' \item \code{always.sample}{: Boolean, activates force observation even if it leads to poor conditioning}
##' }
##' 
##' @param  fun scalar function to be minimized,
##' @param  budget total number of calls to the objective and constraint functions,
##' @param  lower vector of lower bounds for the variables to be optimized over,
##' @param  upper vector of upper bounds for the variables to be optimized over,
##' @param  X initial design of experiments. If not provided, X is taken as a maximin LHD with budget/3 points
##' @param  y initial set of objective observations \eqn{f(X)}. Computed if not provided.
##' @param control an optional list of control parameters. See "Details".
##' @param n.cores number of cores for parallel computation
##' @param ... additional parameters to be given to \code{fun}
##' @export
##' @return
##' A list with components:
##' \itemize{
##' \item{\code{par}}{: the best feasible point}
##' \item{\code{values}}{: a vector of the objective and constraints at the point given in \code{par},}
##' \item{\code{history}}{: a list containing all the points visited by the algorithm (\code{X}) and their corresponding objectives (\code{y}).}
##' \item{\code{model}}{: the last GP model, class \code{\link[DiceKriging]{km}}}
##' \item{\code{control}}{: full list of control values, see "Details"}
##' \item{\code{res}}{: the output of either \code{\link[DiceOptim]{fastEGO.nsteps}} or \code{\link[DiceOptim]{TREGO.nsteps}}}
##' }
##' 
##' 
##' @author
##' Victor Picheny 
##' 
##' @references
##' D.R. Jones, M. Schonlau, and W.J. Welch (1998), Efficient global
##' optimization of expensive black-box functions, \emph{Journal of Global
##' Optimization}, 13, 455-492.
##' 
##' @importFrom DiceDesign maximinESE_LHS lhsDesign
##' 
##' @examples
##' library(parallel)
##' library(DiceOptim)
##' set.seed(123)
##' 
##' #########################################################
##' ### 	10 ITERATIONS OF TREGO ON THE BRANIN FUNCTION, ####
##' ###	 STARTING FROM A 9-POINTS FACTORIAL DESIGN       ####
##' ########################################################
##' 
##' # a 9-points factorial design, and the corresponding response
##' ylim=NULL
##' fun <- branin; d <- 2
##' budget <- 5*d
##' lower <- rep(0,d)
##' upper <- rep(1,d)
##' n.init <- 2*d
##' 
##' control <- list(n.init=2*d, TR=TRUE, nugget=1e-5, trcontrol=list(algo="TREGO"), multistart=1)
##' 
##' res1 <- easyEGO(fun=fun, budget=budget, lower=lower, upper=upper, control=control, n.cores=1)
##' 
##' par(mfrow=c(3,1))
##' y <- res1$history$y
##' steps <- res1$res$all.steps
##' success <- res1$res$all.success
##' sigma <- res1$res$all.sigma
##' ymin <- cummin(y)
##' pch <- rep(1, length(sigma))
##' col <- rep("red", length(sigma))
##' pch[which(!steps)] <- 2
##' col[which(success)] <- "darkgreen"
##' 
##' pch2 <- c(rep(3, n.init), pch)
##' col2 <- c(rep("black", n.init), col)
##' plot(y, col=col2, ylim=ylim, pch=pch2, lwd=2, xlim=c(0, budget))
##' lines(ymin, col="darkgreen")
##' abline(v=n.init+.5)
##' 
##' plot((n.init+1):budget, sigma, xlim=c(0, budget), ylim=c(0, max(sigma)), 
##' pch=pch, col=col, lwd=2, main="TR size") 
##' lines((n.init+1):budget, sigma, xlim=c(0, budget)) 
##' abline(v=n.init+.5)
##' 
##' plot(NA, xlim=c(0, budget), ylim=c(0, 1), main="x0 (coordinates)")
##' for (i in 1:d) {
##'   lines((n.init+1):budget, res1$res$all.x0[,i]) 
##'   points((n.init+1):budget, res1$res$all.x0[,i], pch=pch, col=col, lwd=2) 
##' }
##' abline(v=n.init+.5)
##' 
##' par(mfrow=c(1,1))
##' pairs(res1$model@X, pch=pch2, col=col2)
##' 
easyEGO <- function (fun, budget, lower, upper, X=NULL, y=NULL, control=list(trace=1, seed=42), n.cores=1, ...) {
  
  if (is.null(control$trace)) control$trace   <- 1
  if (is.null(control$seed)) control$seed   <- 42
  
  n.init <- control$n.init
  d <- length(lower)
  
  if (length(lower) != length(upper)) {
    cat("Bound values lower and upper are not consistent. Both should be vectors of size d.")
    return(0)
  }
    
  if (!is.null(X)) {
    design.init <- data.frame(x=X)
    temp <- dim(design.init)
    
    if (temp[2] != d) {
      cat("Bound values (lower and upper) and initial DOE (X) are not consistent. \n 
          lower and upper should be vectors of size d and \n 
          X either a matrix with d columns or a data frame of d variables.")
      return(0)
    }
    n.init <- temp[1]
  } else {
    if (is.null(n.init)) {
      n.init <- min(10*d, max(3*d, round(budget/4)))
    }
    design.init <- data.frame(x=lower + (upper-lower)*maximinESE_LHS(lhsDesign(n.init, d, seed=control$seed)$design)$design)
  }
  #----------------------------------------
  if (!is.null(X) && !is.null(y)) {
    obs.obj.init <- as.numeric(y)
    if (length(obs.obj.init) != n.init) {
      cat("Initial DOE (X) and objective (y) are not consistent.")
      return(0)
    }
  } else {
    obs.obj.init <- apply(design.init, 1, fun, ...)
  }
  #----------------------------------------
  if (!is.null(X) && !is.null(y)) {
    n.ite <- budget
  } else {
    n.ite <- budget - n.init
  }
  #----------------------------------------
  # Check control
  if (!is.null(control$cov.reestim)) cov.reestim <- control$cov.reestim else cov.reestim <- TRUE
  if (!is.null(control$model.trend)) model.trend <- control$model.trend else model.trend <- ~1
  if (!is.null(control$lb)) kmlb <- control$lb else kmlb <- (upper - lower)*sqrt(d)/100
  if (!is.null(control$ub)) kmub <- control$ub else kmub <- (upper - lower)*sqrt(d)
  if (!is.null(control$nugget)) nugget <- control$nugget else nugget <- NULL
  if (!is.null(control$covtype)) covtype <- control$covtype else covtype <- "matern5_2"
  if (!is.null(control$optim.method)) optim.method <- control$optim.method else optim.method <- "BFGS"
  if (!is.null(control$multistart)) multistart <- control$multistart else multistart <- 2*d
  if (is.null(control$gpmean.trick)) control$gpmean.trick <- FALSE
  if (is.null(control$gpmean.freq)) control$gpmean.freq <- 1e4
  if (!is.null(control$scaling)) scaling <- control$scaling else scaling <- FALSE
  if (is.null(control$warping)) control$warping <- FALSE
  if (!is.null(control$TR)) TR <- control$TR else TR <- FALSE
  if (!is.null(control$trcontrol)) trcontrol <- control$trcontrol else trcontrol <- NULL
  if (is.null(control$always.sample)) control$always.sample <- FALSE
  
  if (model.trend == "quad") {
    model.trend <- as.formula(paste0("y ~ .", paste0(paste0("+ I(", colnames(design.init)), "^2)", collapse = " "), collapse = " "))
  }
  
  if (control$warping) {
    o <- optim(par=runif(6), fn=logLw, lower=rep(c(rep(0,2), -2),each=2), upper=rep(5,6), y=obs.obj.init, method="L-BFGS-B")
    init.obs <- obs.obj.init
    wparam <- c(o$par, mean(obs.obj.init), sd(obs.obj.init))
    obs.obj.init <- warp(obs.obj.init, wparam=wparam)
    control$wparam <- wparam
  } else {
    control$warping <- TRUE
    init.obs <- obs.obj.init
    wparam <- c(mean(obs.obj.init), sd(obs.obj.init))
    obs.obj.init <- warp(obs.obj.init, wparam=wparam)
    control$wparam <- wparam
  }
  
  model <- km(formula=model.trend, design = design.init, response = obs.obj.init, covtype=covtype, 
                  optim.method=optim.method, multistart=multistart, nugget=nugget,
                  control=list(trace=FALSE), lower=kmlb, upper=kmub, scaling=scaling)
  
  #----------------------------------------
  optimcontrol = list(threshold = 1e-5, distance = "euclidean", notrace = !control$trace>0)
  
  if (TR) {
    res <- TREGO.nsteps(model = model, fun = fun, nsteps=n.ite, lower=lower, upper=upper, control=control, trace=control$trace, 
                    n.cores=n.cores, trcontrol=trcontrol, ...)
  } else {
    res <- fastEGO.nsteps(model = model, fun = fun, nsteps=n.ite, lower=lower, upper=upper, control=control, trace=control$trace, 
                    n.cores=n.cores, ...)
  }
  colnames(res$par) <- colnames(design.init)
  allX <- rbind(design.init, res$par)
  ally <- c(init.obs, res$value)

  par   <- allX[which.min(ally),]
  value <- min(ally)
  
  return(list(par=par, value = value, model=res$lastmodel, history=list(X=allX, y=ally), control=control, res=res))
}
