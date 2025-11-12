#' Sequential EI maximization and model re-estimation, with a number of
#' iterations fixed in advance by the user
#' 
#' Executes \emph{nsteps} iterations of the EGO method to an object of class
#' \code{\link[DiceKriging]{km}}.  At each step, a kriging model is
#' re-estimated (including covariance parameters re-estimation) based on the
#' initial design points plus the points visited during all previous
#' iterations; then a new point is obtained by maximizing the Expected
#' Improvement criterion (\code{\link{EI}}).
#' 
#' 
#' @param model an object of class \code{\link[DiceKriging]{km}} ,
#' @param fun the objective function to be minimized,
#' @param nsteps an integer representing the desired number of iterations,
#' @param lower vector of lower bounds for the variables to be optimized over,
#' @param upper vector of upper bounds for the variables to be optimized over,
#' @param control an optional list of control parameters for EGO. One
#' can control
#' 
#' \code{"warping"} whether or not a warping is applied to the outputs (default FALSE)
#' 
#' \code{"cov.reestim"} whether or not the covariance parameters are estimated at
#' each step (default TRUE)
#' \code{"gpmean.trick"} whether or not EI should be replaced periodically by the GP mean
#'  (default FALSE)
#'  
#' \code{"gpmean.freq"} frequency at which EI is replaced by the GP mean (default 1e4)
#' 
#' \code{"always.sample"} if TRUE, forces observation even if it creates poor conditioning
#' 
#' @param trace between -1 (no trace) and 3 (full messages)
#' @param n.cores number of cores used for EI maximisation
#' @param ... additional parameters to be given to \code{fun}
#' @return A list with components:
#' 
#' \item{par}{a data frame representing the additional points visited during
#' the algorithm,}
#' 
#' \item{value}{a data frame representing the response values at the points
#' given in \code{par},}
#' 
#' \item{npoints}{an integer representing the number of parallel computations
#' (=1 here),}
#' 
#' \item{nsteps}{an integer representing the desired number of iterations
#' (given in argument),}
#' 
#' \item{lastmodel}{an object of class \code{\link[DiceKriging]{km}}
#' corresponding to the last kriging model fitted. If warping is true, 
#' \code{y} values are normalized (warped) and will not match \code{value}.}
#' @author Victor Picheny
#' 
#' @seealso \code{\link{EI}}, \code{\link{max_crit}}, \code{\link{EI.grad}}
#' @references
#' 
#' D.R. Jones, M. Schonlau, and W.J. Welch (1998), Efficient global
#' optimization of expensive black-box functions, \emph{Journal of Global
#' Optimization}, 13, 455-492.
#' 
#' J. Mockus (1988), \emph{Bayesian Approach to Global Optimization}. Kluwer
#' academic publishers.
#' 
#' T.J. Santner, B.J. Williams, and W.J. Notz (2003), \emph{The design and
#' analysis of computer experiments}, Springer.
#' 
#' M. Schonlau (1997), \emph{Computer experiments and global optimization},
#' Ph.D. thesis, University of Waterloo.
#' @keywords optimize
#' @examples
#' 
#' set.seed(123)
#' ###############################################################
#' ### 	10 ITERATIONS OF EGO ON THE BRANIN FUNCTION, 	   ####
#' ###	 STARTING FROM A 9-POINTS FACTORIAL DESIGN         ####
#' ###############################################################
#' 
#' # a 9-points factorial design, and the corresponding response
#' d <- 2
#' n <- 9
#' design.fact <- expand.grid(seq(0,1,length=3), seq(0,1,length=3)) 
#' names(design.fact)<-c("x1", "x2")
#' design.fact <- data.frame(design.fact) 
#' names(design.fact)<-c("x1", "x2")
#' response.branin <- apply(design.fact, 1, branin)
#' response.branin <- data.frame(response.branin) 
#' names(response.branin) <- "y" 
#' 
#' # model identification
#' fitted.model1 <- km(~1, design=design.fact, response=response.branin, 
#' covtype="gauss", control=list(pop.size=50,trace=FALSE), parinit=c(0.5, 0.5))
#' 
#' # EGO n steps
#' nsteps <- 5 
#' lower <- rep(0,d) 
#' upper <- rep(1,d)     
#' oEGO <- fastEGO.nsteps(model=fitted.model1, fun=branin, nsteps=nsteps, lower=lower, upper=upper)
#' print(oEGO$par)
#' print(oEGO$value)
#' 
#' # graphics
#' n.grid <- 15 # Was 20, reduced to 15 for speeding up compilation
#' x.grid <- y.grid <- seq(0,1,length=n.grid)
#' design.grid <- expand.grid(x.grid, y.grid)
#' response.grid <- apply(design.grid, 1, branin)
#' z.grid <- matrix(response.grid, n.grid, n.grid)
#' contour(x.grid, y.grid, z.grid, 40)
#' title("Branin function")
#' points(design.fact[,1], design.fact[,2], pch=17, col="blue")
#' points(oEGO$par, pch=19, col="red")
#' text(oEGO$par[,1], oEGO$par[,2], labels=1:nsteps, pos=3)
#' 
#' @export fastEGO.nsteps
fastEGO.nsteps <-function(model, fun, nsteps, lower, upper, control=NULL, trace=0, n.cores=1, ...) {
  
  n <- nrow(model@X)
  n.proxy <- n.unif <- 0
  d <- model@d
  added.obs <- c()
  added.X <- c()
  
  if (is.null(control$warping)) control$warping <- FALSE
  if (is.null(control$cov.reestim)) control$cov.reestim <- TRUE
  if (!is.null(control$gpmean.freq)) gpmean.freq <- control$gpmean.freq else gpmean.freq <- 1e4
  if (!is.null(control$gpmean.trick)) gpmean.trick <- control$gpmean.trick else gpmean.trick <- FALSE
  if (!is.null(control$always.sample)) always.sample <- control$always.sample else always.sample <- FALSE
  
  observations <- model@y
  
  if (trace > 1) message("Ite | max(EI) ||  ", colnames(model@X),"|  obj \n")
  
  for (i in 1:nsteps) {
    
    if (((i %% gpmean.freq)==0) && (gpmean.trick)) {
      # Bypass EI maximization to switch to gpmean
      sol <- list(value=0)
    } else {
      sol <- try(max_crit(model=model, lower=lower, upper=upper, control=control, n.cores=n.cores))  
    }
    
    ## Exit if optimization failed
    if (typeof(sol) == "character") {
      warning("Unable to maximize EI at iteration ", i, "- optimization stopped \n Last model returned \n")
      
      par <- values <- c()
      if (i > 1) {
        par <- model@X[(n+1):model@n,, drop = FALSE]
        values <- model@y[(n+1):model@n]
      }
      return(list(par=par, values=values, nsteps = i-1, lastmodel = model))
    }
    
    if (sol$value < max(1e-16, sqrt(model@covariance@sd2/1e16)) && gpmean.trick) {
      # Restart optimization with a proxy criterion
      if (trace>0)
        if ((i %% gpmean.freq)==0) {
          message("Forcing proxy acquisition criterion.\n")
        } else {
          message("No point found with non-null EI, restarted with proxy criterion.\n")
          n.proxy <- n.proxy+1
        }                
      
      sol <- try(max_crit(model=model, lower=lower, upper=upper, control=control, proxy=TRUE, n.cores=n.cores))

      ## Exit if optimization failed
      if (typeof(sol) == "character") {
        warning("Unable to maximize proxy criterion - optimization stopped \n Last model returned \n")
        
        par <- values <- c()
        if (i > 1) {
          par <- model@X[(n+1):model@n,, drop = FALSE]
          values <- model@y[(n+1):model@n]
        }
        return(list(par=par, values=values, nsteps = i-1, lastmodel = model))
      }
    } else {
      n.proxy <- 0
    }
    
    ## Check if optimization does not return already known point
    replace.point = FALSE
    
    if (!always.sample) {
      if (checkPredict(sol$par, model= list(model), threshold = 1e-6)){
        if (trace>0)
          message("Optimization failed (an existing point was selected), so a random point is taken instead.\n")
        sol$par <- matrix(runif(d) * (upper - lower) + lower, nrow = 1) 
        n.unif <- n.unif + 1
      } else {
        n.unif <- 0
      }
    } else {
      if (checkPredict(sol$par, model= list(model), threshold = 1e-6)){
        if (trace>0)
          message("Point too close to an existing one - replacement scheme activated.\n")
        replace.point = TRUE
      }
    }
    
    ## Update
    X.new <- matrix(as.numeric(sol$par), nrow=1, ncol=d)
    Y.new <- try(fun(as.numeric(sol$par),...))

    if (typeof(Y.new) == "character") {
      if (trace>0) {
        warning("Unable to compute objective function at iteration ", i, "- optimization stopped \n")
        warning("Problem occured for the design: ", X.new, "\n Last model returned \n")
      }
      
      par <- values <- c()
      if (i > 1) {
        par <- added.X 
        values <- added.obs
      }
      return(list(par=par, values=values, nsteps = i-1, lastmodel = model))
    }
    if (trace > 1) message( i, "|", signif(sol$val,3), "||", signif(X.new,3), "|", signif(Y.new,3), "\n")
    
    # Transform Ynew if needed
    added.obs <- c(added.obs, Y.new)
    added.X <- rbind(added.X, X.new)
    if (control$warping) {
      Y.new <- warp(Y.new, control$wparam)
    }
    
    # Update models
    observations <- rbind(observations, Y.new)
    newmodel <- model
    
    if (!replace.point) {
      if (n.proxy < 5 && n.unif <3) {
        newmodel <- try(update(object = model, newX = X.new, newy=Y.new, newX.alreadyExist=FALSE,
                               cov.reestim = control$cov.reestim), silent = TRUE)
  
        if (typeof(newmodel) == "character" && control$cov.reestim) {
          newmodel <- try(update(object = model, newX = X.new, newy=Y.new, newX.alreadyExist=FALSE, cov.reestim = FALSE), silent = TRUE)
          if (typeof(newmodel) != "character") {
            if (trace > 0) message("Error in hyperparameter estimation - old hyperparameter values used instead for the objective model \n")
          } else {
            if (length(model@covariance@nugget) == 0) {
              newmodel <- try(km(formula=model@trend.formula, design=rbind(model@X, X.new), response=c(model@y, Y.new), 
                                 covtype=model@covariance@name, coef.trend=model@trend.coef, nugget= model@covariance@sd2/1e12,
                                 multistart=model@control$multistart, 
                                 coef.cov=covparam2vect(model@covariance), coef.var=model@covariance@sd2, iso=is(model@covariance,"covIso")))
              if (typeof(newmodel) != "character") {
                if (trace > 0) message("Error in hyperparameter estimation - nugget added to the model \n")
              }
            }
          }
        }
      } else {
        # If too many proxy or random iterations => update model with smaller range
        newmodel <- try(km(formula=model@trend.formula, design=rbind(model@X, X.new), response=c(model@y, Y.new), 
                           covtype=model@covariance@name, coef.trend=model@trend.coef, nugget= model@covariance@nugget,
                           multistart=model@control$multistart, 
                           coef.cov=covparam2vect(model@covariance)/1.5, coef.var=model@covariance@sd2, iso=is(model@covariance,"covIso")))
  
        if (typeof(newmodel) == "character" && length(model@covariance@nugget) == 0) {
          newmodel <- try(km(formula=model@trend.formula, design=rbind(model@X, X.new), response=c(model@y, Y.new), 
                             covtype=model@covariance@name, coef.trend=model@trend.coef, nugget= model@covariance@sd2/1e12,
                             multistart=model@control$multistart, 
                             coef.cov=covparam2vect(model@covariance)/1.5, coef.var=model@covariance@sd2, iso=is(model@covariance,"covIso")))
          if (typeof(newmodel) != "character") {
            if (trace > 0) message("Too many proxy/random iterations - model updated with smaller covariance ranges and nugget \n")
            model@lower <- model@lower/1.5
            model@upper <- model@upper/1.5
          }
        } else {
          if (trace > 0) message("Too many proxy/random iterations - model updated with smaller covariance ranges \n")
        }
      }
    } else {
      # Replace current point in the model by a new point if better
      if (Y.new < min(model@y)) {
        ibest <- which.min(model@y)
        model@y[ibest] <- Y.new
        model@X[ibest,] <- X.new
        newmodel <- computeAuxVariables(model)
      }
    }
    
    if (typeof(newmodel) == "character") {
      warning("Unable to udpate kriging model at iteration", i-1, "- optimization stopped \n")
      warning("lastmodel is the model at iteration", i-1, "\n")
      warning("par and values contain the ",i, "th observation \n \n")
      if (i > 1) allX.new <- rbind(model@X[(n+1):model@n,, drop=FALSE], X.new)
      return(list(
        par    = added.X,
        values = added.obs,
        nsteps = i, 
        lastmodel = model))
    } else {
      newmodel@known.param <- model@known.param
      newmodel@case <- model@case
      newmodel@param.estim <- model@param.estim
      newmodel@method <- model@method
      newmodel@penalty <- model@penalty
      newmodel@optim.method <- model@optim.method
      newmodel@lower <- model@lower
      newmodel@upper <- model@upper
      newmodel@control <- model@control
      
      model <- newmodel
    }
    if (trace > 2) message("Covariance parameters: ", model@covariance@sd2, "|", covparam2vect(model@covariance), "\n")
  }
  
  if(trace>0) message("\n")
  #### End of main loop ################
  
  return(list(par=added.X, value=added.obs, npoints=1, nsteps=nsteps, lastmodel=model))
}