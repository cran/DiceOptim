#' Trust-region based EGO algorithm. 
#' 
#' Executes \emph{nsteps} iterations of the TREGO method to an object of class
#' \code{\link[DiceKriging]{km}}.  At each step, a kriging model is
#' re-estimated (including covariance parameters re-estimation) based on the
#' initial design points plus the points visited during all previous
#' iterations; then a new point is obtained by maximizing the Expected
#' Improvement criterion (\code{\link{EI}}) over either the entire search space
#' or restricted to a trust region. The trust region is updated at each iteration
#' based on a sufficient decrease condition.
#' 
#' @param model an object of class \code{\link[DiceKriging]{km}},
#' @param fun the objective function to be minimized,
#' @param nsteps an integer representing the desired number of iterations,
#' @param lower,upper vector of lower and upper bounds for the variables to be optimized over,
#' @param control an optional list of control parameters for optimization.  
#' For now only the number of \code{restarts} can be set.
#' @param kmcontrol an optional list representing the control variables for
#' the re-estimation of the kriging model.
#' @param trcontrol an optional list of control parameters for the trust-region scheme: 
#' \code{sigma} the initial size of the trust region,
#' \code{x0} its initial center, 
#' \code{beta} the contraction factor,
#' \code{alpha} its dilatation factor, 
#' \code{kappa} the forcing factor, 
#' \code{crit} the criterion used inside the TR (either "EI" or "gpmean"),
#' \code{GLratio} number of consecutive global and local steps,
#' \code{algo} either "TREGO" or "TRIKE",
#' \code{minsigma} minimal sigma value,
#' \code{maxsigma} maximal sigma value,
#' \code{minEI} stopping criterion for TRIKE,
#' \code{local.model} Boolean; if TRUE, a local model is used within the trust region,
#' \code{local.trend, local.covtype} trend and covariance for the local model,
#' \code{n.local.min} minimal number of points used to build the local model,
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
#' 
#' \item{all.success}{ a vector of Boolean indicating the successful steps according
#' to the sufficient decrease condtion}
#' 
#' \item{all.steps}{ a vector of Boolean indicating which steps were global}
#' 
#' \item{all.sigma}{ history of trust region size}
#' 
#' \item{all.x0}{ history of trust region centers}
#' 
#' \item{local.model}{ if trcontrol$local.model=TRUE, the latest local model}
#' @author Victor Picheny
#' 
#' @seealso \code{\link{EI}}, \code{\link{max_crit}}, \code{\link{EI.grad}}
#' @references
#' 
#' Diouane, Picheny, Le Riche, Scotto Di Perrotolo (2021),
#' \emph{TREGO: a Trust-Region Framework for Efficient Global Optimization}, ArXiv 
#' 
#' @keywords optimize
#' @importFrom stats as.formula var runif
#' @examples
#' 
#' 
#' set.seed(123)
#' ###############################################################
#' ### 	10 ITERATIONS OF TREGO ON THE BRANIN FUNCTION, 	 ####
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
#' # TREGO n steps
#' nsteps <- 5
#' lower <- rep(0, d) 
#' upper <- rep(1, d)     
#' oEGO <- TREGO.nsteps(model=fitted.model1, fun=branin, nsteps=nsteps, 
#' lower=lower, upper=upper)
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
#' @export TREGO.nsteps
TREGO.nsteps <- function(model, fun, nsteps, lower, upper, control=NULL, kmcontrol=NULL, trcontrol=NULL, trace=0, n.cores=1, ...) {
  
  n <- nrow(model@X)
  d <- model@d
  added.obs <- c()
  
  # Check control
  if (is.null(control$warping)) control$warping <- FALSE
  
  # Check kmcontrol
  if (is.null(kmcontrol$penalty)) kmcontrol$penalty <- model@penalty
  if (length(model@penalty==0)) kmcontrol$penalty <- NULL 
  if (is.null(kmcontrol$optim.method)) kmcontrol$optim.method <- model@optim.method 
  if (is.null(kmcontrol$parinit)) kmcontrol$parinit <- model@parinit
  if (is.null(kmcontrol$control)) kmcontrol$control <- model@control
  if (is.null(kmcontrol$cov.reestim)) kmcontrol$cov.reestim <- TRUE
  
  # Check trcontrol
  if (is.null(trcontrol$sigma)) trcontrol$sigma <- 0.5*(upper - lower) / (5^(1/model@d))
  if (is.null(trcontrol$x0)) trcontrol$x0 <- model@X[which.min(model@y),]
  if (is.null(trcontrol$beta)) trcontrol$beta <- 0.9
  if (is.null(trcontrol$alpha)) trcontrol$alpha <-  1/trcontrol$beta
  if (is.null(trcontrol$kappa)) trcontrol$kappa <- 1e-4
  if (is.null(trcontrol$crit)) trcontrol$crit <- "EI"
  if (is.null(trcontrol$GLratio)) trcontrol$GLratio <- c(1,1)
  if (is.null(trcontrol$algo))   trcontrol$algo <- "TREGO"
  if (is.null(trcontrol$minsigma)) trcontrol$minsigma <-1/64*min(upper - lower)
  if (is.null(trcontrol$maxsigma)) trcontrol$maxsigma <- 2*min(upper - lower)
  if (is.null(trcontrol$minEI)) trcontrol$minEI <- 1e-6
  if (is.null(trcontrol$local.model)) trcontrol$local.model <- FALSE
  if (is.null(trcontrol$local.trend)) trcontrol$local.trend <- ~1
  if (trcontrol$local.trend == "quad") {
    trcontrol$local.trend <- as.formula(paste0("y ~ .", paste0(paste0("+ I(", colnames(model@X)), "^2)", collapse = " "), collapse = " "))
  }
  if (is.null(trcontrol$local.covtype)) trcontrol$local.covtype <- model@covariance@name
  
  if (is.null(trcontrol$n.local.min)) {
    if (length(as.formula(trcontrol$local.trend)) == 2) {
      if (trcontrol$local.trend == ~1) {
        trcontrol$n.local.min <- d+5
      } else if (trcontrol$local.trend == ~.) {
        trcontrol$n.local.min <- 2*d+5
      }
    } else {
      trcontrol$n.local.min <- 2*d+1
    }
  }
  
  max.global <- trcontrol$GLratio[1]
  max.local <- trcontrol$GLratio[2]
  n.global <- 0
  n.local <- 0
  
  observations <- model@y
  
  if (trace > 1) message("It |  Step  | Success | Crit  | maxcrit ||  ", colnames(model@X),"||  obj  \n")
  
  if (max.global == 0) {
    global.step <- FALSE
  } else if (max.local == 0){
    global.step <- TRUE
  } else {
    global.step <- TRUE
  }
  
  ### MAIN LOOP STARTS HERE ###
  all.steps <- all.success <- all.sigma <- all.x0 <- c()
  
  for (i in 1:nsteps) {
    
    all.steps <- c(all.steps, global.step)
    all.sigma <- c(all.sigma, mean(trcontrol$sigma))
    all.x0 <- rbind(all.x0, trcontrol$x0)
    if (global.step){
      trcontrol$inTR <- FALSE 
      optim.model <- model
    } else {
      trcontrol$inTR <- TRUE
      
      # Build local model if needed
      if (trcontrol$local.model) {
        TRlower <- trcontrol$x0 - trcontrol$sigma
        TRupper <- trcontrol$x0 + trcontrol$sigma
        I1 <- which( (rowSums(model@X > TRlower) + rowSums(model@X < TRupper)) == 2*model@d)
        
        if (length(I1) < trcontrol$n.local.min) {
          tp1 <- as.numeric(t(model@X))
          tp2 <- matrix(tp1 - as.numeric(trcontrol$x0), nrow = model@n, byrow = TRUE)^2
          I2 <- order(rowSums(tp2))
          I1 <- unique(c(I1, I2))[1:trcontrol$n.local.min]
        }
        
        X.local <- model@X[I1,]
        y.local <- model@y[I1]
        
        local.model <- try(km(formula=trcontrol$local.trend, design=X.local, response=y.local, 
                              covtype=trcontrol$local.covtype, nugget= trcontrol$local.nugget,
                              multistart=model@control$multistart, control=model@control))
        
        if (typeof(local.model) == "character") {
          if (is.null(trcontrol$local.nugget)) {
            local.model <- try(km(formula=trcontrol$local.trend, design=X.local, response=y.local, 
                                  covtype=trcontrol$local.covtype, nugget= var(y.local)*1e-12,
                                  multistart=model@control$multistart, control=model@control))
            
            if (typeof(local.model) == "character") {
              if (trace > 0) message("Error in local model hyperparameter estimation - nugget added to the model \n")
            }
          }
        }
        if (typeof(local.model) == "character") {
          warning("Unable to update local model at iteration", i, "- optimization stopped \n")
          return(list(par=model@X[(n+1):(n+i-1),, drop=FALSE], 
                      value=added.obs, npoints=1, nsteps=nsteps, lastmodel=model, local.model=optim.model,
                      all.success=all.success, all.steps=all.steps, all.sigma=all.sigma, all.x0=all.x0))
        } else {
          optim.model <- local.model
        }
      } else {
        optim.model <- model
      }
    }
    
    ## Change the seeds for genoud to avoid selecting always the same initial values
    if(is.null(control$unif.seed))      control$unif.seed <- runif(1)
    
    sol <- try(max_crit(model=optim.model, lower=lower, upper=upper, control=control,
                        trcontrol=trcontrol, n.cores=n.cores))
    
    crit.used <- "  EI  "
    
    ## Exit if optimization failed
    if (typeof(sol) == "character") {
      warning("Unable to maximize EI at iteration ", i, "- optimization stopped \n Last model returned \n")
      
      par <- values <- c()
      if (i > 1) {
        par <- model@X[(n+1):(n+i-1),, drop = FALSE]
        values <- model@y[(n+1):(n+i-1)] 
      }
      return(list(par=par, 
                  value=values, npoints=1, nsteps=i-1, lastmodel=model, local.model=optim.model,
                  all.success=all.success, all.steps=all.steps, all.sigma=all.sigma, all.x0=all.x0))
    }
    
    # Early stop for TRIKE
    if (sol$value < trcontrol$minEI && trcontrol$algo == "TRIKE") {
      if (trace>0) message("TRIKE terminated due to too small EI \n")
      
      par <- values <- c()
      if (i > 1) {
        par <- model@X[(n+1):(n+i-1),, drop = FALSE]
        values <- model@y[(n+1):(n+i-1)]
      }
      return(list(par=par, 
                  value=values, npoints=1, nsteps=i-1, lastmodel=model, local.model=optim.model,
                  all.success=all.success, all.steps=all.steps, all.sigma=all.sigma, all.x0=all.x0))
    }
    
    if (sol$value < max(1e-16, optim.model@covariance@sd2/1e12) && (!trcontrol$inTR || (trcontrol$inTR && trcontrol$crit=="EI") )  ) {
      # Restart optimization with a proxy criterion
      sol <- try(max_crit(model=optim.model, lower=lower, upper=upper, control=control, proxy=TRUE,
                          trcontrol=trcontrol, n.cores=n.cores))
      
      crit.used <- "gpmean"
      
      ## Exit if optimization failed
      if (typeof(sol) == "character") {
        warning("Unable to maximize proxy criterion - optimization stopped \n Last model returned \n")
        
        par <- values <- c()
        if (i > 1) {
          par <- model@X[(n+1):(n+i-1),, drop = FALSE]
          values <- model@y[(n+1):(n+i-1)]
        }
        return(list(par=par, 
                    value=values, npoints=1, nsteps=i-1, lastmodel=model, local.model=optim.model,
                    all.success=all.success, all.steps=all.steps, all.sigma=all.sigma, all.x0=all.x0))
      }
    }
    
    ## Check if optimization do not return already known point
    if(checkPredict(sol$par, model= list(optim.model), threshold = 1e-6)){
      sol$par <- matrix(runif(d), nrow = 1)
      crit.used <- "random"
    }
    
    ## Update
    X.new <- matrix(as.numeric(sol$par), nrow=1, ncol=d)
    Y.new <- try(fun(as.numeric(sol$par),...))
    
    if (typeof(Y.new) == "character") {
      warning("Unable to compute objective function at iteration ", i, "- optimization stopped \n")
      warning("Problem occured for the design: ", X.new, "\n Last model returned \n")
      
      par <- values <- c()
      if (i > 1) {
        par <- model@X[(n+1):model@n,, drop=FALSE]
        values <- added.obs
      }
      return(list(par=par, 
                  value=values, npoints=1, nsteps=i-1, lastmodel=model, local.model=optim.model,
                  all.success=all.success, all.steps=all.steps, all.sigma=all.sigma, all.x0=all.x0))
    }
    
    # Transform Ynew if needed
    added.obs <- c(added.obs, Y.new)
    Y.new.or <- Y.new
    if (control$warping) {
      Y.new <- warp(Y.new, control$wparam)
    }
    
    
    observations <- rbind(observations, Y.new)
    newmodel <- model
    
    if (trcontrol$algo == "TRIKE") {
      AI <- max(0, min(observations[-length(observations)]) - observations[length(observations)])
      EI <- sol$value
      success <- (AI / EI) > 1
    } else {
      success <- observations[length(observations)] < (min(observations[-length(observations)]) - trcontrol$kappa*mean(trcontrol$sigma^2))
    }
    if (trace > 1) message( i, "|", switch(2-global.step, "global", "local "),"|", switch(2-success, "success", "failure"), "|", crit.used,  "|", signif(sol$val,3), "||", 
                        signif(X.new,3), "||", signif(Y.new.or,3), "\n")
    
    # Update number of consecutive trial / local steps
    if (global.step) n.global <- n.global + 1 else n.local <- n.local + 1
    
    all.success <- c(all.success, success)
    
    # Update models
    a1 <- Sys.time()
    
    newmodel <- try(update(object = model, newX = X.new, newy=Y.new, newX.alreadyExist=FALSE,
                           cov.reestim = kmcontrol$cov.reestim), silent = TRUE)
    
    if (typeof(newmodel) == "character" && kmcontrol$cov.reestim) {
      newmodel <- try(update(object = model, newX = X.new, newy=Y.new, newX.alreadyExist=FALSE, cov.reestim = FALSE), silent = TRUE)
      if (typeof(newmodel) != "character") {
        if (trace > 0) message("Error in hyperparameter estimation - old hyperparameter values used instead for the objective model \n")
      } else {
        if (length(model@covariance@nugget) == 0) {
          newmodel <- try(km(formula=model@trend.formula, design=rbind(model@X, X.new), response=c(model@y, Y.new), 
                             covtype=model@covariance@name, coef.trend=model@trend.coef, nugget= model@covariance@sd2/1e12,
                             coef.cov=covparam2vect(model@covariance), coef.var=model@covariance@sd2, iso=is(model@covariance,"covIso")))
          if (typeof(newmodel) != "character") {
            if (trace > 0) message("Error in hyperparameter estimation - nugget added to the model \n")
          }
        }
      }
      newmodel@known.param <- model@known.param
      newmodel@case <- model@case
      newmodel@param.estim <- model@param.estim
      newmodel@method <- model@method
      newmodel@penalty <- model@penalty
      newmodel@optim.method <- model@optim.method
      newmodel@lower <- model@lower
      newmodel@upper <- model@upper
      newmodel@control <- model@control
    }
    
    if (typeof(newmodel) == "character") {
      warning("Unable to update kriging model at iteration", i-1, "- optimization stopped \n")
      warning("lastmodel is the model at iteration", i-1, "\n")
      warning("par and values contain the ",i, "th observation \n \n")
      if (i > 1) allX.new <- rbind(model@X[(n+1):(n+i-1),, drop=FALSE], X.new)
      return(list(
        par    = allX.new,
        values = added.obs,
        nsteps = i, 
        npoints=1, 
        nsteps=i-1, 
        lastmodel=model, local.model=optim.model,
        all.success=all.success, all.steps=all.steps, all.sigma=all.sigma, all.x0=all.x0))
    } else {
      model <- newmodel
    }
    if (trace > 2) message("Covariance parameters: ", model@covariance@sd2, "|", covparam2vect(model@covariance), "\n")
    if (trace > 2) message(Sys.time() - a1)
    
    # Update Trust-region scheme
    if (trcontrol$algo == "TRIKE") {
      trcontrol$x0 <- model@X[which.min(model@y),]
      if (success) {
        trcontrol$sigma <- pmin(trcontrol$alpha*trcontrol$sigma, trcontrol$maxsigma)
      } else {
        if (AI == 0){
          trcontrol$sigma <- pmax(trcontrol$beta*trcontrol$sigma, trcontrol$minsigma)
        }
      }
      global.step <- FALSE
    } else {
      if (success) {
        ### Successfull step: update TR and go back to EGO step if needed
        trcontrol$x0 <- model@X[which.min(model@y),]
        
        if (global.step) {
          # If successful EGO step: keep global and reinitialize counters
          global.step <- TRUE
          n.global <- n.local <- 0
        } else {
          # If successful local step: switch to global only if local steps exhausted
          trcontrol$sigma <- trcontrol$alpha*trcontrol$sigma
          if (n.local >= max.local) {
            global.step <- TRUE
            n.global <- n.local <- 0
          }
        }
      } else {
        ### Unsuccessfull step
        if (global.step) {
          # Unsuccessfull EGO step: go to TR step if max.global attained
          if (n.global >= max.global) {
            global.step <- FALSE
            n.global <- n.local <- 0
          }
        } else {
          # Unsuccessfull TR step: update TR and go back to EGO if max.local attained
          trcontrol$sigma <- trcontrol$beta*trcontrol$sigma
          if (n.local >= max.local) {
            global.step <- TRUE
            n.global <- n.local <- 0
          }
        }
      }
      if (max.global == 0) {
        global.step <- FALSE
      } else if (max.local == 0){
        global.step <- TRUE
      }
    } # Closes the update if TR
  } # Closes the main loop
  
  if(trace>0) message("\n")
  #### End of main loop ################
  if (trcontrol$local.model) {
    return(list(par=model@X[(n+1):(n+nsteps),, drop=FALSE], 
                value=added.obs, npoints=1, nsteps=nsteps, lastmodel=model, local.model=local.model,
                all.success=all.success, all.steps=all.steps, all.sigma=all.sigma, all.x0=all.x0))
  } else {
    return(list(par=model@X[(n+1):(n+nsteps),, drop=FALSE], 
                value=added.obs, npoints=1, nsteps=nsteps, lastmodel=model, 
                all.success=all.success, all.steps=all.steps, all.sigma=all.sigma, all.x0=all.x0))
  }
}
