##' Executes \code{nsteps} iterations of EGO methods integrating constraints, based on objects of class \code{\link[DiceKriging]{km}}.
##' At each step, kriging models are re-estimated (including covariance parameters re-estimation)
##'  based on the initial design points plus the points visited during all previous iterations;
##'  then a new point is obtained by maximizing one of the constrained Expected Improvement criteria available. 
##' @title Sequential constrained Expected Improvement maximization and model re-estimation,
##'  with a number of iterations fixed in advance by the user
##' @details Extension of the function \code{\link[DiceOptim]{EGO.nsteps}} to constrained optimization.\cr
##' 
##' The problem considered is of the form: \eqn{min f(x)} s.t. \eqn{g(x) \le 0}, 
##' \eqn{g} having a vectorial output. 
##' By default all its components are supposed to be inequalities, but one can use a boolean vector in \code{equality} to specify which are equality constraints.
##' In this case one can modify the tolerance on the constraints using the \code{tolConstraints} component of \code{critcontrol}:
##' an optional vector giving a tolerance for each of the constraints (equality or inequality). 
##' It is highly recommended to use it when there are equality constraints since the default tolerance of \code{0.05} in such case might not be suited.\cr
##' 
##' Available infill criteria with \code{crit} are: \cr
##' \itemize{
##' \item Expected Probability of Feasibily (\code{EFI}) \code{\link[DOlab]{crit_EFI}},
##' \item Augmented Lagrangian (\code{AL}) \code{\link[DOlab]{crit_AL}},
##' \item Stepwise Uncertainty Reduction of the excursion volume (\code{SUR}) \code{\link[DOlab]{crit_SUR_cst}}.
##' }
##' Depending on the selected criterion, various parameters are available.
##' More precisions are given in the corresponding help pages.\cr 
##' 
##' It is possible to consider a cheap to evaluate objective function submitted to expensive constraints. 
##' In this case, provide only a function in \code{cheapfun}, with both \code{model.fun} and \code{fun} to NULL, see examples below.
##' 
##' @param model.fun object of class \code{\link[DiceKriging]{km}} corresponding to the objective function,
##' @param fun scalar function to be minimized, corresponding to \code{model.fun} found by a call to \code{\link[base]{match.fun}},
##' @param cheapfun optional scalar function to use if the objective is a fast-to-evaluate function (handled next with class \code{\link[DOlab]{fastfun}},
##' through the use of \code{\link[base]{match.fun}}), 
##' which does not need a kriging model, see details below,
##' @param constraint vectorial function corresponding to the constraints, see details below,
##' @param equality either \code{FALSE} if all constraints are for inequalities, else a vector of boolean indicating which are equalities 
##' @param model.constraint either one or a list of models of class \code{\link[DiceKriging]{km}}, one per constraint, 
##' @param crit choice of constrained improvement function: "\code{AL}", "\code{EFI}" or "\code{SUR}",
##' see details below,
##' @param nsteps an integer representing the desired number of iterations,
##' @param lower vector of lower bounds for the variables to be optimized over,
##' @param upper vector of upper bounds for the variables to be optimized over,
##' @param type "\code{SK}" or "\code{UK}" (by default), depending whether uncertainty related to trend estimation has to be taken into account,
##' @param cov.reestim optional boolean specifying if the kriging hyperparameters should be re-estimated at each iteration,
##' @param critcontrol optional list of parameters for criterion \code{crit}, see details,
##' @param optimcontrol an optional list of control parameters for optimization of the selected infill criterion:
##' \itemize{
##' \item{\code{method} can be set to "\code{discrete}" or "\code{genoud}". For "\code{discrete}", a matrix \code{candidate.points} must be given,
##' For "\code{genoud}", specific parameters to the chosen method can also be specified  (see \code{\link[rgenoud]{genoud}}).}
##' \item{Options for the \code{\link[DOlab]{checkPredict}} function: \code{threshold} (\code{1e-4}) and \code{distance} (\code{covdist}) are used to avoid 
##' numerical issues occuring when adding points too close to the existing ones.}
##' \item{\code{notrace} can be set to \code{TRUE} to suppress printing of the optimization progresses.}
##' }
##'  
##' 
##' @param ... additional parameters to be given to the objective \code{fun} and \code{constraint}.
##' @export
##' @import DiceKriging
##' @return
##' A list with components:
##' \itemize{
##' \item{\code{par}}{: a matrix representing the additional points visited during the algorithm,}
##' \item{\code{values}}{: a vector representing the response (objective) values at the points given in \code{par},}
##' \item{\code{constraint}}{: a matrix representing the constraints values at the points given in \code{par},}
##' \item{\code{feasibility}}{: a boolean vector saying if points given in \code{par} respect the constraints,}
##' \item{\code{nsteps}}{: an integer representing the desired number of iterations (given in argument),}
##' \item{\code{lastmodel.fun}}{: an object of class \code{\link[DiceKriging]{km}} corresponding to the objective function,}
##' \item{\code{lastmodel.constraint}}{: one or a list of objects of class \code{\link[DiceKriging]{km}} corresponding to the last kriging models fitted to the constraints.}\cr \cr
##' If a problem occurs during either model updates or criterion maximization, the last working model and corresponding values are returned.
##' }
##' 
##' @author
##' Victor Picheny 
##' 
##' Mickael Binois 
##' 
##' @seealso \code{\link{critcst_optimizer}}, \code{\link{crit_EFI}}, \code{\link{crit_AL}},
##' \code{\link{crit_SUR_cst}}, \code{\link{easyEGO.cst}}
##' 
##' @references
##' D.R. Jones, M. Schonlau, and W.J. Welch (1998), Efficient global
##' optimization of expensive black-box functions, \emph{Journal of Global
##' Optimization}, 13, 455-492.
##' 
##' 
##' M. Schonlau, W.J. Welch, and D.R. Jones (1998),
##'  Global versus local search in constrained optimization of computer models,
##'  \emph{Lecture Notes-Monograph Series}, 11-25.
##'  
##' M.J. Sasena, P. Papalambros, and P.Goovaerts (2002),
##'  Exploration of metamodeling sampling criteria for constrained global optimization,
##'  \emph{Engineering optimization}, 34, 263-278.
##'  
##' R.B. Gramacy, G.A. Gray, S. Le Digabel, H.K.H Lee, P. Ranjan, G. Wells, Garth, and S.M. Wild (2014+),
##' Modeling an augmented Lagrangian for improved blackbox constrained optimization,
##' \emph{arXiv preprint arXiv:1403.4890}.
##' 
##' J.M. Parr (2012),
##' \emph{Improvement criteria for constraint handling and multiobjective optimization},
##' University of Southampton.
##' 
##' V. Picheny (2014),
##' A stepwise uncertainty reduction approach to constrained global optimization,
##' \emph{Proceedings of the 17th International Conference on Artificial Intelligence and Statistics},  JMLR W&CP 33, 787-795.
##' 
##' @examples
##' #---------------------------------------------------------------------------
##' # 2D objective function, 3 cases
##' #---------------------------------------------------------------------------
##' \donttest{
##' set.seed(25468)
##' library(DiceDesign)
##' 
##' n_var <- 2 
##' fun <- goldsteinprice
##' fun1.cst <- function(x){return(-branin(x) + 25)}
##' fun2.cst <- function(x){return(3/2 - x[1] - 2*x[2] - .5*sin(2*pi*(x[1]^2 - 2*x[2])))}
##' 
##' # Constraint function with vectorial output
##' cstfun <- function(x){
##'   return(c(fun1.cst(x), fun2.cst(x)))
##' }
##' 
##' # For illustration purposes
##' n.grid <- 31
##' test.grid <- expand.grid(X1 = seq(0, 1, length.out = n.grid), X2 = seq(0, 1, length.out = n.grid))
##' obj.grid <- apply(test.grid, 1, fun)
##' cst1.grid <- apply(test.grid, 1, fun1.cst)
##' cst2.grid <- apply(test.grid, 1, fun2.cst)
##' 
##' # Initial set of observations and models
##' n.init <- 12 
##' design.grid <- round(maximinESE_LHS(lhsDesign(n.init, n_var, seed = 42)$design)$design, 1)
##' obj.init <- apply(design.grid, 1, fun)
##' cst1.init <- apply(design.grid, 1, fun1.cst)
##' cst2.init <- apply(design.grid, 1, fun2.cst)
##' model.fun <- km(~., design = design.grid, response = obj.init)
##' model.constraint1 <- km(~., design = design.grid, response = cst1.init, lower=c(.2,.2))
##' model.constraint2 <- km(~., design = design.grid, response = cst2.init, lower=c(.2,.2))
##' model.constraint <- list(model.constraint1, model.constraint2)
##' 
##' lower <- rep(0, n_var)
##' upper <- rep(1, n_var)
##'
##' #---------------------------------------------------------------------------
##' # 1- Expected Feasible Improvement criterion, expensive objective function,
##' # two inequality constraints, 5 iterations, using genoud
##' #---------------------------------------------------------------------------
##' 
##' cstEGO <- EGO.cst(model.fun = model.fun, fun = fun, model.constraint = model.constraint,
##'                   crit = "EFI", constraint = cstfun, equality = FALSE, lower = lower, 
##'                   upper = upper, nsteps = 5, optimcontrol = list(method = "genoud", maxit = 20))
##' 
##' # Plots: objective function in colour, constraint boundaries in red
##' # Initial DoE: white circles, added points: blue crosses, best solution: red cross
##' 
##' filled.contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), nlevels = 50,
##'                matrix(obj.grid, n.grid), main = "Two inequality constraints",
##'                xlab = expression(x[1]), ylab = expression(x[2]), color = terrain.colors, 
##'                plot.axes = {axis(1); axis(2);
##'                             points(design.grid[,1], design.grid[,2], pch = 21, bg = "white")
##'                             contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), 
##'                                     matrix(cst1.grid, n.grid), level = 0, add=TRUE,drawlabels=FALSE,
##'                                     lwd=1.5, col = "red")
##'                             contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), 
##'                                     matrix(cst2.grid, n.grid), level = 0, add=TRUE,drawlabels=FALSE,
##'                                     lwd=1.5, col = "red")
##'                             points(cstEGO$par, col = "blue", pch = 4, lwd = 2)
##'                }
##' )
##' 
##' #---------------------------------------------------------------------------
##' # 2- Augmented Lagrangian Improvement criterion, expensive objective function,
##' # one inequality and one equality constraint, using a discrete set of candidates (grid)
##' #---------------------------------------------------------------------------
##' cstEGO2 <- EGO.cst(model.fun = model.fun, fun = fun, model.constraint = model.constraint,
##'                    crit = "AL", constraint = cstfun, equality = c(TRUE, FALSE), lower = lower,  
##'                    upper = upper, nsteps = 10,
##'                    critcontrol = list(tolConstraints = c(2, 0), always.update=TRUE),
##'                    optimcontrol=list(method="discrete", candidate.points=test.grid))
##' 
##' # Plots: objective function in colour, inequality constraint boundary in red,
##' # equality constraint in orange
##' # Initial DoE: white circles, added points: blue crosses, best solution: red cross
##' 
##' filled.contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), nlevels = 50,
##'                matrix(obj.grid, n.grid),
##'                main = "Inequality (red) and equality (orange) constraints",
##'                xlab = expression(x[1]), ylab = expression(x[2]), color = terrain.colors, 
##'                plot.axes = {axis(1); axis(2);
##'                             points(design.grid[,1], design.grid[,2], pch = 21, bg = "white")
##'                             contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), 
##'                                     matrix(cst1.grid, n.grid), level = 0, add=TRUE,
##'                                     drawlabels=FALSE,lwd=1.5, col = "orange")
##'                             contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid),
##'                                     matrix(cst2.grid, n.grid), level = 0, add=TRUE,
##'                                     drawlabels=FALSE,lwd=1.5, col = "red")
##'                             points(cstEGO2$par, col = "blue", pch = 4, lwd = 2)
##'                }
##' )
##' 
##' #---------------------------------------------------------------------------
##' # 3- Stepwise Uncertainty Reduction criterion, fast objective function,
##' # single inequality constraint, 5 steps, importance sampling scheme
##' #---------------------------------------------------------------------------
##' 
##' cstEGO3 <- EGO.cst(model.fun = NULL, fun = NULL, cheapfun = fun,
##'                    model.constraint = model.constraint2, constraint = fun2.cst,
##'                    crit = "SUR", lower = lower, upper = upper,
##'                    nsteps =5, critcontrol=list(distrib="SUR"))
##' 
##' # Plots: objective function in colour, inequality constraint boundary in red,
##' # Initial DoE: white circles, added points: blue crosses, best solution: red cross
##' 
##' filled.contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), nlevels = 50,
##'                matrix(obj.grid, n.grid), main = "Single constraint, fast objective",
##'                xlab = expression(x[1]), ylab = expression(x[2]), color = terrain.colors, 
##'                plot.axes = {axis(1); axis(2);
##'                             points(design.grid[,1], design.grid[,2], pch = 21, bg = "white")
##'                             contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), 
##'                                     matrix(obj.grid, n.grid), nlevels = 10, add = TRUE,
##'                                     drawlabels = TRUE)
##'                             contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), 
##'                                     matrix(cst2.grid, n.grid), level = 0, add=TRUE,
##'                                     drawlabels=FALSE,lwd=1.5, col = "black")
##'                             points(cstEGO3$par, col = "blue", pch = 4, lwd = 2)
##'                                            }
##'                 )
##' }
##' 

EGO.cst <- function (model.fun=NULL, fun, cheapfun = NULL, model.constraint, constraint, equality = FALSE, crit="EFI", nsteps, lower, upper, type="UK", cov.reestim=TRUE,
                     critcontrol = NULL, optimcontrol = list(method="genoud", threshold = 1e-5, distance = "euclidean", notrace = FALSE), ...){
  ##########################################################################################
  # Inputs :
  # model.fun : km or fastfun object
  # fun: scalar objective function
  # constraint: vectorial function (vectorial output)
  # nsteps: number of iterations
  # lower, upper: design region
  # optimcontrol, type, CovReEstimate: parameters as in DiceOptim & KrigInv
  ##########################################################################################

  if(!is.null(model.fun) && !is.null(cheapfun)){
    cat("Either model.fun or cheapfun should be set to null.")
    return(0)
  }
  
  n.cst <- length(model.constraint)
  
  if (n.cst!=length(equality)) equality <- rep(FALSE, n.cst)
  n.ineq <- sum(!equality)
  
  if(crit == "SUR" && n.cst > 3){
    cat("crit_SUR_cst does not take more than 3 constraints \n")
    return(NA)
  }
  
  if(n.cst == 1 && class(model.constraint) != "list") model.constraint <- list(model.constraint)
  
  # Build fastfun if necessary  
  if (!is.null(cheapfun)){
    if(!is.null(model.fun) && !is.null(fun)){
      cat("model.fun and fun should be set to NULL for fast-to-evaluate objective mode \n")
      return(NA)
    }else{
      cat("Fastfun-Mode on \n")
      cheapfun <- match.fun(cheapfun)
      fastobs <- apply(model.constraint[[1]]@X, 1, cheapfun)
      model.fun <- fastfun(fn = cheapfun, design = model.constraint[[1]]@X, response = fastobs)
      fun <- match.fun(cheapfun)
    }
  } else {
    fun <- match.fun(fun)
  }
  
  n     <- nrow(model.fun@X)
  d     <- model.fun@d
  
  # Store objective values of observations
  observations_obj <- model.fun@y
  
  # Regroup all constraints values from observations
  observations_cst <- c()
  for (i in 1:n.cst) observations_cst <- cbind(observations_cst, model.constraint[[i]]@y)
  
  # Feasibility of the added observations
#   feasibility <- rep(TRUE, n)
#   feasibility[which(observations_cst > 0, arr.ind = T)[,1]] <- FALSE
  feasibility <- test_feas_vec(cst=observations_cst, equality=equality, tolConstraints = critcontrol$tolConstraints)
  
  if(is.null(optimcontrol$notrace)) notrace <- FALSE
  else                              notrace <- optimcontrol$notrace
  
  # Initialize critcontrol
  if (crit=="AL") {
    ## smallest sum of squared invalids
    Civ <- observations_cst[!feasibility,,drop=FALSE]
    if (length(equality) !=1 && sum(equality) !=0) {
      Civ[,which(!equality)] <- pmax(Civ[,which(!equality)], 0)
    } else {
      Civ[Civ <= 0] <- 0
    }
    cm2 <- min(rowSums(Civ^2))
    ## smallest valid objective
    fun.feas <- observations_obj[feasibility]
    ## adapt rho
    if(length(fun.feas) > 0) fm <- max(1e-4, abs(min(fun.feas)))
    else                     fm <- abs(median(observations_obj))
    rho <- abs(cm2/(2*fm))
    
    
    if (is.null(critcontrol)) critcontrol <- list()
    if (is.null(critcontrol$lambda)) critcontrol$lambda <- matrix(rep(1, n.cst),nrow=1)
    if (is.null(dim(critcontrol$lambda))) critcontrol$lambda <- matrix(critcontrol$lambda, nrow=1)
    if (is.null(critcontrol$rho)) critcontrol$rho <- rho
    if (is.null(critcontrol$elit)) critcontrol$elit <- FALSE
    if (is.null(critcontrol$n.mc)) critcontrol$n.mc <- 50
    if (!is.null(critcontrol$slack))      if (critcontrol$slack!=TRUE) critcontrol$slack <- NULL
    mu <- 1/2/critcontrol$rho
    if (!is.null(critcontrol$optimslack)){
      if (critcontrol$optimslack==TRUE && optimcontrol$method=="discrete") {
        cat("optimslack option incompatible with discrete optimization - set to FALSE \n")
        critcontrol$optimslack <- NULL
      }
      if (critcontrol$optimslack!=TRUE) critcontrol$optimslack <- NULL
    } 
    
#     if (!is.null(critcontrol$slack)) {
#       mu <- 1/2/critcontrol$rho
#       critcontrol$slack <- pmax( -matrix(rep(critcontrol$lambda/2/mu, model.fun@n), ncol=n.cst, byrow=TRUE) - observations_cst, 0)
#       critcontrol$slack[,equality] <- 0
#     }
  }
  
  if(!notrace){
    cat("----------------------------\n")
    cat("Starting optimization with : \n The criterion", crit, "\n The solver",  optimcontrol$method, "\n")
    cat("----------------------------\n")
    if (crit=="AL"){
      if (is.null(critcontrol$slack)) cat("Ite | Crit ||  x |  obj |  cst || rho | lambda \n")
      else                            cat("Ite | Crit ||  x |  obj |  cst ||  slack | rho | lambda \n")
    } else {
      cat("Ite | Crit ||  x |  obj |  cst \n")
    }
  }
  
  #### Main loop starts here #############################################
  for (i in 1:nsteps) {
    
    # Update critcontrol for AL
    if (crit=="AL" && i>1) {
      # Update critcontrol if needed
      if (is.null(critcontrol$slack)) {
        lagrang <- model.fun@y + observations_cst%*%t(critcontrol$lambda) + 
          1/2/critcontrol$rho*rowSums(pmax(observations_cst, t(matrix(rep(0,model.fun@n), ncol=model.fun@n)))^2)
        i.best <- which.min(lagrang)
        cst.best <- observations_cst[i.best,]
        critcontrol$lambda <- pmax(critcontrol$lambda + 1/critcontrol$rho*cst.best, 0)
        if (max(cst.best)>0) critcontrol$rho <- critcontrol$rho/2
        
      } else {
        mu <- 1/2/critcontrol$rho
        slack <- pmax( -matrix(rep(critcontrol$lambda/2/mu, model.fun@n), ncol=n.cst, byrow=TRUE) - observations_cst, 0)
        slack[,equality] <- 0
        
        lagrang <- model.fun@y + observations_cst%*%t(critcontrol$lambda) + 
          1/2/critcontrol$rho*rowSums((observations_cst + slack)^2)
        i.best <- which.min(lagrang)
        cst.best <- observations_cst[i.best,]
        critcontrol$lambda <- pmax(critcontrol$lambda + 1/critcontrol$rho*cst.best, 0) #critcontrol$lambda + 1/critcontrol$rho*cst.best
        if (max(cst.best)>0) critcontrol$rho <- critcontrol$rho/2
      }
    }
    
    ## Change the seeds for genoud to avoid selecting always the same initial values
    if(optimcontrol$method == "genoud" & is.null(optimcontrol$unif.seed)){
      optimcontrol$unif.seed <- runif(1)
    }
    sol <- try(critcst_optimizer(crit=crit, model.fun = model.fun, model.constraint = model.constraint,
                                 equality = equality, lower=lower, upper=upper,
                                 optimcontrol=optimcontrol, type=type, 
                                 critcontrol=critcontrol))
    ## Exit if optimization failed
    if (typeof(sol) == "character") {
      if(!notrace){
        cat("Unable to maximize criterion at iteration ", i, "- optimization stopped \n")
        cat("Last model returned \n")
      }
      
      par <- values <- c()
      if (i > 1) {
        par <- model.fun@X[(n+1):(n+i-1),, drop = FALSE]
        values <- model.fun@y[(n+1):(n+i-1)] 
        constraint <- observations_cst[(n+1):(n+i-1),, drop=FALSE]
        feasibility <- feasibility[(n+1):(n+i-1),drop=FALSE]
      }
      return(list(par=par, values=values, constraint=constraint, feasibility=feasibility, nsteps = i-1, lastmodel.fun = model.fun,
                  lastmodel.constraint = model.constraint))
    }
    
    if (sol$value==0 && crit=="AL") {
      # Restart optimization with a proxy criterion
      if(!notrace)
        cat("Optimization failed, restarted with proxy criterion.\n")
     
      sol <- try(critcst_optimizer(crit=crit, model.fun = model.fun, model.constraint = model.constraint,
                                   equality = equality, lower=lower, upper=upper,
                                   optimcontrol=optimcontrol, type=type, 
                                   critcontrol=c(critcontrol, proxy=TRUE)))
      ## Exit if optimization failed
      if (typeof(sol) == "character") {
        if(!notrace){
          cat("Unable to maximize criterion at iteration ", i, "- optimization stopped \n")
          cat("Last model returned \n")
        }
        
        par <- values <- c()
        if (i > 1) {
          par <- model.fun@X[(n+1):(n+i-1),, drop = FALSE]
          values <- model.fun@y[(n+1):(n+i-1)] 
          constraint <- observations_cst[(n+1):(n+i-1),, drop=FALSE]
          feasibility <- feasibility[(n+1):(n+i-1),drop=FALSE]
        }
        return(list(par=par, values=values, constraint=constraint, feasibility=feasibility, nsteps = i-1, lastmodel.fun = model.fun,
                    lastmodel.constraint = model.constraint))
      }
    }
    
    ## Check if optimization do not return already known point
    if(checkPredict(sol$par, model= c(model.constraint, model.fun), type = type, 
                    distance = critcontrol$distance, threshold = critcontrol$threshold) || sol$value <=0){
      if(!notrace)
        cat("Optimization failed, so a random point is selected (consider increasing the inner optimization budget).\n")
      sol$par <- matrix(runif(d), nrow = 1)
    }
    
    ## Update
    X.new <- matrix(as.numeric(sol$par), nrow=1, ncol=d)
    Y.new <- try(fun(as.numeric(sol$par), ...))
    Y.cst.new <- try(constraint(as.numeric(sol$par), ...))
    
    if (typeof(Y.new) == "character") {
      if(!notrace){
        cat("Unable to compute objective function at iteration ", i, "- optimization stopped \n")
        cat("Problem occured for the design: ", X.new, "\n")
        cat("Last model returned \n")
      }
      
      par <- values <- c()
      if (i > 1) {
        par <- model.fun@X[(n+1):model.fun@n,, drop=FALSE]
        values <- model.fun@y[(n+1):model.fun@n]
        constraint <- observations_cst[(n+1):(n+i-1),, drop=FALSE]
        feasibility <- feasibility[(n+1):(n+i-1),drop=FALSE]
      }
      
      return(list(par=par, values=values, constraint=constraint, feasibility=feasibility, nsteps = i-1, lastmodel.fun = model.fun,
                  lastmodel.constraint = model.constraint))
    }
    if (typeof(Y.cst.new) == "character") {
      if(!notrace){
        cat("Unable to compute constraint function at iteration ", i, "- optimization stopped \n")
        cat("Problem occured for the design: ", X.new, "\n")
        cat("Last model returned \n")
      }
      
      par <- values <- c()
      if (i > 1) {
        par <- model.fun@X[(n+1):model.fun@n,, drop=FALSE]
        values <- model.fun@y[(n+1):model.fun@n]
        constraint <- observations_cst[(n+1):(n+i-1),, drop=FALSE]
        feasibility <- feasibility[(n+1):(n+i-1),drop=FALSE]
      }
      return(list(par=par, values=values, constraint=constraint, feasibility=feasibility, nsteps = i-1, lastmodel.fun = model.fun,
                  lastmodel.constraint = model.constraint))
    }
    if(!notrace) {
      if (crit=="AL"){
        if (is.null(sol$slack)) {
          sol$slack <- pmax( -critcontrol$lambda/2/mu - Y.cst.new, 0)
          sol$slack[,equality] <- 0 
        }
        
        if (!is.null(critcontrol$slack) && !is.null(sol$slack)) {
          cat( i, "|", signif(sol$val,3), "||", signif(X.new,3), "|", signif(Y.new,3), "|", signif(Y.cst.new,3), 
               "||", signif(sol$slack,3), "|", signif(critcontrol$rho,3), "|", signif(critcontrol$lambda,3),  "\n")
        } else {
          cat( i, "|", signif(sol$val,3), "||", signif(X.new,3), "|", signif(Y.new,3), "|", signif(Y.cst.new,3), 
               "||", signif(critcontrol$rho,3), "|", signif(critcontrol$lambda,3),  "\n")
        }                           
      } else {
        cat( i, "|", signif(sol$val,3), "||", signif(X.new,3), "|", signif(Y.new,3), "|", signif(Y.cst.new,3), "\n")
      }
    }
    
    # Remove new observation from integration points if discrete case is used
    if (optimcontrol$method=="discrete") {
      optimcontrol$candidate.points <- optimcontrol$candidate.points[-sol$index,,drop=FALSE]
      if (crit=="SUR") { 
        critcontrol$integration.points <- critcontrol$integration.points[-sol$index,,drop=FALSE]
      }
    }
    
    # Update models
    observations_obj <- rbind(observations_obj, Y.new)
    observations_cst <- rbind(observations_cst, Y.cst.new)
    newmodel.fun <- model.fun
    if(any(Y.cst.new > 0)){
      feasibility <- c(feasibility, FALSE)
    }else{
      feasibility <- c(feasibility, TRUE)
    }
    
    #First the objective
    newmodel.fun <- try(update(object = model.fun, newX = X.new, newy=Y.new, newX.alreadyExist=FALSE,
                               cov.reestim = cov.reestim, kmcontrol = list(control = list(trace = FALSE))), silent = TRUE)
    
    if (typeof(newmodel.fun) == "character" && cov.reestim) {
      cat("Error in hyperparameter estimation - old hyperparameter values used instead for the objective model \n")
      newmodel.fun <- try(update(object = model.fun, newX = X.new, newy=Y.new, newX.alreadyExist=FALSE, cov.reestim = FALSE), silent = TRUE)
    }
    if (typeof(newmodel.fun) == "character") {
      cat("Unable to udpate kriging model at iteration", i-1, "- optimization stopped \n")
      cat("lastmodel is the model at iteration", i-1, "\n")
      cat("par and values contain the ",i, "th observation \n \n")
      if (i > 1) allX.new <- rbind(model.fun@X[(n+1):(n+i-1),, drop=FALSE], X.new)
      return(list(
        par    = allX.new,
        values = observations_obj[(n+1):(n+i),drop=FALSE],
        constraint = observations_cst[(n+1):(n+i),, drop=FALSE],
        feasibility = feasibility[(n+1):(n+i),drop=FALSE],
        nsteps = i, 
        lastmodel.fun = model.fun,
        lastmodel.constraint = model.constraint))
    } else {
      model.fun<- newmodel.fun
    }
    
    #then constraints
    newmodel.constraint <- model.constraint
    for (j in 1:n.cst) {
      newmodel.constraint[[j]] <- try(update(object = model.constraint[[j]], newX = X.new, newy=Y.cst.new[j], newX.alreadyExist=FALSE,
                                             cov.reestim = cov.reestim, kmcontrol = list(control = list(trace = FALSE))), silent = TRUE)
      if (typeof(newmodel.constraint[[j]]) == "character" && cov.reestim) {
        cat("Error in hyperparameter estimation - old hyperparameter values used instead for model ", j, "\n")
        newmodel.constraint[[j]] <- try(update(object = model.constraint[[j]], newX = X.new, newy=Y.cst.new[j], newX.alreadyExist=FALSE, cov.reestim = FALSE), silent = TRUE)
      }
      if (typeof(newmodel.constraint[[j]]) == "character") {
        cat("Unable to udpate constraint kriging model at iteration", i-1, "- optimization stopped \n")
        cat("lastmodel.constraint is the model at iteration", i-1, "\n")
        cat("par and values contain the ",i, "th observation \n \n")
        if (i > 1) allX.new <- rbind(model.constraint[[1]]@X[(n+1):(n+i-1),, drop=FALSE], X.new)
        else       allX.new <- X.new
        return(list(
          par    = allX.new,
          values = observations_obj[(n+1):(n+i),drop=FALSE],
          constraint = observations_cst[(n+1):(n+i),, drop=FALSE],
          feasibility = feasibility[(n+1):(n+i),drop=FALSE],
          nsteps = i, 
          lastmodel.fun = model.fun,
          lastmodel.constraint = model.constraint))
      } else {
        model.constraint[[j]] <- newmodel.constraint[[j]]
      }
    }
  }
  
  if(!notrace) cat("\n")
  #### End of main loop ################
  
  
  return(list(
    par=model.fun@X[(n+1):(n+i),, drop=FALSE], 
    values = observations_obj[(n+1):(n+i),drop=FALSE],
    constraint = observations_cst[(n+1):(n+i),, drop=FALSE],
    feasibility = feasibility[(n+1):(n+i),drop=FALSE],
    nsteps=i, 
    lastmodel.fun=model.fun,
    lastmodel.constraint = model.constraint))
}
