##' Given objects of class \code{\link[DiceKriging]{km}} for the objective and constraints,
##' and a set of tuning parameters (\code{lower, upper and critcontrol}), \code{critcst_optimizer} performs
##' the maximization of a constrained Expected Improvement or SUR criterion and delivers
##' the next point to be visited in an EGO-like procedure. \cr \cr
##' The latter maximization relies either on a genetic algorithm using derivatives,
##' \code{\link[rgenoud]{genoud}} or exhaustive search at pre-specified points. 
##' It is important to remark that the information
##' needed about the objective and constraint functions reduces here to the vector of response values
##' embedded in the models (no call to the objective/constraint functions or simulators (except possibly for the objective)).
##'
##' @title Maximization of constrained Expected Improvement criteria
##' 
##' @param crit sampling criterion. Three choices are available : "\code{EFI}", "\code{AL}" and "\code{SUR}",
##' @param model.fun object of class \code{\link[DiceKriging]{km}} corresponding to the objective function,
##' or, if the objective function is fast-to-evaluate, either the objective function to be minimized or a 
##' \code{\link[DOlab]{fastfun}} object, see details and examples below, 
##' @param model.constraint either one or a list of models of class \code{\link[DiceKriging]{km}}, one for each constraint, 
##' @param equality either \code{FALSE} if all constraints are inequalities, or a Boolean vector indicating which are equalities, 
##' @param lower vector of lower bounds for the variables to be optimized over,
##' @param upper vector of upper bounds for the variables to be optimized over,
##' @param type "\code{SK}" or "\code{UK}" (default), depending whether uncertainty related to trend estimation has to be taken into account.
##' @param critcontrol optional list of control parameters for criterion \code{crit}, see details.\cr
##' Options for the \code{\link[DOlab]{checkPredict}} function: \code{threshold} (\code{1e-4}) and \code{distance} (\code{covdist}) are used to avoid numerical issues occuring when adding points too close to the existing ones.
##' @param optimcontrol optional list of control parameters for optimization of the selected infill criterion. 
##'       "\code{method}" set the optimization method; one can 
##'        choose between "\code{discrete}" and "\code{genoud}". For each method, further parameters can be set.\cr 
##'        For "\code{discrete}", one has to provide the argument "\code{candidate.points}". \cr
##'        For "\code{genoud}", one can control, among others, "\code{pop.size}" (default :  \code{[N = 3*2^dim} for \code{dim < 6} and  \code{N = 32*dim} otherwise]),
##' "\code{max.generations}" (\code{12}), "\code{wait.generations}" (\code{2})),
##'  see \code{\link[rgenoud]{genoud}}. Numbers into brackets are the default values.
##'  @return A list with components: 
##'  \itemize{
##'  \item{\code{par}}{: The best set of parameters found,}
##'  \item{\code{value}}{: The value of expected improvement at \code{par}.}
##'  }
##' 
##' @details
##' Extension of the function \code{\link[DiceOptim]{max_EI}} for constrained optimization.\cr
##' 
##' Available infill criteria with \code{crit} are : \cr
##' \itemize{
##' \item Expected Probability of Feasibily (\code{EFI}) \code{\link[DOlab]{crit_EFI}},
##' \item Augmented Lagrangian (\code{AL}) \code{\link[DOlab]{crit_AL}},
##' \item Stepwise Uncertainty Reduction of the excursion volume (\code{SUR}) \code{\link[DOlab]{crit_SUR_cst}}.
##' }
##' Depending on the selected criterion, parameters  can be given with \code{critcontrol}.
##' Also options for \code{\link[DOlab]{checkPredict}} are available.
##' More precisions are given in the corresponding help pages. \cr
##' 
##' If the objective function to minimize is inexpensive, i.e. no need of a kriging model,
##'  then one can provide it in \code{model.obj}, which is handled next with class \code{\link[DOlab]{fastfun}} (or directly as a \code{\link[DOlab]{fastfun}} object). 
##' See example below.
##' 
##' In the case of equality constraints, it is possible to define them with \code{equality}.
##' Additionally, one can modify the tolerance on the constraints using the \code{tolConstraints} component of \code{critcontrol}:
##' an optional vector giving a tolerance for each of the constraints (equality or inequality). 
##' It is highly recommended to use it when there are equality constraints since the default tolerance of 0.05 (resp. 0 for inequality constraints) 
##' in such case might not be suited.\cr
##' 
##' 
##' @author
##' Victor Picheny 
##' 
##' Mickael Binois  
##' 
##' @seealso \code{\link{critcst_optimizer}}, \code{\link{crit_EFI}}, \code{\link{crit_AL}},
##' \code{\link{crit_SUR_cst}}
##' 
##' @references
##' W.R. Jr. Mebane and J.S. Sekhon (2011), Genetic optimization using derivatives: The rgenoud package for R, \emph{Journal of Statistical Software}, 42(11), 1-26 \cr
##' 
##' D.R. Jones, M. Schonlau, and W.J. Welch (1998), Efficient global optimization of expensive black-box functions, \emph{Journal of Global Optimization}, 13, 455-492.
##' @export
##' @importFrom rgenoud genoud
##' @examples
##' #---------------------------------------------------------------------------
##' # 2D objective function, 2 cases
##' #---------------------------------------------------------------------------
##' \donttest{
##' set.seed(2546)
##' library(DiceDesign)
##' 
##' n_var <- 2
##' fun <- branin
##' 
##' fun1.cst <- function(x){return(goldsteinprice(x)+.5)}
##' fun2.cst <- function(x){return(3/2 - x[1] - 2*x[2] - .5*sin(2*pi*(x[1]^2 - 2*x[2])))}
##' 
##' # Constraint function with vectorial output
##' cstfun <- function(x){return(c(fun1.cst(x), fun2.cst(x)))}
##' 
##' n.grid <- 31
##' test.grid <- expand.grid(X1 = seq(0, 1, length.out = n.grid), X2 = seq(0, 1, length.out = n.grid))
##' obj.grid <- apply(test.grid, 1, fun)
##' cst1.grid <- apply(test.grid, 1, fun1.cst)
##' cst2.grid <- apply(test.grid, 1, fun2.cst)
##' 
##' n_appr <- 12 
##' design.grid <- round(maximinESE_LHS(lhsDesign(n_appr, n_var, seed = 2)$design)$design, 1)
##' obj.init <- apply(design.grid, 1, fun)
##' cst1.init <- apply(design.grid, 1, fun1.cst)
##' cst2.init <- apply(design.grid, 1, fun2.cst)
##' model.fun <- km(~., design = design.grid, response = obj.init)
##' model.constraint1 <- km(~., design = design.grid, response = cst1.init, lower=c(.2,.2))
##' model.constraint2 <- km(~., design = design.grid, response = cst2.init, lower=c(.2,.2))
##' models.cst <- list(model.constraint1, model.constraint2)
##' 
##' lower <- rep(0, n_var)
##' upper <- rep(1, n_var)
##' 
##' #---------------------------------------------------------------------------
##' # Augmented Lagrangian Improvement, fast objective function, two ineq constraints,
##' # optimization with genoud
##' #---------------------------------------------------------------------------
##' critcontrol <- list(lambda=c(.5,2), rho=.5)
##' optimcontrol <- list(method = "genoud", max.generations=10, pop.size=20)
##' 
##' AL_grid <- apply(test.grid, 1, crit_AL, model.fun = fastfun(fun, design.grid),
##'                  model.constraint = models.cst, critcontrol=critcontrol)
##' 
##' cstEGO1 <- critcst_optimizer(crit = "AL", model.fun = fun,
##'                              model.constraint = models.cst, equality = FALSE,
##'                              lower = lower, upper = upper, 
##'                              optimcontrol = optimcontrol, critcontrol=critcontrol)
##' 
##' filled.contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), nlevels = 50,
##'                matrix(AL_grid, n.grid), main = "AL map and its maximizer (blue)",
##'                xlab = expression(x[1]), ylab = expression(x[2]), color = terrain.colors, 
##'                plot.axes = {axis(1); axis(2);
##'                             points(design.grid[,1], design.grid[,2], pch = 21, bg = "white")
##'                             contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), 
##'                             matrix(obj.grid, n.grid), nlevels = 10, add=TRUE,drawlabels=TRUE,
##'                                    col = "black")
##'                             contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), 
##'                             matrix(cst1.grid, n.grid), level = 0, add=TRUE,drawlabels=FALSE,
##'                                    lwd=1.5, col = "red")
##'                             contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), 
##'                             matrix(cst2.grid, n.grid), level = 0, add=TRUE,drawlabels=FALSE,
##'                                    lwd=1.5, col = "red")
##'                             points(cstEGO1$par, col = "blue", pch = 4, lwd = 2)
##'                             }
##'               )
##' #---------------------------------------------------------------------------
##' # SUR, expensive objective function, one equality constraint,
##' # optimization with genoud, integration on a regular grid
##' #---------------------------------------------------------------------------
##' optimcontrol <- list(method = "genoud", s = 40, maxit = 40)
##' critcontrol  <- list(tolConstraints = .15, integration.points=as.matrix(test.grid))
##' 
##' SUR_grid <- apply(test.grid, 1, crit_SUR_cst, model.fun = model.fun,
##'                   model.constraint = model.constraint1, equality = TRUE, critcontrol = critcontrol)
##' 
##' cstEGO2 <- critcst_optimizer(crit = "SUR", model.fun = model.fun,
##'                              model.constraint = model.constraint1, equality = TRUE,
##'                              lower = lower, upper = upper, 
##'                              optimcontrol = optimcontrol, critcontrol = critcontrol)
##' 
##' filled.contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), nlevels = 50,
##'                matrix(SUR_grid, n.grid), main = "SUR map and its maximizer (blue)",
##'                xlab = expression(x[1]), ylab = expression(x[2]), color = terrain.colors, 
##'                plot.axes = {axis(1); axis(2);
##'                             points(design.grid[,1], design.grid[,2], pch = 21, bg = "white")
##'                             contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), 
##'                             matrix(obj.grid, n.grid), nlevels = 10, add=TRUE,
##'                             drawlabels=TRUE, col = "black")
##'                             contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), 
##'                             matrix(cst1.grid, n.grid), level = c(-critcontrol$tolConstraints,
##'                             critcontrol$tolConstraints), 
##'                             add=TRUE, drawlabels=FALSE,lwd=1.5, col = "orange")
##'                             points(cstEGO2$par, col = "blue", pch = 4, lwd = 2)
##'                             }
##'               )
##' }
critcst_optimizer <- function(crit = "EFI", model.fun, model.constraint, equality = FALSE, lower, upper, type = "UK", 
                              critcontrol = NULL, optimcontrol = NULL){
  ###########################################################################################
  # Finds the maximizer of the criterion
  ###########################################################################################
  if(is.null(optimcontrol$method)) optimcontrol$method <- "genoud"

  n.cst <- length(model.constraint)
  if (n.cst!=length(equality)) equality <- rep(FALSE, n.cst)
  n.ineq <- sum(!equality)
  
  if(crit == "SUR" && n.cst > 3){
    cat("crit_SUR_cst does not take more than 3 constraints \n")
    return(NA)
  }
  
  if(n.cst == 1 && class(model.constraint) != "list"){
    model.constraint <- list(model.constraint)
  }
  
  if(class(model.fun) == "km"){
    
  }else if(class(model.fun) == "function"){
    #     fastobs <- apply(model.constraint[[1]]@X, 1, cheapfun)
    model.fun <- fastfun(fn = model.fun, design = model.constraint[[1]]@X)
    # cat("Fastfun-Mode on \n")
  }else if(class(model.fun) == "fastfun"){
    # cat("Fastfun-Mode on \n")
  }else{
    cat("model.fun should be either a km, fonction or fastfun object \n")
    return(NA)
  }
  #   if((!is.null(model.fun) && class(model.fun) != "fastfun") && !is.null(cheapfun)){
  #     cat("Either model.fun or cheapfun should be set to null (or model.fun can be a fastfun).")
  #     return(0)
  #   }
  #   
  #   if (!is.null(cheapfun)){
  #     fastobs <- apply(model.constraint[[1]]@X, 1, cheapfun)
  #     if(is.null(model.fun))
  #       model.fun <- fastfun(fn = cheapfun, design = model.constraint[[1]]@X, response = fastobs)
  #   }
  
  if(!is.null(critcontrol) && !is.null(critcontrol$tolConstraints)){
    if(min(critcontrol$tolConstraints) < 0){
      cat("tolConstraints has negative components, this may not work with equality contraints \n")
    }
  }
  
  d <- model.fun@d
  
  # Different calls for crit_optimizer, depending on the criterion chosen
  if (crit=="EFI"){
    #-------------------------------------------------------
    criterion <- crit_EFI
    #-------------------------------------------------------
  } else if (crit=="AL"){
    if (is.null(critcontrol))             critcontrol <- list()
    if (is.null(critcontrol$lambda))      critcontrol$lambda <- matrix(rep(1,n.cst), nrow=1)
    if (is.null(dim(critcontrol$lambda))) critcontrol$lambda <- matrix(critcontrol$lambda, nrow=1)
    if (is.null(critcontrol$rho))         critcontrol$rho <- .5
    if (is.null(critcontrol$elit))        critcontrol$elit <- FALSE
    if (!is.null(critcontrol$slack))      if (critcontrol$slack!=TRUE) critcontrol$slack <- NULL
    if (!is.null(critcontrol$optimslack)) if (critcontrol$optimslack!=TRUE) critcontrol$optimslack <- NULL
    
    criterion <- crit_AL
    #     if (!is.null(critcontrol$slack)) {
    #       observations_cst <- c()
    #       for (i in 1:n.cst) observations_cst <- cbind(observations_cst, model.constraint[[i]]@y)
    #       mu <- 1/2/critcontrol$rho
    #       critcontrol$slack <- pmax( -matrix(rep(critcontrol$lambda/2/mu, model.fun@n), ncol=n.cst, byrow=TRUE) - observations_cst, 0)
    #       critcontrol$slack[,equality] <- 0
    #       
    #       criterion <- crit_AL
    #     }
    
    #-------------------------------------------------------
  } else if (crit=="SUR"){
    criterion <- crit_SUR_cst
    
    # Format integration points and weights if needed
    integration.param <- integration_design_cst(integcontrol=critcontrol, lower=lower, upper=upper, model.fun=model.fun,
                                                model.constraint=model.constraint)
    critcontrol$integration.points  <- (integration.param$integration.points)
    critcontrol$integration.weights <- (integration.param$integration.weights)
    
    # Precompute necessary quantities if not provided
    if(is.null(critcontrol$precalc.data.cst) || is.null(critcontrol$mn.X.cst) || is.null(critcontrol$sn.X.cst) ||
       (is.null(critcontrol$precalc.data.obj)&&class(model.fun)=="km") || is.null(critcontrol$mn.X.obj) || is.null(critcontrol$sn.X.obj) ){
      precalc.data.cst <- vector("list", n.cst)
      intpoints.oldmean.cst <- intpoints.oldsd.cst <- matrix(0, n.cst, nrow(critcontrol$integration.points))
      for (i in 1:n.cst){
        p.tst <- predict(model.constraint[[i]], newdata=critcontrol$integration.points, type=type, checkNames=FALSE)
        intpoints.oldmean.cst[i,] <- p.tst$mean
        intpoints.oldsd.cst[i,]   <- p.tst$sd
        precalc.data.cst[[i]] <- precomputeUpdateData(model.constraint[[i]], critcontrol$integration.points)
      }
      p.tst <- predict(model.fun, newdata=critcontrol$integration.points, type=type, checkNames=FALSE)
      intpoints.oldmean.obj <- p.tst$mean
      intpoints.oldsd.obj   <- p.tst$sd
      if (class(model.fun)=="km") precalc.data.obj <- precomputeUpdateData(model.fun, critcontrol$integration.points)
      else                        precalc.data.obj <- NULL
      critcontrol <- c(critcontrol, list(mn.X.cst=intpoints.oldmean.cst, sn.X.cst=intpoints.oldsd.cst, precalc.data.cst=precalc.data.cst,
                                         mn.X.obj=intpoints.oldmean.obj, sn.X.obj=intpoints.oldsd.obj, precalc.data.obj=precalc.data.obj))
    }
  }

  ########################################################################################
  ## Discrete Optimisation
  ########################################################################################
  
  if(optimcontrol$method=="discrete"){
    optim.points <- optimcontrol$candidate.points
    colnames(optim.points) <- colnames(model.fun@X)
    n.optim.points <- nrow(optim.points)
    all.crit <- rep(0, n.optim.points)
    critcontroldiscrete <- critcontrol
    
    for (k in 1:n.optim.points){
      if (crit=="SUR"){
        critcontroldiscrete$integration.points <- critcontrol$integration.points[-k,,drop=FALSE]
        if (!is.null(critcontroldiscrete$integration.weights)) critcontroldiscrete$integration.weights <- critcontrol$integration.weights[-k]
        critcontroldiscrete$mn.X.cst <- critcontrol$mn.X.cst[,-k,drop=FALSE]
        critcontroldiscrete$sn.X.cst <- critcontrol$sn.X.cst[,-k,drop=FALSE]
        for (i in 1:n.cst){
          #                 if (!is.null(precalc.data.cst[[i]])) {
          critcontroldiscrete$precalc.data.cst[[i]]$Kinv.c.olddata <- critcontrol$precalc.data.cst[[i]]$Kinv.c.olddata[,-k,drop=FALSE]
          critcontroldiscrete$precalc.data.cst[[i]]$first.member   <- critcontrol$precalc.data.cst[[i]]$first.member[-k]
          #                 }
        }
        
        critcontroldiscrete$mn.X.obj <- critcontrol$mn.X.obj[-k,drop=FALSE]
        critcontroldiscrete$sn.X.obj <- critcontrol$sn.X.obj[-k,drop=FALSE]
        if (class(model.fun)=="km") {
          critcontroldiscrete$precalc.data.obj$Kinv.c.olddata <- critcontrol$precalc.data.obj$Kinv.c.olddata[,-k,drop=FALSE]
          critcontroldiscrete$precalc.data.obj$first.member   <- critcontrol$precalc.data.obj$first.member[-k]
        }
      }
      all.crit[k] <- criterion(x=as.numeric(optim.points[k,]), model.fun=model.fun,
                               model.constraint = model.constraint, equality = equality, type=type,
                               critcontrol=critcontroldiscrete)
    }
    all.crit[is.na(all.crit)] <- 0
    value  <- max(all.crit)
    i.best <- which(all.crit==value)[sample.int(length(which(all.crit==value)), 1)]
    par    <- matrix(optim.points[i.best,], nrow=1, ncol=model.fun@d)
    colnames(par) <- colnames(model.fun@X)
    
    return(list(par=par, value=value, index=i.best, slack=NULL))
  }
  
  ########################################################################################
  # Special parameters for the slack variables
  if (!is.null(critcontrol$slack) && !is.null(critcontrol$optimslack)) {
    d <- d + n.ineq
    upper   <- c(upper, rep(1, n.ineq))
    lower   <- c(lower, rep(0, n.ineq))
  }
  ########################################################################################
  ## Optimization with Genoud
  ########################################################################################
  if(optimcontrol$method=="genoud"){
    
    if (d <= 6) N <- 3 * 2^d else N <- 32 * d
    
    if (is.null(optimcontrol$pop.size))         optimcontrol$pop.size <- N
    if (is.null(optimcontrol$max.generations))  optimcontrol$max.generations <- 12
    if (is.null(optimcontrol$wait.generations)) optimcontrol$wait.generations <- 2
    if (is.null(optimcontrol$BFGSburnin))       optimcontrol$BFGSburnin <- 2
    if (is.null(optimcontrol$parinit))          optimcontrol$parinit <- lower + runif(d) * (upper - lower)
    if (is.null(optimcontrol$unif.seed))        optimcontrol$unif.seed <- 1
    if (is.null(optimcontrol$int.seed))         optimcontrol$int.seed <- 1
    if (is.null(optimcontrol$print.level))           optimcontrol$print.level <- 0
    if (is.null(optimcontrol$BFGSmaxit))             optimcontrol$BFGSmaxit <- N
    if (is.null(optimcontrol$solution.tolerance))    optimcontrol$solution.tolerance <- 1e-21
    
    # Mutations
    if (is.null(optimcontrol$P1)) optimcontrol$P1<-50
    if (is.null(optimcontrol$P2)) optimcontrol$P2<-50
    if (is.null(optimcontrol$P3)) optimcontrol$P3<-50
    if (is.null(optimcontrol$P4)) optimcontrol$P4<-50
    if (is.null(optimcontrol$P5)) optimcontrol$P5<-50
    if (is.null(optimcontrol$P6)) optimcontrol$P6<-50
    if (is.null(optimcontrol$P7)) optimcontrol$P7<-50
    if (is.null(optimcontrol$P8)) optimcontrol$P8<-50
    if (is.null(optimcontrol$P9)) optimcontrol$P9<-0
    
    domaine <- cbind(lower, upper)

    o <- genoud(fn=criterion, nvars=d, max=TRUE, pop.size=optimcontrol$pop.size,
                max.generations=optimcontrol$max.generations,wait.generations=optimcontrol$wait.generations,
                hard.generation.limit=TRUE, starting.values=optimcontrol$parinit, MemoryMatrix=TRUE,
                Domains=domaine, default.domains=10, solution.tolerance=optimcontrol$solution.tolerance,
                boundary.enforcement=2, lexical=FALSE, gradient.check=FALSE, BFGS=TRUE,
                data.type.int=FALSE, hessian=FALSE, unif.seed=optimcontrol$unif.seed, 
                int.seed=optimcontrol$int.seed,print.level=optimcontrol$print.level, share.type=0, instance.number=0,
                output.path="stdout", output.append=FALSE, project.path=NULL,
                P1=optimcontrol$P1, P2=optimcontrol$P2, P3=optimcontrol$P3, 
                P4=optimcontrol$P4, P5=optimcontrol$P5, P6=optimcontrol$P6,
                P7=optimcontrol$P7, P8=optimcontrol$P8, P9=optimcontrol$P9,
                P9mix=NULL, BFGSburnin=optimcontrol$BFGSburnin,BFGSfn=NULL, BFGShelp=NULL,
                control = list(maxit = optimcontrol$BFGSmaxit),
                cluster=FALSE, balance=FALSE, debug=FALSE,
                model.fun = model.fun, model.constraint = model.constraint, equality = equality, type=type, 
                critcontrol=critcontrol)
    
    par <- t(as.matrix(o$par))
    slack <- NULL
    if (!is.null(critcontrol$slack) && !is.null(critcontrol$optimslack)) {
      slack <- par[1, model.fun@d + (1:n.ineq)]
      par   <- par[1, 1:model.fun@d,drop=FALSE]
    }
    colnames(par) <- colnames(model.fun@X)
    value <- as.matrix(o$value)
    return(list(par=par, value=value, slack=slack))
  }
#   ########################################################################################
#   ## Optimization with PSO
#   ########################################################################################
#   if(optimcontrol$method=="pso"){
#     
#     control <- list(fnscale=-1, maxit=optimcontrol$maxit, s = optimcontrol$s)
#     if (is.null(control$maxit))   control$maxit=400
#     if (is.null(control$s)) control$s = floor(10+2*sqrt(length(d)))
#     
#     o <- psoptim(par = rep(NA, d) , criterion, lower = lower, upper = upper, control = control,
#                  model.fun = model.fun, model.constraint = model.constraint, equality = equality, type=type, 
#                  critcontrol=critcontrol)
#     par <- t(as.matrix(o$par))
#     slack <- NULL
#     if (!is.null(critcontrol$slack)) {
#       if (!is.null(critcontrol$optimslack)) {
#         slack <- par[1, model.fun@d + (1:n.ineq)]
#         par   <- par[1, 1:model.fun@d,drop=FALSE]
#       }
#     }
#     colnames(par) <- colnames(model.fun@X)
#     value <- as.matrix(o$value)
#     return(list(par=par, value=-value, slack=slack))
#   }
}