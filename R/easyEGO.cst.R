##' User-friendly wrapper of the function \code{\link[DOlab]{EGO.cst}}
##' Generates initial DOEs and kriging models (objects of class \code{\link[DiceKriging]{km}}), 
##' and executes \code{nsteps} iterations of EGO methods integrating constraints.
##' @title EGO algorithm with constraints
##' @details Does not require knowledge on kriging models (objects of class \code{\link[DiceKriging]{km}})  \cr
##' 
##' The problem considered is of the form: \eqn{min f(x)} s.t. \eqn{g(x) \le 0}, 
##' \eqn{g} having a vectorial output. 
##' By default all its components are supposed to be inequalities, but one can use a Boolean vector in \code{equality} 
##' to specify which are equality constraints, hence of the type \eqn{g(x) = 0}.
##' The \code{control} argument is a list that can supply any of the following components: \cr
##' \itemize{
##' \item \code{method}{: choice of constrained improvement function: "\code{AL}", "\code{EFI}" or "\code{SUR}" 
##' (see \code{\link[DOlab]{crit_EFI}}, \code{\link[DOlab]{crit_AL}}, \code{\link[DOlab]{crit_SUR_cst}})}
##' \item \code{trace}{:  if positive, tracing information on the progress of the optimization is produced.}
##' \item \code{inneroptim}{: choice of the inner optimization algorithm: "\code{genoud}" or "\code{random}"
##'  (see \code{\link[rgenoud]{genoud}}).}
##' \item \code{maxit}{: maximum number of iterations of the inner loop. }
##' \item \code{seed}{: to fix the random variable generator}
##' }
##' For additional details, see \code{\link[DOlab]{EGO.cst}}.
##' 
##' @param  fun scalar function to be minimized,
##' @param  constraint vectorial function corresponding to the constraints, see details below,
##' @param  n.cst number of constraints,
##' @param  budget total number of calls to the objective and constraint functions,
##' @param  lower vector of lower bounds for the variables to be optimized over,
##' @param  upper vector of upper bounds for the variables to be optimized over,
##' @param  cheapfun optional boolean, \code{TRUE} if the objective is a fast-to-evaluate function that does not need a kriging model
##' @param  equality either \code{FALSE} if all constraints are inequalities, else a Boolean vector indicating which are equalities 
##' @param  X initial design of experiments. If not provided, X is taken as a maximin LHD with budget/3 points
##' @param  y initial set of objective observations \eqn{f(X)}. Computed if not provided.
##' @param  C initial set of constraint observations \eqn{g(X)}. Computed if not provided.
##' @param control an optional list of control parameters. See "Details".
##' @param ... additional parameters to be given to BOTH the objective \code{fun} and \code{constraints}.
##' @export
##' @return
##' A list with components:
##' \itemize{
##' \item{\code{par}}{: the best feasible point}
##' \item{\code{values}}{: a vector of the objective and constraints at the point given in \code{par},}
##' \item{\code{history}}{: a list containing all the points visited by the algorithm (\code{X}) and their corresponding objectives (\code{y}) and constraints (\code{C}) \cr \cr
##' If no feasible point is found, \code{par} returns the most feasible point (in the least square sense).}
##' }
##' 
##' 
##' @author
##' Victor Picheny 
##' 
##' Mickael Binois 
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
##' @importFrom DiceDesign maximinESE_LHS lhsDesign
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
##' constraint <- function(x){return(c(fun1.cst(x), fun2.cst(x)))}
##' 
##' # For illustration purposes
##' n.grid <- 31
##' test.grid <- expand.grid(X1 = seq(0, 1, length.out = n.grid), X2 = seq(0, 1, length.out = n.grid))
##' obj.grid <- apply(test.grid, 1, fun)
##' cst1.grid <- apply(test.grid, 1, fun1.cst)
##' cst2.grid <- apply(test.grid, 1, fun2.cst)
##' 
##' lower <- rep(0, n_var)
##' upper <- rep(1, n_var)
##' 
##' #---------------------------------------------------------------------------
##' # 1- Expected Feasible Improvement criterion, expensive objective function,
##' # two inequality constraints, 15 observations budget, using genoud
##' #---------------------------------------------------------------------------
##' res <- easyEGO.cst(fun=fun, constraint=constraint, n.cst=2, lower=lower, upper=upper, budget=15, 
##'                    control=list(method="EFI", inneroptim="genoud", maxit=20))
##' 
##' cat("best design found:", res$par, "\n")
##' cat("corresponding objective and constraints:", res$value, "\n")
##' 
##' # Objective function in colour, constraint boundaries in red
##' # Initial DoE: white circles, added points: blue crosses, best solution: red cross
##' 
##' filled.contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), nlevels = 50,
##'                matrix(obj.grid, n.grid), main = "Two inequality constraints",
##'                xlab = expression(x[1]), ylab = expression(x[2]), color = terrain.colors, 
##'                plot.axes = {axis(1); axis(2);
##'                             points(design.grid[,1], design.grid[,2], pch = 21, bg = "white")
##'                             contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), 
##'                                     matrix(cst1.grid, n.grid), level = 0, add=TRUE,
##'                                     drawlabels=FALSE, lwd=1.5, col = "red")
##'                             contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), 
##'                                     matrix(cst2.grid, n.grid), level = 0, add=TRUE,drawlabels=FALSE,
##'                                     lwd=1.5, col = "red")
##'                             points(res$history$X, col = "blue", pch = 4, lwd = 2)       
##'                             points(res$par[1], res$par[2], col = "red", pch = 4, lwd = 2, cex=2) 
##'                }
##' )
##' 
##' #---------------------------------------------------------------------------
##' # 2- Augmented Lagrangian Improvement criterion, expensive objective function,
##' # one inequality and one equality constraint, 25 observations budget, using random search
##' #---------------------------------------------------------------------------
##' res2 <- easyEGO.cst(fun=fun, constraint=constraint, n.cst=2, lower=lower, upper=upper, budget=25,
##'                    equality = c(TRUE, FALSE),
##'                    control=list(method="AL", inneroptim="random", maxit=100))
##' 
##' cat("best design found:", res2$par, "\n")
##' cat("corresponding objective and constraints:", res2$value, "\n")
##' 
##' # Objective function in colour, inequality constraint boundary in red, equality
##' # constraint in orange
##' # Initial DoE: white circles, added points: blue crosses, best solution: red cross
##' 
##' filled.contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), nlevels = 50,
##'                matrix(obj.grid, n.grid), xlab = expression(x[1]), ylab = expression(x[2]),
##'                main = "Inequality (red) and equality (orange) constraints", color = terrain.colors, 
##'                plot.axes = {axis(1); axis(2);
##'                             points(design.grid[,1], design.grid[,2], pch = 21, bg = "white")
##'                             contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), 
##'                                     matrix(cst1.grid, n.grid), level = 0, add=TRUE,
##'                                     drawlabels=FALSE,lwd=1.5, col = "orange")
##'                             contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid),
##'                                     matrix(cst2.grid, n.grid), level = 0, add=TRUE,
##'                                     drawlabels=FALSE,lwd=1.5, col = "red")
##'                             points(res2$history$X, col = "blue", pch = 4, lwd = 2)
##'                             points(res2$par[1], res2$par[2], col = "red", pch = 4, lwd = 2, cex=2) 
##'                }
##' )
##' 
##' #---------------------------------------------------------------------------
##' # 3- Stepwise Uncertainty Reduction criterion, fast objective function,
##' # single inequality constraint, with initial DOE given + 10 observations,
##' # using genoud
##' #---------------------------------------------------------------------------
##' n.init <- 12 
##' design.grid <- round(maximinESE_LHS(lhsDesign(n.init, n_var, seed = 42)$design)$design, 1)
##' cst2.init <- apply(design.grid, 1, fun2.cst)
##' 
##' res3 <- easyEGO.cst(fun=fun, constraint=fun2.cst, n.cst=1, lower=lower, upper=upper, budget=10,
##'                     X=design.grid, C=cst2.init,
##'                    cheapfun=TRUE, control=list(method="SUR", inneroptim="genoud", maxit=20))
##'         
##' cat("best design found:", res3$par, "\n")
##' cat("corresponding objective and constraint:", res3$value, "\n")
##' 
##' # Objective function in colour, inequality constraint boundary in red
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
##'                                     matrix(cst2.grid, n.grid), level = 0, add = TRUE,
##'                                     drawlabels = FALSE,lwd = 1.5, col = "red")
##'                             points(res3$history$X, col = "blue", pch = 4, lwd = 2)
##'                             points(res3$par[1], res3$par[2], col = "red", pch = 4, lwd = 2, cex=2) 
##'                                            }
##'                 )
##' }
##' 
easyEGO.cst <- function (fun, constraint, n.cst=1, budget, lower, upper, 
                         cheapfun = FALSE, equality = FALSE, X=NULL, y=NULL, C=NULL,
                         control=list(method="EFI", trace=1, inneroptim="genoud", maxit=100, seed=42), ...) {
  
  if (is.null(control$method)) control$method <- "EFI"
  if (is.null(control$trace)) control$trace   <- 1
  if (is.null(control$inneroptim))  control$inneroptim <- "genoud"
  if (is.null(control$maxit)) control$maxit   <- 100
  if (is.null(control$seed)) control$seed   <- 100
  
  dim <- length(lower)
  
  if (length(lower) != length(upper)) {
    cat("Bound values lower and upper are not consistent. Both should be vectors of size d.")
    return(0)
  }
  
  if (!is.null(X)) {
    design.init <- data.frame(x=X) #, ncol=dim, byrow=TRUE)
    temp <- dim(design.init)
    
    if (temp[2] != dim) {
      cat("Bound values (lower and upper) and initial DOE (X) are not consistent. \n 
          lower and upper should be vectors of size d and \n 
          X either a matrix with d columns or a data frame of d variables.")
      return(0)
    }
    n.init <- temp[1]
  } else {
    n.init <- max(4*dim, round(budget/3))
    design.init <- data.frame(x=lower + (upper-lower)*maximinESE_LHS(lhsDesign(n.init, dim, seed=control$seed)$design)$design)
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
  if (!is.null(X) && !is.null(C)) {
    obs.cst.init <- matrix(C, ncol=n.cst)
    if (nrow(obs.cst.init) != n.init) {
      cat("Initial DOE (X) and constraint (C) are not consistent, or n.cst is incorrect.")
      return(0)
    }
  } else {
    if (n.cst == 1) obs.cst.init <- matrix(apply(design.init, 1, constraint, ...),ncol=1)
    else           obs.cst.init <- t(apply(design.init, 1, constraint, ...))
    if (ncol(obs.cst.init) != n.cst || nrow(obs.cst.init) != n.init) {
      cat("n.cst does not match the number of outputs of the constraint function.")
      return(0)
    }
  }
  #----------------------------------------
  if (!is.null(X) && !is.null(y) && !is.null(C)) {
    n.ite <- budget
  } else {
    n.ite <- budget - n.init
  }
  #----------------------------------------
  if (!cheapfun) {
    model.fun <- km(~., design = design.init, response = obs.obj.init, control=list(trace=FALSE), lower=rep(.1,dim), upper=rep(1,dim))
    cheapfun  <- NULL
  } else {
    model.fun <- NULL
    cheapfun  <- fun
  }
  
  model.constraint <- vector("list", n.cst)
  for (j in 1:(n.cst)) {
    model.constraint[[j]] <- km(~., design = design.init, response = obs.cst.init[,j], control=list(trace=FALSE), lower=rep(.1,dim), upper=rep(1,dim))
  }
  #----------------------------------------
  critcontrol <- NULL
  if (any(equality)) {
    tolConstraints <- rep(0, n.cst)
    for (j in which(equality)) tolConstraints[j] <- 0.05*sd(model.constraint[[j]]@y)
    critcontrol <- list(tolConstraints = tolConstraints)
  }
  if (control$inneroptim=="genoud") optimcontrol <- list(method="genoud", max.generations=control$maxit, threshold = 1e-5, distance = "euclidean", notrace = !control$trace>0)
  # if (control$inneroptim=="pso")    optimcontrol <- list(method="pso", maxit=control$maxit, threshold = 1e-5, distance = "euclidean", notrace = !control$trace>0)
  if (control$inneroptim=="random") {
    candidate.points <- matrix(rep(lower,control$maxit)+rep(upper-lower,control$maxit)*runif(control$maxit*dim),byrow=T,ncol=dim)
    optimcontrol = list(method="discrete", candidate.points=candidate.points, threshold = 1e-5, distance = "euclidean", notrace = !control$trace>0)
  }
  
  if (control$method=="SUR") {
    critcontrol$distrib <- "SUR"
    critcontrol$n.points <- 50*dim
  }
  if (control$method=="AL") {
    critcontrol$always.update <- TRUE
  }
  cstEGO <- EGO.cst(model.fun = model.fun, fun = fun, model.constraint = model.constraint,
                    constraint = constraint, crit=control$method, nsteps=n.ite,
                    equality = equality, cheapfun = cheapfun,
                    lower=lower, upper=upper, cov.reestim=TRUE, optimcontrol=optimcontrol,
                    critcontrol = critcontrol, ...)

  allX <- cstEGO$lastmodel.fun@X
  ally <- cstEGO$lastmodel.fun@y
  
  if (n.cst==1) allC <- matrix(cstEGO$lastmodel.constraint[[1]]@y,ncol=1) #c(obs.cst.init, cstEGO$constraints)
  else          allC <- rbind(obs.cst.init, cstEGO$constraint)

  # Compute current best
  feasibility <- test_feas_vec(cst=allC, equality=equality, tolConstraints = critcontrol$tolConstraints)

  if (any(feasibility)) {
    Xtemp <- allX[feasibility,,drop=FALSE]
    ytemp <- ally[feasibility]
    Ctemp <- allC[feasibility,,drop=FALSE]
    par   <- Xtemp[which.min(ytemp),]
    value <- c(min(ytemp), Ctemp[which.min(ytemp),])
  } else {
    # Return least infeasible points
    if (control$trace>0) cat("No feasible point found - least violating point returned instead \n")
    I <- which.min(colSums(allC^2))
    par   <- allX[I,]
    value <- c(ally[I], allC[I,])
  }

  return(list(par=par, value = value, history=list(X=allX, y=ally, C=allC)))
  }