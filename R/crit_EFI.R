##' Computes the Expected Feasible Improvement at current location. 
##' The current feasible minimum of the observations can be replaced by an arbitrary value (plugin), which is usefull in particular in noisy frameworks.
##' 
##' @title Expected Feasible Improvement
##' 
##' @param x a vector representing the input for which one wishes to calculate \code{EFI},
##' @param model.fun object of class \code{\link[DiceKriging]{km}} corresponding to the objective function,
##' or, if the objective function is fast-to-evaluate, a \code{\link[DOlab]{fastfun}} object,
##' @param model.constraint either one or a list of objects of class \code{\link[DiceKriging]{km}}, one for each constraint function,
##' @param equality either \code{FALSE} if all constraints are for inequalities, else a vector of boolean indicating which are equalities, 
##' @param critcontrol optional list with argument \code{tolConstraints}, an optional vector giving a tolerance (> 0) for each of the constraints (equality or inequality).
##' It is highly recommended to use it when there are equality constraints since the default tolerance of 0.05 in such case might not be suited.\cr
##' 
##' Options for the \code{\link[DOlab]{checkPredict}} function: \code{threshold} (\code{1e-4}) and \code{distance} (\code{covdist}) are used to avoid numerical issues occuring when adding points too close to the existing ones.
##' @param plugin optional scalar: if provided, it replaces the feasible minimum of the current observations.
##' If set to \code{Inf}, e.g. when their is no feasible solution, then the criterion is equal to the probability of feasibility,
##' @param type "\code{SK}" or "\code{UK}" (by default), depending whether uncertainty related to trend estimation 
##'        has to be taken into account. 
##' @return The Expected Feasible Improvement at \code{x}.
##' @seealso \code{\link[DiceOptim]{EI}} from package DiceOptim, \code{\link[DOlab]{crit_AL}}, \code{\link[DOlab]{crit_SUR_cst}}.
##' 
##' @export
##' 
##' @author
##' Victor Picheny 
##' 
##' Mickael Binois 
##' 
##' @references 
##' M. Schonlau, W.J. Welch, and D.R. Jones (1998),
##'  Global versus local search in constrained optimization of computer models,
##'  \emph{Lecture Notes-Monograph Series}, 11-25.
##'  
##' M.J. Sasena, P. Papalambros, and P.Goovaerts (2002),
##'  Exploration of metamodeling sampling criteria for constrained global optimization,
##'  \emph{Engineering optimization}, 34, 263-278.
##' 
##' @examples
##' #---------------------------------------------------------------------------
##' # Expected Feasible Improvement surface with one inequality constraint
##' #---------------------------------------------------------------------------
##' \donttest{
##' set.seed(25468)
##' library(DiceDesign)
##' 
##' n_var <- 2 
##' fun.obj <- goldsteinprice
##' fun.cst <- function(x){return(-branin(x) + 25)}
##' n.grid <- 51
##' test.grid <- expand.grid(X1 = seq(0, 1, length.out = n.grid), X2 = seq(0, 1, length.out = n.grid))
##' obj.grid <- apply(test.grid, 1, fun.obj)
##' cst.grid <- apply(test.grid, 1, fun.cst)
##' n.init <- 15 
##' design.grid <- round(maximinESE_LHS(lhsDesign(n.init, n_var, seed = 42)$design)$design, 1)
##' obj.init <- apply(design.grid, 1, fun.obj)
##' cst.init <- apply(design.grid, 1, fun.cst)
##' model.fun <- km(~., design = design.grid, response = obj.init)
##' model.constraint <- km(~., design = design.grid, response = cst.init)
##' 
##' EFI_grid <- apply(test.grid, 1, crit_EFI, model.fun = model.fun,
##'                   model.constraint = model.constraint)
##' 
##' filled.contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), nlevels = 50,
##'                matrix(EFI_grid, n.grid), main = "Expected Feasible Improvement",
##'                xlab = expression(x[1]), ylab = expression(x[2]), color = terrain.colors, 
##'                plot.axes = {axis(1); axis(2);
##'                             points(design.grid[,1], design.grid[,2], pch = 21, bg = "white")
##'                             contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), 
##'                             matrix(obj.grid, n.grid), nlevels = 10,
##'                                    add=TRUE,drawlabels=TRUE, col = "black")
##'                             contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), 
##'                             matrix(cst.grid, n.grid), level = 0, add=TRUE,
##'                                    drawlabels=FALSE,lwd=1.5, col = "red")
##'                             }
##'               )
##' }
##' #---------------------------------------------------------------------------
##' # Expected Feasible Improvement surface with one inequality and one equality constraint
##' #---------------------------------------------------------------------------
##' \donttest{
##' set.seed(25468)
##' library(DiceDesign)
##' 
##' n_var <- 2 
##' fun.obj <- goldsteinprice
##' fun.cstineq <- function(x){return(3/2 - x[1] - 2*x[2] - .5*sin(2*pi*(x[1]^2 - 2*x[2])))}
##' fun.csteq <- function(x){return(branin(x) - 25)}
##' n.grid <- 51
##' test.grid <- expand.grid(X1 = seq(0, 1, length.out = n.grid), X2 = seq(0, 1, length.out = n.grid))
##' obj.grid <- apply(test.grid, 1, fun.obj)
##' cstineq.grid <- apply(test.grid, 1, fun.cstineq)
##' csteq.grid <- apply(test.grid, 1, fun.csteq)
##' n.init <- 25 
##' design.grid <- round(maximinESE_LHS(lhsDesign(n.init, n_var, seed = 42)$design)$design, 1)
##' obj.init <- apply(design.grid, 1, fun.obj)
##' cstineq.init <- apply(design.grid, 1, fun.cstineq)
##' csteq.init <- apply(design.grid, 1, fun.csteq)
##' model.fun <- km(~., design = design.grid, response = obj.init)
##' model.constraintineq <- km(~., design = design.grid, response = cstineq.init)
##' model.constrainteq <- km(~., design = design.grid, response = csteq.init)
##' 
##' models.cst <- list(model.constraintineq, model.constrainteq)
##'  
##' EFI_grid <- apply(test.grid, 1, crit_EFI, model.fun = model.fun, model.constraint = models.cst,
##'                   equality = c(FALSE, TRUE), critcontrol = list(tolConstraints = c(0.05, 3)))
##' 
##' filled.contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), nlevels = 50,
##'                matrix(EFI_grid, n.grid), main = "Expected Feasible Improvement",
##'                xlab = expression(x[1]), ylab = expression(x[2]), color = terrain.colors, 
##'                plot.axes = {axis(1); axis(2);
##'                             points(design.grid[,1], design.grid[,2], pch = 21, bg = "white")
##'                             contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), 
##'                             matrix(obj.grid, n.grid), nlevels = 10,
##'                                    add=TRUE,drawlabels=TRUE, col = "black")
##'                             contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), 
##'                             matrix(cstineq.grid, n.grid), level = 0, add=TRUE,
##'                                    drawlabels=FALSE,lwd=1.5, col = "red")
##'                             contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid),
##'                             matrix(csteq.grid, n.grid), level = 0, add=TRUE,
##'                                    drawlabels=FALSE,lwd=1.5, col = "orange")
##'                             }
##'               )
##' }

crit_EFI <- function(x, model.fun, model.constraint, equality = FALSE,
                     critcontrol = NULL, plugin = NULL, type = "UK"){
  n.cst <- length(model.constraint)
  
  if(n.cst == 1 && class(model.constraint) != "list"){
    model.constraint <- list(model.constraint)
  }
  
  # Regroup all observations
  observations <- c()
  for (i in 1:n.cst) observations <- cbind(observations, model.constraint[[i]]@y)
  
  
  if(is.null(plugin)){
    # Current feasible best
    feasibility <- test_feas_vec(cst=observations, equality=equality, tolConstraints = critcontrol$tolConstraints)
    if (any(feasibility)) plugin <- min(model.fun@y[feasibility])
    else                         plugin <- Inf ## When no feasible point is known 
  }
  
  if(n.cst == 1 && class(model.constraint) != "list"){
    model.constraint <- list(model.constraint)
  }
  
  ## For now obj and constraints are separated
  if(checkPredict(x, c(model.constraint, model.fun), type = type,
                  distance = critcontrol$distance, threshold = critcontrol$threshold
                  ))
    #|| checkPredict(x, model.constraint, type = type, distance = critcontrol$distance,
                              #      threshold = critcontrol$threshold))
    {
    return(0)
  }else{
    # Case when no feasible point is known
    if(plugin == Inf){
      Ei <- 1
    }else if(class(model.fun) == "km"){
      # Expected Improvement at x
      Ei <- EI(x = x, model = model.fun, plugin = plugin, type = type)
    }else{
      # Case of a fast-to-evaluate function
      x.data <- data.frame(t(as.numeric(x)))
      Ei <- -min(0, predict(object=model.fun, newdata=x.data, type=type, checkNames = FALSE)$mean - plugin)
    }
    
    if(length(equality) == 1 && equality == FALSE){
      equalityvec <- rep(FALSE, n.cst) 
    }else{
      equalityvec <- equality
    }
    
    if(is.null(critcontrol$tolConstraints)){
      tolvec <- rep(0, n.cst)
      if(any(equality != FALSE)){
        tolvec[equality] <- 0.05
        cat("No tolerance for equality constraints provided, 0.05 is used \n")
      }
    }else{
      tolvec <- critcontrol$tolConstraints
    }
    
    # Expected Feasibility
    if(class(model.constraint) == "km"){
      Pf <- Feas(model.constraint, newdata = t(x), type = "UK", equality = equalityvec,
                 threshold = tolvec)
    }else{
#       Pf <- prod(apply(model.constraint, Feas, newdata = t(x), type = "UK"))
      Pf <- prod(mapply(Feas, model = model.constraint, equality = equalityvec, threshold = tolvec,
                        MoreArgs = list(newdata = data.frame(t(x)), type = "UK")))
    }
    
    return(Pf*Ei) 
  } 
}

# Probability of feasibility for inequality or equality constraints
Feas <- function(model, newdata, type, equality = FALSE, threshold = 0){
  tmp <- predict(model, newdata = newdata, type = "UK", light.return = TRUE, checkNames = FALSE)
  if(equality == FALSE){
    return(pnorm(q = threshold, mean = tmp$mean, sd = tmp$sd))
  }else{
    return(pnorm(q = threshold, mean = tmp$mean, sd = tmp$sd) - pnorm(q = -threshold, mean = tmp$mean, sd = tmp$sd))
  }
}



