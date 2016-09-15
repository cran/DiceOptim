##' Computes the Stepwise Uncertainty Reduction (SUR) criterion at current location
##' 
##' @title Stepwise Uncertainty Reduction criterion
##' 
##' @param x a vector representing the input for which one wishes to calculate \code{SUR},
##' @param model.fun object of class \code{\link[DiceKriging]{km}} corresponding to the objective function,
##' or, if the objective function is fast-to-evaluate, a \code{\link[DOlab]{fastfun}} object,
##' @param model.constraint either one or a list of objects of class \code{\link[DiceKriging]{km}}, one for each constraint function,
##' @param equality either \code{FALSE} if all constraints are for inequalities, else a vector of boolean indicating which are equalities 
##' @param critcontrol optional list with arguments:
##' \itemize{
##' \item \code{tolConstraints} optional vector giving a tolerance (> 0) for each of the constraints (equality or inequality).
##' It is highly recommended to use it when there are equality constraints since the default tolerance of 0.05 in such case might not be suited;
##' \item \code{integration.points} and \code{integration.weights}: optional matrix and vector of integration points;
##' \item \code{precalc.data.cst, precalc.data.obj, mn.X.cst, sn.X.cst, mn.X.obj, sn.X.obj}: useful quantities for the 
##' fast evaluation of the criterion. 
##' \item Options for the \code{\link[DOlab]{checkPredict}} function: \code{threshold} (\code{1e-4}) and \code{distance} (\code{covdist}) 
##' are used to avoid numerical issues occuring when adding points too close to the existing ones.
##' }
##' @param type "\code{SK}" or "\code{UK}" (by default), depending whether uncertainty related to trend estimation 
##'        has to be taken into account. 
##' @return The Stepwise Uncertainty Reduction criterion at \code{x}.
##' @seealso \code{\link[DiceOptim]{EI}} from package DiceOptim, \code{\link[DOlab]{crit_EFI}}, \code{\link[DOlab]{crit_AL}}.
##' 
##' @export
##' @importFrom pbivnorm pbivnorm
##' 
##' @author
##' Victor Picheny 
##' 
##' Mickael Binois 
##' 
##' @references 
##' V. Picheny (2014),
##' A stepwise uncertainty reduction approach to constrained global optimization,
##' \emph{Proceedings of the 17th International Conference on Artificial Intelligence and Statistics}, JMLR W&CP 33, 787-795.
##' 
##' @examples
##' #---------------------------------------------------------------------------
##' # Stepwise Uncertainty Reduction criterion surface with one inequality constraint
##' #---------------------------------------------------------------------------
##' \donttest{
##' set.seed(25468)
##' library(DiceDesign)
##' 
##' n_var <- 2 
##' fun.obj <- goldsteinprice
##' fun.cst <- function(x){return(-branin(x) + 25)}
##' n.grid <- 21
##' test.grid <- expand.grid(X1 = seq(0, 1, length.out = n.grid), X2 = seq(0, 1, length.out = n.grid))
##' obj.grid <- apply(test.grid, 1, fun.obj)
##' cst.grid <- apply(test.grid, 1, fun.cst)
##' 
##' n_appr <- 15 
##' design.grid <- round(maximinESE_LHS(lhsDesign(n_appr, n_var, seed = 42)$design)$design, 1)
##' obj.init <- apply(design.grid, 1, fun.obj)
##' cst.init <- apply(design.grid, 1, fun.cst)
##' model.fun <- km(~., design = design.grid, response = obj.init)
##' model.constraint <- km(~., design = design.grid, response = cst.init)
##' 
##' integration.param <- integration_design_cst(integcontrol =list(integration.points = test.grid),
##'                                             lower = rep(0, n_var), upper = rep(1, n_var))
##' 
##' SUR_grid <- apply(test.grid, 1, crit_SUR_cst, model.fun = model.fun,
##'                   model.constraint = model.constraint, critcontrol=integration.param)
##' 
##' filled.contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), nlevels = 50,
##'                matrix(SUR_grid, n.grid), main = "SUR criterion",
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
##' # SUR with one inequality and one equality constraint
##' #---------------------------------------------------------------------------
##' \donttest{
##' set.seed(25468)
##' library(DiceDesign)
##' 
##' n_var <- 2 
##' fun.obj <- goldsteinprice
##' fun.cstineq <- function(x){return(3/2 - x[1] - 2*x[2] - .5*sin(2*pi*(x[1]^2 - 2*x[2])))}
##' fun.csteq <- function(x){return(branin(x) - 25)}
##' n.grid <- 21
##' test.grid <- expand.grid(X1 = seq(0, 1, length.out = n.grid), X2 = seq(0, 1, length.out = n.grid))
##' obj.grid <- apply(test.grid, 1, fun.obj)
##' cstineq.grid <- apply(test.grid, 1, fun.cstineq)
##' csteq.grid <- apply(test.grid, 1, fun.csteq)
##' n_appr <- 25 
##' design.grid <- round(maximinESE_LHS(lhsDesign(n_appr, n_var, seed = 42)$design)$design, 1)
##' obj.init <- apply(design.grid, 1, fun.obj)
##' cstineq.init <- apply(design.grid, 1, fun.cstineq)
##' csteq.init <- apply(design.grid, 1, fun.csteq)
##' model.fun <- km(~., design = design.grid, response = obj.init)
##' model.constraintineq <- km(~., design = design.grid, response = cstineq.init)
##' model.constrainteq <- km(~., design = design.grid, response = csteq.init)
##' 
##' models.cst <- list(model.constraintineq, model.constrainteq)
##'  
##' SUR_grid <- apply(test.grid, 1, crit_SUR_cst, model.fun = model.fun, model.constraint = models.cst,
##'                   equality = c(FALSE, TRUE), critcontrol = list(tolConstraints = c(0.05, 3), 
##'                   integration.points=integration.param$integration.points))
##' 
##' filled.contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), nlevels = 50,
##'                matrix(SUR_grid, n.grid), main = "SUR criterion",
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
##' 
crit_SUR_cst <- function(x, model.fun, model.constraint, equality = FALSE, critcontrol = NULL, type = "UK")
{
  n.cst <- length(model.constraint)
  
  if(n.cst > 3){
    cat("crit_SUR_cst does not take more than 3 constraints \n")
    return(NA)
  }
  
  if (is.null(dim(x))) x <- matrix(x, nrow=1)
  if(n.cst == 1 && class(model.constraint) != "list")  model.constraint <- list(model.constraint)
  
  # Check if x is not singular
  if(checkPredict(x, c(model.constraint, model.fun), type = type, distance = critcontrol$distance, threshold = critcontrol$threshold))
  {
    return(0)
  } else {
    
    #--- Complete optional inputs if missing --------------------------------
    if(is.null(critcontrol$integration.points)){
      integration.param   <- integration_design_cst(integcontrol=critcontrol, lower=critcontrol$lower, upper=critcontrol$upper,
                                                    model.fun=model.fun, model.constraint=model.constraint)
      integration.points  <- integration.param$integration.points
    } else {
      integration.points  <- critcontrol$integration.points
      integration.weights <- critcontrol$integration.weights
    }
    
    if (is.null(nrow(integration.points))) {integration.points <- as.matrix(integration.points)}
    n.integration.points <- nrow(integration.points)
    if (is.null(integration.weights)) {integration.weights <- rep(1/n.integration.points, n.integration.points)}
    
    if(is.null(critcontrol$precalc.data.cst) || is.null(critcontrol$mn.X.cst) || is.null(critcontrol$sn.X.cst) ||
         (is.null(critcontrol$precalc.data.obj)&&class(model.fun)=="km") || is.null(critcontrol$mn.X.obj) || is.null(critcontrol$sn.X.obj) ){
      precalc.data.cst <- vector("list", n.cst)
      mn.X.cst <- sn.X.cst <- matrix(0, n.cst, nrow(integration.points))
      
      for (i in 1:n.cst){
        p.tst <- predict(model.constraint[[i]], newdata=integration.points, type=type, checkNames=FALSE)
        mn.X.cst[i,] <- p.tst$mean
        sn.X.cst[i,]   <- p.tst$sd
        precalc.data.cst[[i]] <- precomputeUpdateData(model.constraint[[i]], integration.points)
      }
      
      p.tst <- predict(model.fun, newdata=integration.points, type=type, checkNames=FALSE)
      mn.X.obj <- p.tst$mean
      sn.X.obj   <- p.tst$sd
      if (class(model.fun)=="km") precalc.data.obj <- precomputeUpdateData(model.fun, integration.points)
      else                        precalc.data.obj <- NULL
      
    } else {
      mn.X.cst <- critcontrol$mn.X.cst
      sn.X.cst <- critcontrol$sn.X.cst
      precalc.data.cst <- critcontrol$precalc.data.cst
      mn.X.obj <- critcontrol$mn.X.obj
      sn.X.obj <- critcontrol$sn.X.obj
      precalc.data.obj <- critcontrol$precalc.data.obj
    }
    
    #--- Compute current best ----------------------------------------------------
    obs.cst <- c()
    for (i in 1:n.cst) obs.cst <- cbind(obs.cst, model.constraint[[i]]@y)
    feasibility <- test_feas_vec(cst=obs.cst, equality=equality, tolConstraints = critcontrol$tolConstraints)
    
    if (any(feasibility)) obj.min <- min(model.fun@y[feasibility])
    else                  obj.min <- Inf
#     print(obj.min)
#     observations <- c()    
#     for (i in 1:n.cst) observations <- cbind(observations, model.constraint[[i]]@y)
#     
#     if (sum(apply(observations, 1, max)<=0) > 0) obj.min <- min(model.fun@y[apply(observations, 1, max)<=0])
#     else                                         obj.min <- Inf
    #     else                                         obj.min <- (2*max(model.fun@y) + min(model.fun@y))/2
    #--- Remove integration points with zero variance -----------------------
    A <- matrix( rep(sqrt(sapply(model.constraint, function(model) return(model@covariance@sd2)))/1e4, n.integration.points), ncol=n.integration.points)
    I <- which(apply(sn.X.cst < A, 2, prod)!=0)
    
    if (length(I) > 0){
      integration.points <- integration.points[-I,,drop=FALSE]
      mn.X.obj  <- mn.X.obj[-I]
      sn.X.obj  <- sn.X.obj[-I]
      precalc.data.obj$Kinv.c.olddata <- precalc.data.obj$Kinv.c.olddata[,-I,drop=FALSE]
      precalc.data.obj$first.member   <- precalc.data.obj$first.member[-I]
      mn.X.cst  <- mn.X.cst[,-I, drop=FALSE]
      sn.X.cst  <- sn.X.cst[,-I, drop=FALSE]
      for (i in 1:n.cst) {
        precalc.data.cst[[i]]$Kinv.c.olddata <- precalc.data.cst[[i]]$Kinv.c.olddata[,-I,drop=FALSE]
        precalc.data.cst[[i]]$first.member   <- precalc.data.cst[[i]]$first.member[-I]
      }
      if (length(integration.weights) > 1){integration.weights <- integration.weights[-I]}
    }
    n.integration.points <- nrow(integration.points)
    X.new <- as.matrix(x)
    if (is.null(n.integration.points)) {
      cat("No valid integration point \n")
      return(0)
    }
    
    #--- Check equality setting ---------------------------------------------------
    if(is.null(critcontrol$tolConstraints)){
      tolvec <- rep(0, n.cst)
      if(any(equality != FALSE)){
        tolvec[equality] <- 0.05
      }
    } else {
      tolvec <- critcontrol$tolConstraints
    }
    
    #--- Update model.fun ---------------------------------------------------
    krig.obj <- predict(object=model.fun, newdata=x, type=type, se.compute=TRUE, cov.compute=FALSE, checkNames=FALSE)   
    mk.obj <- krig.obj$mean
    sk.obj <- krig.obj$sd
    if (class(model.fun)=="km") {
#       if (obj.min!=Inf) {
        kn.obj <- computeQuickKrigcov2(model.fun,integration.points=(integration.points),X.new=(X.new),
                                                precalc.data=precalc.data.obj, F.newdata=krig.obj$F.newdata, c.newdata=krig.obj$c)
        
        rho.obj <- kn.obj / (sk.obj*sn.X.obj)
        eta.obj <- (mk.obj - mn.X.obj) / sqrt( sk.obj^2 + sn.X.obj^2 - 2*kn.obj)
        nu.obj  <- (kn.obj - sk.obj^2) / sk.obj / sqrt( sk.obj^2 + sn.X.obj^2 - 2*kn.obj )
        obj.min.tilde <- (obj.min - mn.X.obj) / sn.X.obj
        p.obj.min.tilde <- pnorm(obj.min.tilde)
        obj.min.bar <- (obj.min - mk.obj)/sk.obj
        
        w1 <- w2 <- rep(0, n.integration.points)
        if (obj.min.bar > -10){
          if (obj.min.bar==Inf) {
            w1 <- pnorm(eta.obj)
          } else {
            J  <- which(eta.obj > -10)
            w1[J] <- pbivnorm(rep(obj.min.bar,  length(J)), eta.obj[J], nu.obj[J])
          }
        }
        if (-obj.min.bar > -10){
          if (-obj.min.bar==Inf) {
            w2 <- p.obj.min.tilde
          } else {
            J  <- which(obj.min.tilde > -10)
            w2[J] <- pbivnorm(rep(-obj.min.bar,  length(J)), obj.min.tilde[J], -rho.obj[J])
          }
        }
        pf_minus <- as.numeric(w1 + w2)
        pf_plus  <- as.numeric(p.obj.min.tilde)
#       } else {
#         pf_minus <- pf_plus <- p.obj.min.tilde <- rep(1, n.integration.points)
#       }
    } else {
      p.obj.min.tilde <- as.numeric(mn.X.obj <= obj.min)
      pf_minus <- as.numeric(mn.X.obj <= min(mk.obj, obj.min))
      pf_plus  <- as.numeric(mn.X.obj <= obj.min)
    }
    #--- Update model.constraint ----------------------------------------------------
    Tbar   <- Tbarm <- Tbarp <- rep(0, n.cst)
    Ttilde <- Ttildem <- Ttildep <- rho.cst <- matrix(0, n.cst, n.integration.points)
    
    for (i in 1:n.cst) {
      krig.cst  <- predict(object=model.constraint[[i]], newdata=x, type=type, se.compute=TRUE,cov.compute=FALSE,checkNames=FALSE) 
      mk.cst    <- krig.cst$mean
      sk.cst    <- krig.cst$sd
      kn.cst <- computeQuickKrigcov2(model.constraint[[i]],integration.points=integration.points,X.new=(X.new),
                                              precalc.data=precalc.data.cst[[i]], F.newdata=krig.cst$F.newdata, c.newdata=krig.cst$c)
      
      if (any(equality != FALSE)) {
        # Equality setting
        if (equality[i]) {
          Tbarp[i]     <- (tolvec[i] - mk.cst)/sk.cst
          Ttildep[i,]  <- (tolvec[i] - mn.X.cst[i,] )/sn.X.cst[i,]
          Tbarm[i]     <- (-tolvec[i] - mk.cst)/sk.cst
          Ttildem[i,]  <- (-tolvec[i] - mn.X.cst[i,] )/sn.X.cst[i,]
        } else {
          Tbar[i]     <- (-tolvec[i] - mk.cst)/sk.cst
          Ttilde[i,]  <- (-tolvec[i] - mn.X.cst[i,] )/sn.X.cst[i,]
        }
      } else {
        Tbar[i]     <- (-tolvec[i] - mk.cst)/sk.cst
        Ttilde[i,]  <- (-tolvec[i] - mn.X.cst[i,] )/sn.X.cst[i,] 
      }
      
      rho.cst[i,] <- kn.cst / (sk.cst*sn.X.cst[i,] )
    }
    
    #--- Compute criterion -------------------------------------------------
    
    if (any(equality != FALSE)) {
      p.f <- 1
      for (i in 1:n.cst) {
        if (equality[i]) {
          p.f <- p.f*(pnorm(Ttildep[i,]) - pnorm(Ttildem[i,]))
        } else {
          p.f <- p.f*pnorm(Ttilde[i,])
        }
      }
      oldcrit <- p.obj.min.tilde*p.f
    } else {
      oldcrit <- p.obj.min.tilde*apply(pnorm(Ttilde),2,prod)
    }
    
    pg_minus <- pg_plus <- matrix(0, n.cst, n.integration.points)
    
    for (i in 1:n.cst) {
      
      if (any(equality != FALSE)) {
        #-- Equality setting ------------------------
        if (equality[i]) {
          pg_minus[i,] <- pbivnorm(rep(Tbarp[i],  n.integration.points), Ttildep[i,], rho.cst[i,]) -
            pbivnorm(rep(Tbarp[i],  n.integration.points), Ttildem[i,], rho.cst[i,]) -
            pbivnorm(rep(Tbarm[i],  n.integration.points), Ttildep[i,], rho.cst[i,]) +
            pbivnorm(rep(Tbarm[i],  n.integration.points), Ttildem[i,], rho.cst[i,])
          
          pg_plus[i,] <- pnorm(Ttildep[i,]) - pnorm(Ttildem[i,]) - pg_minus[i,]
          
        } else {
          if (Tbar[i] > -10) {
            J <- which(Ttilde[i,] > -10)
            pg_minus[i,J] <- pbivnorm(rep(Tbar[i],  length(J)), Ttilde[i,J], rho.cst[i,J])
          }
          pg_plus[i,]  <- pnorm(Ttilde[i,]) - pg_minus[i,]
        }
      } else {
        #-- Inequality setting ----------------------
        if (Tbar[i] > -10) {
          J <- which(Ttilde[i,] > -10)
          pg_minus[i,J] <- pbivnorm(rep(Tbar[i],  length(J)), Ttilde[i,J], rho.cst[i,J])
        }
        #       pg_minus[i,] <- pbivnorm(rep(Tbar[i],  n.integration.points), Ttilde[i,], rho.cst[i,])
        pg_plus[i,]  <- pnorm(Ttilde[i,]) - pg_minus[i,]
      }
    }
    
    
    if (n.cst==1){
      newcrit <- pf_minus*pg_minus + pf_plus*pg_plus
    } 
    if (n.cst==2){
      newcrit <- pf_minus*apply(pg_minus, 2, prod) + pf_plus*(pnorm(Ttilde[1,])*pg_plus[2,] + 
                                                                pnorm(Ttilde[2,])*pg_plus[1,] - apply(pg_plus, 2, prod))
    }
    if (n.cst==3){
      newcrit <- pf_minus*apply(pg_minus, 2, prod) + 
        pf_plus*(pnorm(Ttilde[2,])*pnorm(Ttilde[3,])*pg_plus[1,] + 
                   pnorm(Ttilde[1,])*pnorm(Ttilde[3,])*pg_plus[2,] +   
                   pnorm(Ttilde[1,])*pnorm(Ttilde[2,])*pg_plus[3,] -
                   pnorm(Ttilde[3,])*pg_plus[1,]*pg_plus[2,] -
                   pnorm(Ttilde[2,])*pg_plus[1,]*pg_plus[3,] -
                   pnorm(Ttilde[1,])*pg_plus[2,]*pg_plus[3,] +
                   apply(pg_plus, 2, prod))
    }
    a <- oldcrit - newcrit
    crit <- mean(a*integration.weights, na.rm = TRUE)
    
    return(crit)
  }
}