##' Computes the Expected Augmented Lagrangian Improvement at current location, with our without slack variables. Depending on the cases,
##' the computation is either analytical (very fast), based on MC integration (slow) or on the CDF of a weighted sum of non-central
##' chi-square (WNCS) variates (intermediate)
##' 
##' @title Expected Augmented Lagrangian Improvement
##' 
##' @param x either a vector representing the design or the design AND slack variables (see details)
##' @param model.fun object of class \code{\link[DiceKriging]{km}} correspostnding to the objective function,
##' or, if the objective function is fast-to-evaluate, a \code{\link[DOlab]{fastfun}} object,
##' @param model.constraint either one or a list of objects of class \code{\link[DiceKriging]{km}}, one for each constraint function,
##' @param equality either \code{FALSE} if all constraints are for inequalities, or a vector of Booleans indicating which are equalities 
##' @param critcontrol optional list with the following arguments:
##' \itemize{
##'  \item{\code{slack}}{: logical. If TRUE, slack variables are used for inequality constraints (see Details)}
##'  \item{\code{rho}}{: penalty term (scalar),}
##'  \item{\code{lambda}}{: Lagrange multipliers (vector of size the number of constraints),}
##'  \item{\code{elit}}{: logical. If TRUE, sets the criterion to zero for all x's not improving the objective function}
##'  \item{\code{n.mc}}{: number of Monte-Carlo drawings used to evaluate the criterion (see Details)}
##'  \item{\code{nt}}{: number of discretization points for the WNCS distribution (see Details)}
##'  \item \code{tolConstraints}, an optional vector giving a tolerance (> 0) for each of the constraints (equality or inequality).
##' It is highly recommended to use it when there are equality constraints since the default tolerance of 0.05 in such case might not be suited.
##'  }
##'  
##' Options for the \code{\link[DOlab]{checkPredict}} function: \code{threshold} (\code{1e-4}) and \code{distance} (\code{covdist}) 
##' are used to avoid numerical issues occuring when adding points too close to the existing ones.
##' @param type "\code{SK}" or "\code{UK}" (by default), depending whether uncertainty related to trend estimation 
##'        has to be taken into account. 
##' @return The Expected Augmented Lagrangian Improvement at \code{x}.
##' @seealso \code{\link[DiceOptim]{EI}} from package DiceOptim, \code{\link[DOlab]{crit_EFI}}, \code{\link[DOlab]{crit_SUR_cst}}.
##' 
##' @details
##' The AL can be used with or without the help of slack variables for the inequality constraints. 
##' If \code{critcontrol$slack=FALSE}:
##' With a single constraint (inequality or equality) and a fast objective, a very fast formula is
##' used to compute the criterion (recommended setting). 
##' Otherwise, an MC estimator of the criterion is used, which is much more costly. The argument 
##' \code{critcontrol$n.mc} tunes the precision of the estimator.
##' On both cases \code{x} must be of size \code{d}.
##' 
##' If \code{critcontrol$slack=TRUE}:
##' Slack variables are used to handle the inequality constraints. 
##' They can be provided directly through \code{x}, which should be of size \code{d+} the number of inequality constraints.
##' The last values of \code{x} are slack variables scaled to [0,1].
##' 
##' If \code{x} is of size \code{d}, estimates of optimal slack variable are used.\cr
##' 
##' @export
##' 
##' @author
##' Victor Picheny 
##' 
##' Mickael Binois 
##' 
##' @references 
##' R.B. Gramacy, G.A. Gray, S. Le Digabel, H.K.H Lee, P. Ranjan, G. Wells, Garth, and S.M. Wild (2014+),
##' Modeling an augmented Lagrangian for improved blackbox constrained optimization,
##' \emph{arXiv preprint arXiv:1403.4890}.
##' 
##' @examples
##' #---------------------------------------------------------------------------
##' # Expected Augmented Lagrangian Improvement surface with one inequality constraint,
##' # fast objective
##' #---------------------------------------------------------------------------
##' \donttest{
##' set.seed(25468)
##' library(DiceDesign)
##' 
##' n_var <- 2 
##' fun.obj <- goldsteinprice
##' fun.cst <- function(x){return(-branin(x) + 25)}
##' n.grid <- 31
##' test.grid <- expand.grid(X1 = seq(0, 1, length.out = n.grid), X2 = seq(0, 1, length.out = n.grid))
##' obj.grid <- apply(test.grid, 1, fun.obj)
##' cst.grid <- apply(test.grid, 1, fun.cst)
##' n.init <- 15 
##' design.grid <- round(maximinESE_LHS(lhsDesign(n.init, n_var, seed = 42)$design)$design, 1)
##' obj.init <- apply(design.grid, 1, fun.obj)
##' cst.init <- apply(design.grid, 1, fun.cst)
##' model.constraint <- km(~., design = design.grid, response = cst.init)
##' model.fun <- fastfun(fun.obj, design.grid)
##' AL_grid <- apply(test.grid, 1, crit_AL, model.fun = model.fun,
##'                   model.constraint = model.constraint)
##' 
##' filled.contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), nlevels = 50,
##'                matrix(AL_grid, n.grid), main = "Expected AL Improvement",
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
##' # Expected AL Improvement surface with one inequality and one equality constraint,
##' # using slack variables
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
##' AL_grid <- apply(test.grid, 1, crit_AL, model.fun = model.fun, model.constraint = models.cst,
##'                   equality = c(FALSE, TRUE), critcontrol = list(tolConstraints = c(0.05, 3),
##'                   slack=TRUE))
##' 
##' filled.contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), nlevels = 50,
##'                matrix(AL_grid, n.grid), main = "Expected AL Improvement",
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

crit_AL <- function(x, model.fun, model.constraint, equality = FALSE, critcontrol = NULL, type = "UK"){
  
  n.cst  <- length(model.constraint)
  if (n.cst!=length(equality)) equality <- rep(equality, n.cst)
  
  n.ineq <- sum(!equality)
  
  if (is.null(dim(x))) x <- matrix(x, nrow=1)
  
  # Retrieve [scaled] slack variables if provided
  slack.x <- NULL
  if (ncol(x) > model.fun@d) {
    slack.x <- rep(0, n.cst)
    slack.x[!equality] <- as.numeric(x[,model.fun@d+(1:n.ineq),drop=FALSE])
    x <- x[,1:model.fun@d,drop=FALSE]
  }
  
  if(n.cst == 1 && class(model.constraint) != "list"){
    model.constraint <- list(model.constraint)
  }
  
  if (is.null(critcontrol$lambda))      critcontrol$lambda <- matrix(rep(1,n.cst), nrow=1)
  if (is.null(critcontrol$rho))         critcontrol$rho <- .5
  if (is.null(critcontrol$elit))        critcontrol$elit <- FALSE
  if (!is.null(critcontrol$slack))      if (critcontrol$slack!=TRUE) critcontrol$slack <- NULL
  if (!is.null(critcontrol$proxy))      if (critcontrol$proxy!=TRUE) critcontrol$proxy <- NULL
  lambda <- matrix(critcontrol$lambda, nrow=1)
  rho    <- critcontrol$rho
  mu     <- 1/2/rho
  
  # Regroup all observations
  observations <- c()
  for (i in 1:n.cst) observations <- cbind(observations, model.constraint[[i]]@y)
  
  # For AL computation
  Cst.threshold <- rep(0, n.cst)
  Cst.threshold[equality] <- -Inf
  
  # Predict objective
  pred.obj <- predict(object = model.fun, newdata = x, type = type, checkNames = FALSE)
  
  # Predict constraint
  pred.cst.mean <- pred.cst.sd <- rep(0, n.cst)
  for (i in 1:n.cst) {
    pred.cst <- predict(object = model.constraint[[i]], newdata = x, type = type, checkNames = FALSE)
    pred.cst.mean[i] <- pred.cst$mean
    pred.cst.sd[i] <- max(pred.cst$sd, 1e-16)
  }
  
  # Compute DoE slack variables if option is active
  if (!is.null(critcontrol$slack)){
    slack.doe <- pmax( -matrix(rep(as.numeric(lambda)/2/mu, model.fun@n), ncol=n.cst, byrow=TRUE) - observations, 0)
    slack.doe[,equality] <- 0
    Cst.threshold <- rep(-Inf, n.cst)
    
    # Rescale slack.x if provided, otherwise replace with best estimate
    if (!is.null(slack.x)) {
      slack.x <- slack.x*pmax((- pred.cst.mean + 3*pred.cst.sd - critcontrol$lambda/2/mu), 0)
    } else {
      slack.x <- pmax((- pred.cst.mean - critcontrol$lambda/2/mu), 0)
      slack.x[equality] <- 0
    }
  } else {
    slack.doe <- matrix(0, nrow=model.fun@n, ncol=n.cst)
    slack.x   <- rep(0, n.cst)
  }
  
  # Test feasibility of observations
  feasibility <- test_feas_vec(cst=observations, equality=equality, tolConstraints = critcontrol$tolConstraints)
  if (any(feasibility)) fun.feas <- min(model.fun@y[feasibility])
  else                  fun.feas <- Inf
  
  # Get current ALmin
  almin <- min(model.fun@y + (observations + slack.doe)%*%t(lambda) + 
                 1/2/rho*rowSums(pmax((observations + slack.doe), t(matrix(rep(Cst.threshold,model.fun@n), ncol=model.fun@n)))^2))
  #================================================================================
  # CASE 1: x is singular
  #================================================================================
  if(checkPredict(x, c(model.constraint, model.fun), type = type,
                  distance = critcontrol$distance, threshold = critcontrol$threshold))
  {
    if (is.null(critcontrol$slack) || !is.null(critcontrol$proxy)) {
      # Return 0, unless elit case: then penalization
      eiy <- 0
      if (critcontrol$elit) eiy <- min(0, fun.feas - pred.obj$mean)
      if (!is.null(critcontrol$proxy)) eiy <- min(eiy, almin - (min(pred.obj$mean + (pred.cst.mean+slack.x)%*%t(lambda) + 
                                       1/2/rho*sum(pmax(pred.cst.mean+slack.x, Cst.threshold)^2))))
      return(eiy)
    } else {
      # Special penalization for slack - regular criterion only
      alpha <- lambda/mu + 2*slack.x
      # almin <- min(model.fun@y + (observations+slack.doe)%*%t(lambda) + 1/2/rho*rowSums((observations+slack.doe)^2))
      wmin  <-  (almin - pred.obj$mean - mu*sum(slack.x^2) - slack.x%*%t(lambda) + mu*sum(alpha^2/4))
      
      return(min(wmin, 0))
    }
  } else {
    #================================================================================
    ## CASE 2: ELITIST CASE: REPLACE BY PENALTY IF NOT IMPROVING
    #================================================================================
    if (critcontrol$elit &&  pred.obj$mean>fun.feas && class(model.fun) != "km" && sum(equality)<1){
      eiy <- fun.feas - pred.obj$mean
      if (!is.null(critcontrol$proxy)) {
        eiy <- min(eiy, almin - (min(pred.obj$mean + (pred.cst.mean+slack.x)%*%t(lambda) + 
                                       1/2/rho*sum(pmax(pred.cst.mean+slack.x, Cst.threshold)^2))))
      }
    } else {
      #================================================================================
      # CASE 3: PROXY CRITERION (expected value instead of EI)
      #================================================================================
      if (!is.null(critcontrol$proxy)) {
        #         Cst.threshold <- rep(0, n.cst)
        #         Cst.threshold[equality] <- -Inf
        
        # Get current best AL
        # almin <- min(model.fun@y + (observations+slack.doe)%*%t(lambda) + 1/2/rho*rowSums((observations+slack.doe)^2))
        
        # if (is.null(critcontrol$slack)) {
        #--- Subcase 3.1: no slack --------------------------------------------------
        eiy <- almin - (min(pred.obj$mean + (pred.cst.mean+slack.x)%*%t(lambda) + 
                              1/2/rho*sum(pmax(pred.cst.mean+slack.x, Cst.threshold)^2)))
        # } else {
        #--- Subcase 3.2: with slack ------------------------------------------------
        #           
        #           eiy <- almin - min(pred.obj$mean + (pred.cst.mean+slack.x)%*%t(lambda) + 
        #                         1/2/rho*sum((pred.cst.mean+slack.x)^2))
        #         }
      } else {
        #================================================================================
        # CASE 4: NO SLACK
        #================================================================================
        if (is.null(critcontrol$slack)) {
          #--- Subcase 4.1: single constraint, fast objective (analytical) ----------------
          if (n.cst==1 && class(model.fun)!="km") {
            
            lambda <- as.numeric(lambda)
            model.constraint <- model.constraint[[1]]
            #             if (equality) almin <- min( model.fun@y + model.constraint@y*lambda + (1/(2*rho)) * (model.constraint@y)^2 )
            #             else          almin <- min( model.fun@y + model.constraint@y*lambda + (1/(2*rho)) * pmax(model.constraint@y,0)^2 )
            
            # Precompute quantities
            T  <- 2*rho*(almin - pred.obj$mean)
            alpha <- 2*rho*lambda
            T2 <- T + alpha^2/4
            mZ <- pred.cst.mean
            sZ <- pred.cst.sd
            
            if (equality) {
              eiy <- 0
              if (T2 > 0){
                b  <- (sqrt(T2) - alpha/2 - mZ)/sZ
                a  <- (-sqrt(T2) - alpha/2 - mZ)/sZ
                eiy <- T2*(pnorm(b) - pnorm(a)) - 
                  ((alpha/2 + mZ)^2*(pnorm(b) - pnorm(a)) + (alpha+2*mZ)*sZ*(dnorm(a)-dnorm(b)) + 
                     sZ^2*(a*dnorm(a)-b*dnorm(b) + pnorm(b) - pnorm(a)))
              }
            } else {
              u  <- (pmin(T/alpha,0)-mZ)/sZ
              eiy <- T*pnorm(u) - alpha*(mZ*pnorm(u) - sZ*dnorm(u))
              if (T > 0){
                v  <- (sqrt(T2) - alpha/2 - mZ)/sZ
                eiy <- eiy + (T2 - (alpha/2 + mZ)^2-sZ^2)*(pnorm(v) - pnorm(-mZ/sZ)) +
                  (alpha + 2*mZ)*sZ*(dnorm(v) - dnorm(mZ/sZ)) + 
                  sZ^2*( (v)*dnorm(v) + (mZ/sZ)*dnorm(mZ/sZ) )
              }
            }
            eiy <- eiy/(2*rho)
          } else {
            #--- Subcase 4.2: otherwise, MC evaluation ---------------------------------
            if (is.null(critcontrol$n.mc)) critcontrol$n.mc <- 100
            n.mc <- critcontrol$n.mc
            cst.mc <- matrix(0, n.cst, n.mc)
            
            Cst.threshold <- rep(0, n.cst)
            Cst.threshold[equality] <- -Inf
            
            for (i in 1:n.cst) cst.mc[i,]  <- rnorm(n=n.mc, mean=pred.cst.mean[i], sd=pred.cst.sd[i])
            obj.mc   <- rnorm(n=n.mc, mean=pred.obj$mean, sd=pred.obj$sd)
            #             
            #             almin <- min(model.fun@y + observations%*%t(lambda) + 
            #                            1/2/rho*rowSums(pmax(observations, t(matrix(rep(Cst.threshold,model.fun@n), ncol=model.fun@n)))^2))
            
            # Compute MC estimate
            eiy <- mean(pmax(almin - obj.mc - t(cst.mc)%*%t(lambda) - 
                               1/2/rho * colSums(pmax(cst.mc, t(matrix(rep(Cst.threshold,n.mc), ncol=n.mc)))^2), 0))
          }
        } else {
          #===========================================================================
          # CASE 5: WITH SLACK
          #===========================================================================
          
          #           # Rescale slack.x if provided, otherwise replace with best estimate
          #           if (!is.null(slack.x)) {
          #             slack.x <- slack.x*pmax((- pred.cst.mean + 3*pred.cst.sd - critcontrol$lambda/2/mu), 0)
          #           } else {
          #             slack.x <- pmax((- pred.cst.mean - critcontrol$lambda/2/mu), 0)
          #             slack.x[equality] <- 0
          #           }
          
          # Get current best AL
          # almin <- min(model.fun@y + (observations+slack.doe)%*%t(lambda) + 1/2/rho*rowSums((observations+slack.doe)^2))
          
          # Parameters for chi2
          alpha <- lambda/mu + 2*slack.x
          mZ <- pred.cst.mean + alpha/2
          sZ <- pmax(pred.cst.sd, 1e-16)
          if (is.null(critcontrol$nt)) critcontrol$nt  <- 20
          nt <- critcontrol$nt
          eiy <- 0
          if (class(model.fun)!="km") {
            #--- Subcase 5.1: known objective -----------------------------------
            wmin  <-  1/ mu * (almin - pred.obj$mean - mu*sum(slack.x^2) - slack.x%*%t(lambda) + mu*sum(alpha^2/4))
            
            if (wmin > 0) {
              wts <- sZ^2
              ncp <- (mZ/sZ)^2
              t <- seq(0, wmin, length.out=nt)
              #               print(c(wmin, wts, ncp, mu))
              #               print(MYdavies(q=t, lambda=wts, delta = ncp)$Qq)
              eiy <- mu*sum(MYdavies(q=t, lambda=wts, delta = ncp)$Qq)*wmin/(nt-1)
            } else {
              eiy <- wmin
            }
          } else {
            #--- Subcase 5.2: unknown objective ----------------------------
            mF <- pred.obj$mean
            sF <- max(1e-16, pred.obj$sd)
            wmint <- 1/mu*(almin - mu*sum(slack.x^2) - slack.x%*%t(lambda) + mu*sum(alpha^2/4)  )    
            if (wmint > -3*sF/mu) {
              wts <- sZ^2
              ncp <- (mZ/sZ)^2
              t   <- seq(-3*sF/mu, wmint, length.out=nt)
              eiy <- mu*sum(MYdavies(q= t - mF/mu, lambda=wts, delta = ncp, sigma=sF/mu)$Qq)*(wmint+3*sF/mu)/(nt-1)
            } else {
              eiy <- wmint + 3*sF/mu
            }
          }
        }
        #===========================================================================
        if (is.na(eiy)) eiy <- -10
      }
    }
    return(eiy)
  }
}