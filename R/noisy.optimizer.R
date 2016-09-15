##' Optimization of homogenously noisy functions based on Kriging
##' 
##' Sequential optimization of kriging-based criterion conditional on noisy
##' observations, with model update after each evaluation.  Eight criteria are
##' proposed to choose the next observation: random search, sequential
##' parameter optimization (SPO), reinterpolation, Expected Improvement (EI)
##' with plugin, Expected Quantile Improvement (EQI), quantile minimization,
##' Augmented Expected Improvement (AEI) and Approximate Knowledge Gradient
##' (AKG). The criterion optimization is based on the package rgenoud.
##' 
##' 
##' @param optim.crit String defining the criterion to be optimized at each
##' iteration. Possible values are: "random.search", "SPO", "reinterpolation",
##' "EI.plugin", "EQI", "min.quantile", "AEI", "AKG".
##' @param optim.param List of parameters for the chosen criterion.  For
##' "EI.plugin": optim.param$plugin.type is a string defining which plugin is
##' to be used. Possible values are "ytilde", "quantile" and "other".  If
##' "quantile" is chosen, optim.param$quantile defines the quantile level.  If
##' "other" is chosen, optim.param$plugin directly sets the plugin value.
##' 
##' For "EQI": optim.param$quantile defines the quantile level. If not
##' provided, default value is 0.9.
##' 
##' For "min.quantile": optim.param$quantile defines the quantile level. If not
##' provided, default value is 0.1.
##' 
##' For "AEI": optim.param$quantile defines the quantile level to choose the
##' current best point. If not provided, default value is 0.75.
##' @param model a Kriging model of "km" class
##' @param n.ite Number of iterations
##' @param noise.var Noise variance (scalar). If noiseReEstimate=TRUE, it is an
##' initial guess for the unknown variance (used in optimization).
##' @param funnoise objective (noisy) function
##' @param lower vector containing the lower bounds of the variables to be
##' optimized over
##' @param upper vector containing the upper bounds of the variables to be
##' optimized over
##' @param parinit optional vector of initial values for the variables to be
##' optimized over
##' @param control optional list of control parameters for optimization.  One
##' can control \code{"pop.size"} (default : [N=3*2^dim for dim<6 and N=32*dim
##' otherwise]]), \code{"max.generations"} (N), \code{"wait.generations"} (2)
##' and \code{"BFGSburnin"} (0) of function \code{"genoud"} (see
##' \code{\link[rgenoud]{genoud}}).  Numbers into brackets are the default
##' values
##' @param CovReEstimate optional boolean specfiying if the covariance
##' parameters should be re-estimated at every iteration (default value = TRUE)
##' @param NoiseReEstimate optional boolean specfiying if the noise variance
##' should be re-estimated at every iteration (default value = FALSE)
##' @param nugget.LB optional scalar of minimal value for the estimated noise
##' variance. Default value is 1e-5.
##' @param estim.model optional kriging model of "km" class with homogeneous
##' nugget effect (no noise.var). Required if noise variance is reestimated and
##' the initial "model" has heterogenenous noise variances.
##' @param type "SK" or "UK" for Kriging with known or estimated trend
##' @return A list with components: \item{model}{the current (last) kriging
##' model of "km" class} \item{best.x}{ The best design found}
##' \item{best.y}{The objective function value at best.x} \item{best.index}{The
##' index of best.x in the design of experiments}
##' 
##' \item{history.x}{ The added observations} \item{history.y}{The added
##' observation values} \item{history.hyperparam}{The history of the kriging
##' parameters}
##' 
##' \item{estim.model}{If noiseReEstimate=TRUE, the current (last) kriging
##' model of "km" class for estimating the noise variance.}
##' \item{history.noise.var}{If noiseReEstimate=TRUE, the history of the noise
##' variance estimate.}
##' @author Victor Picheny 
##' 
##' @references V. Picheny and D. Ginsbourger (2013), Noisy kriging-based optimization
##' methods: A unified implementation within the DiceOptim package,
##' \emph{Computational Statistics & Data Analysis}
##' @examples
##' 
##' ##########################################################################
##' ### EXAMPLE 1: 3 OPTIMIZATION STEPS USING EQI WITH KNOWN NOISE         ###
##' ### AND KNOWN COVARIANCE PARAMETERS FOR THE BRANIN FUNCTION            ###
##' ##########################################################################
##' 
##' set.seed(10)
##' library(DiceDesign)
##' # Set test problem parameters
##' doe.size <- 9
##' dim <- 2
##' test.function <- get("branin2")
##' lower <- rep(0,1,dim)
##' upper <- rep(1,1,dim)
##' noise.var <- 0.1
##' 
##' # Build noisy simulator
##' funnoise <- function(x)
##' {     f.new <- test.function(x) + sqrt(noise.var)*rnorm(n=1)
##'       return(f.new)}
##' 
##' # Generate DOE and response
##' doe <- as.data.frame(lhsDesign(doe.size, dim)$design)
##' y.tilde <- funnoise(doe)
##' 
##' # Create kriging model
##' model <- km(y~1, design=doe, response=data.frame(y=y.tilde),
##'      covtype="gauss", noise.var=rep(noise.var,1,doe.size), 
##'      lower=rep(.1,dim), upper=rep(1,dim), control=list(trace=FALSE))
##' 
##' # Optimisation with noisy.optimizer (n.ite can be increased)
##' n.ite <- 2 
##' optim.param <- list()
##' optim.param$quantile <- .9
##' optim.result <- noisy.optimizer(optim.crit="EQI", optim.param=optim.param, model=model,
##' 		n.ite=n.ite, noise.var=noise.var, funnoise=funnoise, lower=lower, upper=upper,
##' 		NoiseReEstimate=FALSE, CovReEstimate=FALSE)
##' 
##' new.model <- optim.result$model
##' best.x    <- optim.result$best.x
##' new.doe   <- optim.result$history.x
##' 
##' \dontrun{
##' ##### DRAW RESULTS #####
##' # Compute actual function on a grid
##' n.grid <- 12
##' x.grid <- y.grid <- seq(0,1,length=n.grid)
##' design.grid <- expand.grid(x.grid, y.grid)
##' names(design.grid) <- c("V1","V2")
##' nt <- nrow(design.grid)
##' func.grid <- rep(0,1,nt)
##' 
##' for (i in 1:nt)
##' { func.grid[i] <- test.function(x=design.grid[i,])}
##' 
##' # Compute initial and final kriging on a grid
##' pred <- predict(object=model, newdata=design.grid, type="UK", checkNames = FALSE)
##' mk.grid1 <- pred$m
##' sk.grid1 <- pred$sd
##' 
##' pred <- predict(object=new.model, newdata=design.grid, type="UK", checkNames = FALSE)
##' mk.grid2 <- pred$m
##' sk.grid2 <- pred$sd
##' 
##' # Plot initial kriging mean
##' z.grid <- matrix(mk.grid1, n.grid, n.grid)
##' filled.contour(x.grid,y.grid, z.grid, nlevels=50, color = topo.colors,
##' plot.axes = {title("Initial kriging mean");
##' points(model@@X[,1],model@@X[,2],pch=17,col="black"); 
##' axis(1); axis(2)})
##' 
##' # Plot initial kriging variance
##' z.grid <- matrix(sk.grid1^2, n.grid, n.grid)
##' filled.contour(x.grid,y.grid, z.grid, nlevels=50, color = topo.colors,
##' plot.axes = {title("Initial kriging variance");
##' points(model@@X[,1],model@@X[,2],pch=17,col="black"); 
##' axis(1); axis(2)})
##' 
##' # Plot final kriging mean
##' z.grid <- matrix(mk.grid2, n.grid, n.grid)
##' filled.contour(x.grid,y.grid, z.grid, nlevels=50, color = topo.colors,
##' plot.axes = {title("Final kriging mean");
##' points(new.model@@X[,1],new.model@@X[,2],pch=17,col="black"); 
##' axis(1); axis(2)})
##' 
##' # Plot final kriging variance
##' z.grid <- matrix(sk.grid2^2, n.grid, n.grid)
##' filled.contour(x.grid,y.grid, z.grid, nlevels=50, color = topo.colors,
##' plot.axes = {title("Final kriging variance");
##' points(new.model@@X[,1],new.model@@X[,2],pch=17,col="black"); 
##' axis(1); axis(2)})
##' 
##' # Plot actual function and observations
##' z.grid <- matrix(func.grid, n.grid, n.grid)
##' tit <- "Actual function - Black: initial points; red: added points"
##' filled.contour(x.grid,y.grid, z.grid, nlevels=50, color = topo.colors,
##' plot.axes = {title(tit);points(model@@X[,1],model@@X[,2],pch=17,col="black"); 
##' points(new.doe[1,],new.doe[2,],pch=15,col="red");
##' axis(1); axis(2)})
##' }
##' 
##' ##########################################################################
##' ### EXAMPLE 2: 3 OPTIMIZATION STEPS USING EQI WITH UNKNOWN NOISE       ###
##' ### AND UNKNOWN COVARIANCE PARAMETERS FOR THE BRANIN FUNCTION          ###
##' ##########################################################################
##' # Same initial model and parameters as for example 1
##' n.ite <- 2 # May be changed to a larger value 
##' res <- noisy.optimizer(optim.crit="min.quantile", 
##' optim.param=list(type="quantile",quantile=0.01),
##' model=model, n.ite=n.ite, noise.var=noise.var, funnoise=funnoise, 
##' lower=lower, upper=upper, 
##' control=list(print.level=0),CovReEstimate=TRUE, NoiseReEstimate=TRUE)
##' 
##' # Plot actual function and observations
##' plot(model@@X[,1], model@@X[,2], pch=17,xlim=c(0,1),ylim=c(0,1))
##' points(res$history.x[1,], res$history.x[2,], col="blue")
##' 
##' # Restart: requires the output estim.model of the previous run
##' # to deal with potential repetitions
##' res2 <- noisy.optimizer(optim.crit="min.quantile", 
##' optim.param=list(type="quantile",quantile=0.01), 
##' model=res$model, n.ite=n.ite, noise.var=noise.var, funnoise=funnoise, 
##' lower=lower, upper=upper, estim.model=res$estim.model,
##' control=list(print.level=0),CovReEstimate=TRUE, NoiseReEstimate=TRUE)
##' 
##' # Plot new observations
##' points(res2$history.x[1,], res2$history.x[2,], col="red")
##' 
##' 
##' @export
noisy.optimizer <- function(optim.crit, optim.param=NULL, model, n.ite, noise.var=NULL, funnoise, lower, upper, 
                      parinit=NULL, control=NULL,CovReEstimate=TRUE,NoiseReEstimate=FALSE, 
                      nugget.LB=1e-5, estim.model=NULL, type="UK")
{
  ############################################################################################################
  # A switch is made at every iteration to choose the corresponding infill criterion, then the model is updated 
  # similarly for all methods.
  #
  # All methods return a final noisy kriging model. In addition, the quantities "all.X", "all.y", "all.thetas"
  # (size=n.ite) are saved for analysis.
  #
  # This version allows sequential re-estimation of the noise variance. An additional model is created specifically
  # for parameter estimation (estim.model), with uniform nugget (repeated experiments are not aggregated).
  # "model" and "estim.model" are equivalent, but "model" has less observations ( hence goes faster).
  #
  # Additional (optional) arguments:
  #    - "NoiseReEstimate": TRUE/FALSE
  #    - "noise.var": initial guess for the unknown variance (used in optimization)
  #    - "obs.n.rep": number of repetitions per observation point. Required if "model" has heterogeneous variances
  #    - "estim.model": model with homogeneous nugget effect (no noise.var). Required only for restarting algorithms
  #
  ############################################################################################################
  
  cat("Starting noisy optimization with the following criterion and parameters \n")
  cat(optim.crit, "\n")
  if(!is.null(optim.param))  print(optim.param,quote=FALSE)
  cat("----------------- \n")
  optim.result <- list()
  
  #---------------------------------------------------------------------------------------------------------
  # Initialization for the unknown noise variance case
  #---------------------------------------------------------------------------------------------------------
  if (NoiseReEstimate)
  { if (!CovReEstimate)
  { cat("Noise variance cannot be re-estimated without re-estimating the covariance parameters \n")
    cat("covReEstimate switched to TRUE \n")
    CovReEstimate = TRUE
  }
    
    obs.n.rep <- pmax( 1, round(noise.var / as.numeric(model@noise.var)) ) 
    
    if (!is.null(estim.model))
      #-- Case 1.1: estim.model is provided (noise.var is set to the nugget of estim.model) -----------------
      #-- All the data is assumed to be in the appropriate format -------------------------------------------
    { noise.var <- estim.model@covariance@nugget
    } else
      
      #-- Case 1.2: estim.model needs to be built (noise.var is estimated) ----------------------------------
    { if ( max(obs.n.rep)==1 )
      { estim.model <- km(formula=model@trend.formula, covtype=model@covariance@name,
                         design=model@X, response=model@y,
                         lower=model@lower, upper=model@upper,
                          coef.cov=covparam2vect(model@covariance), coef.var=model@covariance@sd2,
                          optim.method=model@optim.method, 
                         nugget=noise.var, control=model@control)
        if (type=="SK"){ estim.model@trend.coef=model@trend.coef }
        estim.model@covariance@nugget.estim =TRUE
      } else 
      { cat("Unless estim.model is provided, noisy.EGO cannot estimate the noise variance when the initial 
            model has heterogeneous noise variance. NoiseReEstimate switched to FALSE. \n")
        NoiseReEstimate <- FALSE
      }
    }
  }

  ##########################################################################################################
  ### MAIN LOOP STARTS #####################################################################################
  for (i.time.steps in 1:n.ite)
  {
    add.obs <- TRUE

    # Update number of repetitions in the unknown noise variance case
    if (NoiseReEstimate) {  obs.n.rep <- pmax( 1, round(noise.var / as.numeric(model@noise.var)) )  }

    #-------------------------------------------------------------------------------------------------------
    #-- Find the new observation based on the REINTERPOLATION criterion ------------------------------------
    #-------------------------------------------------------------------------------------------------------
    if (optim.crit == "reinterpolation")
    {
      # Build interpolating model
      pred <- predict(object=model, newdata=model@X, type=type, checkNames = FALSE)
      mk <- pred$mean
      optim.param$ymin <- min(mk)
      if(exists("optim.model")) rm(optim.model)
      
      try(optim.model <- km(formula=model@trend.formula, design=model@X, response=mk,
    	            covtype=model@covariance@name, coef.trend=model@trend.coef, coef.cov=covparam2vect(model@covariance),
    	            coef.var=model@covariance@sd2,control=model@control))

      # If km crashes: add nugget to the interpolating model
      if(!exists("optim.model"))
      { print("Error occured during model update - small nugget added to the reinterpolating model")
        optim.model <- km(formula=model@trend.formula, design=model@X, response=mk,
    	            covtype=model@covariance@name, coef.trend=model@trend.coef, coef.cov=covparam2vect(model@covariance),
    	            coef.var=model@covariance@sd2, nugget=1e-8,control=model@control)
      }

      # Choose new observation based on EI
      oEGO <- max_EI(model=optim.model, plugin=min(optim.model@y), upper=upper, type=type, lower=lower, parinit=parinit, control=control)
      x.new <- oEGO$par
      i.best <- 0
    }
    #-------------------------------------------------------------------------------------------------------
    #-- Find the new observation based on the EI WITH PLUGIN criterion -------------------------------------
    #-------------------------------------------------------------------------------------------------------
    else if (optim.crit =="EI.plugin")
    {
      # Set plugin value (depending on "optim.param$plugin.type": "quantile", "ytilde" or "other")
      if (is.null(optim.param$plugin.type))
      { cat("Plugin type not provided: default value ytilde is used \n")
        optim.param$plugin.type <- "ytilde"}
      
      if (optim.param$plugin.type=="quantile")
      {  pred <- predict(object=model, newdata=model@X, type=type, checkNames = FALSE)
         mk <- pred$mean
         sk <- pred$sd
         if (is.null(optim.param$quantile))
         { cat("Quantile level not provided: default value 0.5 is used \n")
           optim.param$quantile <- 0.5
         }
         qk <- mk + qnorm(optim.param$quantile)*sk
         plugin <- min(qk)
      } else if (optim.param$plugin.type=="ytilde")
      {  plugin <- min(model@y)
      } else if (optim.param$plugin.type=="fixed")
      {  plugin <- optim.param$plugin
      } else 
      { cat("Unknown plugin type: default value ytilde is used \n")
        optim.param$plugin.type="ytilde"
        plugin <- min(model@y)
      }

      # Choose new observation based on EI.plugin
      oEGO <- max_EI(model=model, plugin=plugin, type=type, lower=lower, upper=upper, parinit=parinit, control=control)
      x.new <- oEGO$par
      EI <- oEGO$val

      # Compare with values at DoE
      EI.doe <- rep(0, 1, model@n)
      for (i in 1:model@n)
      {  EI.doe[i] <- EI(x=model@X[i,], model=model, type=type, plugin=plugin)}
      i.best <- which.max(EI.doe)

      # Change new observation if DoE is better
      if (EI.doe[i.best] > EI)
      {  x.new <- t(as.numeric(model@X[i.best,]))
         add.obs <- FALSE        }
    }
    #-------------------------------------------------------------------------------------------------------
    #-- Find the new observation based on the EXPECTED QUANTILE IMPROVEMENT criterion ----------------------
    #-------------------------------------------------------------------------------------------------------
    else if (optim.crit =="EQI")
    {
      # Set future noise variance, compute current minimal quantile
      if (is.null(optim.param))
      { beta <- 0.9
      } else { beta <- optim.param$quantile }

      new.noise.var <- noise.var/(n.ite+1-i.time.steps)
      pred <- predict(object=model, newdata=model@X, type=type, checkNames = FALSE)
      qk <- pred$mean + qnorm(beta)*pred$sd
      q.min <- min(qk)

      # Choose new observation based on EQI
      oEGO <- max_EQI(model=model, new.noise.var=new.noise.var, type=type,
                              beta=beta, q.min=q.min, lower=lower, upper=upper, parinit=parinit, control=control)
      x.new <- oEGO$par
      EQI.global.search <- oEGO$val

      # Compare with values at DoE
      EQI.doe <- rep(0, 1, model@n)
      for (i in 1:model@n)
      {  EQI.doe[i] <- EQI(x=model@X[i,], model=model, new.noise.var=new.noise.var, q.min=q.min, beta=beta, type=type)}
      i.best <- which.max(EQI.doe)

      # Change new observation if DoE is better
      if (EQI.doe[i.best] > EQI.global.search)
      {  x.new <- t(as.numeric(model@X[i.best,]))
         add.obs <- FALSE        }
    }
    #-------------------------------------------------------------------------------------------------------
    #-- Find the new observation based on the MINIMAL QUANTILE criterion -----------------------------------
    #-------------------------------------------------------------------------------------------------------
    else if (optim.crit =="min.quantile")
    {
      if (is.null(optim.param))
      { beta <- 0.1 
      } else { beta <- optim.param$quantile }

      # Choose new observation based on minimum quantile
      oEGO <- min_quantile(model=model, beta=beta, type=type, lower=lower, upper=upper, parinit=parinit, control=control)
      x.new <- oEGO$par
      q <- oEGO$val

      # Compare with values at DoE
      pred <- predict(object=model, newdata=model@X, type=type, checkNames = FALSE)
      mk <- pred$mean
      sk <- pred$sd
      q.doe <- mk + qnorm(beta)*sk
      i.best <- which.min(q.doe)

      # Change new observation if DoE is better
      if (q.doe[i.best] < q)
      {  x.new <- t(as.numeric(model@X[i.best,]))
         add.obs <- FALSE        }
    }
    #-------------------------------------------------------------------------------------------------------
    #-- Find the new observation based on the AUGMENTED EXPECTED IMPROVEMENT criterion ---------------------
    #-------------------------------------------------------------------------------------------------------
    else if (optim.crit =="AEI")
    {
      if (is.null(optim.param))
      { beta <- 0.75
      } else { beta <- optim.param$quantile }

      # Find current best and y.min
      pred <- predict(object=model, newdata=model@X, type=type, checkNames = FALSE)
      mk <- pred$mean
      sk <- pred$sd
      qk <- mk + qnorm(beta)*sk
      y.min <- mk[which.min(qk)]

      # Choose new observation based on AEI
      oEGO <- max_AEI(model=model, y.min=y.min, new.noise.var=noise.var, lower=lower, upper=upper, parinit=parinit, control=control, type=type)
      x.new <- oEGO$par
      AEI <- oEGO$val

      # Compare with values at DoE
      AEI.doe <- rep(0, 1, model@n)
      for (i in 1:model@n)
      {  AEI.doe[i] <- AEI(x=model@X[i,], model=model, y.min=y.min, type=type, new.noise.var=noise.var)}
      i.best <- which.max(AEI.doe)

      # Change new observation if DoE is better
      if (AEI.doe[i.best] > AEI)
      {  x.new <- t(as.numeric(model@X[i.best,]))
         add.obs <- FALSE        }
    }

    #-------------------------------------------------------------------------------------------------------
    #-- Find the new observation based on the APPROXIMATE KNOWLEDGE GRADIENT criterion ---------------------
    #-------------------------------------------------------------------------------------------------------
    else if (optim.crit =="AKG")
    {
      # Choose new observation based on AKG
      oEGO <- max_AKG(model=model, new.noise.var=noise.var, type=type, lower=lower, upper=upper, parinit=parinit, control=control)
      x.new <- oEGO$par
      AKG <- oEGO$val

      # Compare with values at DoE
      AKG.doe <- rep(0, 1, model@n)
      for (i in 1:model@n)
      {  AKG.doe[i] <- AKG(x=model@X[i,], model=model, type=type, new.noise.var=noise.var)}
      i.best <- which.max(AKG.doe)

      # Change new observation if DoE is better
      if (AKG.doe[i.best] > AKG)
      {  x.new <- t(as.numeric(model@X[i.best,]))
         add.obs <- FALSE        }
    }
    #-------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------
    #-- Generate observation, update model -----------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------
    # Generate observation
    y.new <- funnoise(x.new)

    # Update model
    if (add.obs)
    { cat(c("Creating obs:", as.numeric(x.new),"\n"))
    } else
    { cat(c("Repeating obs:", as.numeric(x.new),"\n"))}
 
    upmod <- update_km_noisyEGO(model=model, x.new=x.new, y.new=y.new, noise.var=noise.var, type=type,
                                add.obs=add.obs, index.in.DOE=i.best, CovReEstimate=CovReEstimate, NoiseReEstimate=NoiseReEstimate, 
                                estim.model=estim.model, nugget.LB=nugget.LB)
    model <- upmod$model
    if (NoiseReEstimate)
    { estim.model <- upmod$estim.model
      noise.var   <- upmod$noise.var
    }
    
    #-------------------------------------------------------------------------------------------------------
    #-- Save X path, observations and thetas ---------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------
    mu    <- model@trend.coef
    sd2   <- model@covariance@sd2
    range <- model@covariance@range.val
    theta <- c(mu, sd2, range)

    if (i.time.steps==1)
    {  all.X <- t(x.new)
       all.y <- t(y.new)
       all.thetas <- theta
       all.noise.var <- noise.var
    } else
    {  all.X <- cbind(all.X, t(x.new))
       all.y <- c(all.y, t(y.new))
       all.thetas <- cbind(all.thetas, theta)
       all.noise.var <- c(all.noise.var, t(noise.var))
    }
  }
  ### MAIN LOOP ENDS #######################################################################################
  ##########################################################################################################
  
  ############################################################################################################
  ######## COMPUTE BEST POINT ################################################################################
  ############################################################################################################
  # All other config
  
  if (!exists("optim.param$quantile")) { beta <- 0.5 
  } else { beta <- optim.param$quantile }
  
  pred <- predict(object=model, newdata=model@X, type=type, checkNames = FALSE)
  mk <- pred$mean
  sk <- pred$sd
  qk <- mk + qnorm(beta)*sk
  
  i.best <- which.min(qk)
  x.best <- model@X[i.best,]
  y.best <- model@y[i.best]
  
  if (optim.crit=="EI.plugin") 
  { if (optim.param$plugin.type=="ytilde")
    { # EI.plugin with ytilde treated differently
      i.best <- which.min(model@y)
      x.best <- model@X[i.best,]
      y.best <- model@y[i.best]
    }
  } 
  
  ############################################################################################################
  ######## RETURN OPTIMIZATION RESULTS #######################################################################
  ############################################################################################################
  optim.result$model <- model
  
  optim.result$best.x <- x.best
  optim.result$best.y <- y.best
  optim.result$best.index <- i.best
  
  optim.result$history.hyperparam <- all.thetas
  optim.result$history.x <- all.X
  optim.result$history.y <- all.y
  
  if (NoiseReEstimate)
  { optim.result$estim.model <- estim.model
    optim.result$history.noise.var <- all.noise.var}

  return(optim.result)
}
