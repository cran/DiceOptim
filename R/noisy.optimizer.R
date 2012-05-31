#source("noisy.optimizer.R")
noisy.optimizer <- function(optim.crit, optim.param=NULL, model, n.ite, noise.var=NULL, funnoise, lower, upper, parinit=NULL, control=NULL,CovReEstimate=TRUE,
                             noiseReEstimate=FALSE, nugget.LB=1e-5, obs.n.rep=NULL, estim.model=NULL)

############################################################################################################
# The methods "random.search" and "SPO" are treated separately since they do not use
# a sequential update of the noisy kriging model.
#
# The methods "reinterpolation", "EQI", "min.quantile", "EI.plugin" are in the same loop. A switch is made
# at every iteration to choose the corresponding infill criterion, then the model is updated similarly for
# all methods.
#
# All methods return a final noisy kriging model. In addition, the quantities "all.X", "all.y", "all.thetas"
# (size=n.ite) are saved for analysis.
#
# This version allows sequential re-estimation of the noise variance. An additional model is created specifically
# for parameter estimation (estim.model), with uniform nugget (repeated experiments are not aggregated).
# "model" and "estim.model" are equivalent, but "model" has less observations ( hence goes faster).
#
# Additional (optional) arguments:
#    - "noiseReEstimate": TRUE/FALSE
#    - "noise.var": initial guess for the unknown variance (used in optimization)
#    - "obs.n.rep": number of repetitions per observation point. Required if "model" has heterogeneous variances
#    - "estim.model": model with homogeneous nugget effect (no noise.var). Required for restarting algorithms
#
############################################################################################################
{
print("Starting noisy optimization with the following criterion and parameters:")
print(optim.crit)
print(optim.param)
optim.result <- list()

############################################################################################################
########################                                                            ########################
########################            Optimizer 1: RANDOM SEARCH                      ########################
########################                                                            ########################
############################################################################################################
if (optim.crit =="random.search")
{
  ##########################################################################################################
  ### Generate new points randomly  ########################################################################
  ##########################################################################################################
  x.new <- matrix(0, n.ite, model@d)
  y.new <- rep(0, 1, n.ite)

  for (i in 1:n.ite)
  {   x.new[i,] <- runif(model@d)
      y.new[i]  <- funnoise(x.new[i,])
  }

  ##########################################################################################################
  ### Update model #########################################################################################
  ##########################################################################################################

  #-- Case 1: unknown noise variance -----------------------------------------------------------------------
  if (noiseReEstimate==TRUE)
  { model@control$upper.alpha <- model@covariance@sd2 / (model@covariance@sd2 + nugget.LB)

    # Build model for estimating nugget and covariance
    estim.model <- km(formula=model@trend.formula, design=data.frame(rbind(model@X, as.matrix(x.new))), response=data.frame(rbind(model@y, as.matrix(y.new))),
              covtype=model@covariance@name, lower=model@lower, upper=model@upper,
              parinit=covparam2vect(model@covariance), control=model@control,
              nugget=noise.var, nugget.estim=TRUE)

    # Build model without repetitions
    noise.var <- estim.model@covariance@nugget
    model@covariance@sd2 <- estim.model@covariance@sd2
    model@covariance@range.val  <- covparam2vect(estim.model@covariance)
    model@trend.coef <- estim.model@trend.coef
    model@noise.var <- rep(noise.var, model@n)
    model <- computeAuxVariables(model)
    model <- update_km(model=model, NewX=as.matrix(x.new), NewY=as.matrix(y.new), CovReEstimate=FALSE, new.noise.var=rep(noise.var,n.ite))

  } else

  #-- Case 2: known noise variance -----------------------------------------------------------------------
  {  model <- update_km(model=model, NewX=x.new, NewY=as.matrix(y.new), CovReEstimate=CovReEstimate,new.noise.var=rep(noise.var,n.ite))}

  ##########################################################################################################
  ### Save X path, observations and thetas #################################################################
  ##########################################################################################################
  mu    <- model@trend.coef
  sd2   <- model@covariance@sd2
  range <- model@covariance@range.val

  theta <- c(mu, sd2, range)
  all.thetas <- matrix(0,length(theta),n.ite)
  for (i in 1:n.ite) {all.thetas[,i] <- theta}
  all.X <- t(x.new)
  all.y <- t(y.new)
  all.noise.var <- rep(noise.var,n.ite)

############################################################################################################
########################                                                            ########################
########################            Optimizer 2: SPO (v0.4)                         ########################
########################                                                            ########################
############################################################################################################
} else if (optim.crit == "SPO")
{
  #---------------------------------------------------------------------------------------------------------
  #-- STEP 1: Build optim model (noiseless model) ----------------------------------------------------------
  #---------------------------------------------------------------------------------------------------------
  optim.model <- try(km(formula=model@trend.formula, design=model@X, response=model@y, covtype=model@covariance@name, lower=model@lower, upper=model@upper))

  # Add nugget if km crashes
  if (typeof(optim.model)=="character")
  { print("Error occured during SPO optim model construction - small nugget added to the interpolating model")
    optim.model <- km(formula=model@trend.formula, design=model@X, response=model@y, covtype=model@covariance@name, lower=model@lower, upper=model@upper, nugget=1e-6)
  }

  #---------------------------------------------------------------------------------------------------------
  #-- STEP 2: Counters initializations ---------------------------------------------------------------------
  #---------------------------------------------------------------------------------------------------------
  i.time.steps <- 0
  i.best <- which.min(model@y)

  if (is.null(optim.param))
  { nrep <- rep(1,1,model@n)
    history.best <- i.best
  } else
  { nrep <- optim.param$nrep
    history.best <- optim.param$history.best
  }
  r <- max(nrep)

  ##########################################################################################################
  ### MAIN LOOP STARTS #####################################################################################
  while (i.time.steps < n.ite)
  { 
    #-------------------------------------------------------------------------------------------------------
    #-- Choose and generate new observations ---------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------

    #-- Choose new observation based on EI ----------------------------------------------------------------- 
    oEGO <- max_EI(model=optim.model, lower=lower, upper=upper, parinit=parinit, control=control)
    x.new <- oEGO$par

    #-- Generate new observations --------------------------------------------------------------------------
    n.new <- min((n.ite-i.time.steps), r)
    y.new <- rep(0,1,n.new)
    for (i in 1:n.new) {y.new[i] <- funnoise(x.new)}
    i.time.steps <- i.time.steps + n.new
    nrep <- cbind(nrep, n.new)

    #-- Update current best with new repetitions -----------------------------------------------------------
    n.best.add <- 0
    if (i.time.steps < n.ite && nrep[i.best] < r)
    {	n.best.add <- min(n.ite-i.time.steps, r-nrep[i.best])
        y.best.add <- rep(0,1,n.best.add)
        for (i in 1:n.best.add) {y.best.add[i] <- funnoise(optim.model@X[i.best,])}
        i.time.steps <- i.time.steps + n.best.add
    }

    #-------------------------------------------------------------------------------------------------------
    #-- Update optim model ---------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------

    #-- Account for repeated observations ------------------------------------------------------------------
    if (n.best.add > 0)
    {   optim.model@y[i.best] <- (n.best.add + nrep[i.best])^-1 * (mean(y.best.add)*n.best.add + optim.model@y[i.best]*nrep[i.best])
	 nrep[i.best] <- nrep[i.best] + n.best.add  }

    #-- Add new observations, reestimate hyperparameters if specified --------------------------------------
    newmodel <- try(update_km(model=optim.model, NewX=as.matrix(x.new), NewY=as.matrix(mean(y.new)), CovReEstimate=CovReEstimate))

    #-- If km crashes: try to build a model with old hyperparameters ---------------------------------------
    if (typeof(newmodel)=="character")
    { print("Error occured during hyperparameter reestimation")
      newmodel <-  try(update_km(model=optim.model, NewX=as.matrix(x.new), NewY=as.matrix(mean(y.new)), CovReEstimate=FALSE))

      # If km still crashes, add nugget
      if (typeof(newmodel)=="character")
      { print("Error occured during SPO optim model construction - small nugget added to the interpolating model")
        optim.model@X <- rbind(optim.model@X, as.matrix(x.new))
        optim.model@y <- as.matrix(c(optim.model@y, as.matrix(mean(y.new))))
        newmodel <- km(formula=optim.model@trend.formula, design=optim.model@X, response=optim.model@y, covtype=optim.model@covariance@name,
                       coef.cov=covparam2vect(optim.model@covariance), coef.var=optim.model@covariance@sd2, lower=optim.model@lower, upper=optim.model@upper, nugget=1e-6)
      }
    }
    optim.model <- newmodel

    #-------------------------------------------------------------------------------------------------------
    #-- Final updates and savings --------------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------

    #--- Update best ---------------------------------------------------------------------------------------
    i.best <- which.min(optim.model@y)

    if (min(abs(i.best-history.best)) > 1)
    {  history.best <- cbind(history.best, i.best)
    } else
    {  r <- r + 1}

    #--- Save X path, observations and thetas --------------------------------------------------------------
    # New data
    mu    <- optim.model@trend.coef
    sd2   <- optim.model@covariance@sd2
    range <- optim.model@covariance@range.val
    theta <- c(mu, sd2, range)
    X.new <- t(x.new)
    Y.new <- y.new
    Theta <- as.matrix(theta)

    # Duplicate data to account for repetitions
    if (n.new > 1)      { for (i in 2:n.new)       { X.new <- cbind(X.new, t(x.new))} }
    if (n.best.add > 0) { for (i in 1:n.best.add)  { X.new <- cbind(X.new, as.numeric(optim.model@X[i.best,]))}
                          Y.new <- c(Y.new, as.numeric(y.best.add))}
    if (n.new+n.best.add > 1) { for (i in 2:(n.new+n.best.add)) { Theta <- cbind(Theta, as.matrix(theta))} }

    # Update archives
    if (i.time.steps==n.new+n.best.add)
    {
      all.X <- X.new
      all.y <- Y.new
      all.thetas <- Theta
    } else
    {  all.X <- cbind(all.X, X.new)
       all.y <- c(all.y, Y.new)
       all.thetas <- cbind(all.thetas, Theta)
    }
  }
  ### MAIN LOOP ENDS #######################################################################################
  ##########################################################################################################

  #---------------------------------------------------------------------------------------------------------
  #--- Build final model -----------------------------------------------------------------------------------
  #---------------------------------------------------------------------------------------------------------
  
  #-- Case 1: unknown noise variance -----------------------------------------------------------------------
  if (noiseReEstimate==TRUE)
  {
    # Builds model for noise variance estimation
    if (is.null(estim.model))
    { X.estim.model <- rbind(as.matrix(model@X), t(all.X))
      y.estim.model <- c(model@y, all.y)
    } else
    { X.estim.model <- rbind(estim.model@X, t(all.X))
      y.estim.model <- c(estim.model@y, all.y)
    }
    
    # Hyperparameters estimation
    model@control$upper.alpha <- model@covariance@sd2 / (model@covariance@sd2 + nugget.LB)
    estim.model <- km(formula=model@trend.formula, 
                    design=data.frame(X.estim.model),
                    response=data.frame(y.estim.model),
                    covtype=model@covariance@name, nugget=noise.var, nugget.estim=TRUE, lower=model@lower, upper=model@upper, control=model@control)
    
    # Update final model
    noise.var <- estim.model@covariance@nugget
    model <- optim.model
    model@covariance@sd2 <- estim.model@covariance@sd2
    model@covariance@range.val  <- covparam2vect(estim.model@covariance)
    model@trend.coef <- estim.model@trend.coef
    model@noise.var <- as.numeric(noise.var*nrep^-1)
    model <- computeAuxVariables(model)
  } else

  #-- Case 2: known noise variance -------------------------------------------------------------------------
  {  model <- km(formula=optim.model@trend.formula, design=optim.model@X, response=optim.model@y,
    	            covtype=model@covariance@name, noise.var = (noise.var*nrep^-1), lower=model@lower, upper=model@upper)
  }
  all.noise.var <- rep(noise.var,n.ite)

  #---------------------------------------------------------------------------------------------------------
  #-- Add number of repetitions and history of best observation to the output (for potential restart) ------
  #---------------------------------------------------------------------------------------------------------
  optim.result$nrep <- nrep
  optim.result$history.best <- history.best

############################################################################################################
########################        All other methods:                                  ########################
########################        Reinterpolation, EQI, SKO, EI with plugin           ########################
########################                                                            ########################
############################################################################################################
} else
{
  #---------------------------------------------------------------------------------------------------------
  # Initialization for the unknown noise variance case
  #---------------------------------------------------------------------------------------------------------
  if (noiseReEstimate==TRUE)
  {  if (!is.null(estim.model))
     #-- Case 1.1: estim.model is provided (noise.var is set to the nugget of estim.model) -----------------
     #-- All the data is assumed to be in the appropriate format -------------------------------------------
     {   noise.var <- estim.model@covariance@nugget
     } else

     #-- Case 1.2: estim.model needs to be built (noise.var is estimated) ----------------------------------
     {  model@control$upper.alpha <- model@covariance@sd2 / (model@covariance@sd2 + nugget.LB)
        estim.model <- try(km(formula=model@trend.formula, covtype=model@covariance@name,
                          design=model@X, response=model@y,
                          lower=model@lower, upper=model@upper, parinit=covparam2vect(model@covariance),
                          nugget=noise.var, nugget.estim=TRUE, control=model@control))

        # If km crashes: build model with old hyperparameters
        if (typeof(estim.model)=="character")
        { print("Error in hyperparameter estimation of estim.model - old hyperparameter values used instead")
          estim.model <- try(km(formula=model@trend.formula, covtype=model@covariance@name,
                          design=model@X, response=model@y,
                          lower=model@lower, upper=model@upper,
                          coef.cov=covparam2vect(model@covariance), coef.var=model@covariance@sd2,
                          nugget=noise.var, control=model@control)) 
        }

        # If number of repetitions is not provided, default value is 1 repetition for each obs.
        if (is.null(obs.n.rep))    {obs.n.rep <- rep(1,model@n)}

        # Update model with new noise
        noise.var <- as.numeric(estim.model@covariance@nugget)
        model@covariance@sd2 <- as.numeric(estim.model@covariance@sd2)
        model@covariance@range.val  <- covparam2vect(estim.model@covariance)
        model@trend.coef <- estim.model@trend.coef
        model@noise.var <- as.numeric(noise.var*c(1/obs.n.rep))
        newmodel <- try(computeAuxVariables(model))

        # If km crashes: try again by increasing noise variance
        n.try <- 0
        while (typeof(newmodel)=="character" & n.try < 10)
        {  print("Error in model building - noise variance multiplied by 2")
           noise.var <- 2*noise.var
           model@noise.var <- as.numeric(noise.var*c(1/obs.n.rep))
           newmodel <- try(computeAuxVariables(model))
           n.try <- n.try + 1
           print(noise.var)
        }
        model <- newmodel
     }
  }

  ##########################################################################################################
  ### MAIN LOOP STARTS #####################################################################################
  for (i.time.steps in 1:n.ite)
  {
    add.obs <- TRUE

    # Update number of repetitions in the unknown noise variance case
    if (noiseReEstimate==TRUE) {  obs.n.rep <- pmax( 1, round(noise.var / as.numeric(model@noise.var)) )  }

    #-------------------------------------------------------------------------------------------------------
    #-- Find the new observation based on the REINTERPOLATION criterion ------------------------------------
    #-------------------------------------------------------------------------------------------------------
    if (optim.crit == "reinterpolation")
    {
      print("Building interpolating model")
      # Build interpolating model
      pred <- predict(object=model, newdata=model@X, type="UK", checkNames = FALSE)
      mk <- pred$mean
      optim.param$ymin <- min(mk)
      rm(optim.model)
      try(optim.model <- km(formula=model@trend.formula, design=model@X, response=mk,
    	            covtype=model@covariance@name, coef.trend=model@trend.coef, coef.cov=covparam2vect(model@covariance),
    	            coef.var=model@covariance@sd2))

      # If km crashes: add nugget to the interpolating model
      if(!exists("optim.model"))
      {
        print("Error occured during model update - small nugget added to the reinterpolating model")
        optim.model <- km(formula=model@trend.formula, design=model@X, response=mk,
    	            covtype=model@covariance@name, coef.trend=model@trend.coef, coef.cov=covparam2vect(model@covariance),
    	            coef.var=model@covariance@sd2, nugget=1e-8)
      }
      print(optim.model)

      # Choose new observation based on EI
      oEGO <- max_EI(model=optim.model, lower=lower, upper=upper, parinit=parinit, control=control)
      x.new <- oEGO$par
    }
    #-------------------------------------------------------------------------------------------------------
    #-- Find the new observation based on the EI WITH PLUGIN criterion -------------------------------------
    #-------------------------------------------------------------------------------------------------------
    else if (optim.crit =="EI.plugin")
    {
      # Set plugin value (depending on "optim.param$plugin.type": "quantile", "ytilde" or "other")
      if (optim.param$plugin.type=="quantile")
      {  pred <- predict(object=model, newdata=model@X, type="UK", checkNames = FALSE)
         mk <- pred$mean
         sk <- pred$sd
         qk <- mk + qnorm(optim.param$quantile)*sk
         plugin <- min(qk)
      } else if (optim.param$plugin.type=="ytilde")
      {  plugin <- min(model@y)
      } else if (optim.param$plugin.type=="other")
      {  plugin <- optim.param$plugin
      } else {print("Unknown plugin type")}

      # Choose new observation based on EI.plugin
      oEGO <- max_EI.plugin(model=model, plugin=plugin, lower=lower, upper=upper, parinit=parinit, control=control)
      x.new <- oEGO$par
      EI <- oEGO$val

      # Compare with values at DoE
      EI.doe <- rep(0, 1, model@n)
      for (i in 1:model@n)
      {  EI.doe[i] <- EI.plugin(x=model@X[i,],model=model, plugin=plugin)}
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
      { alpha <- 0.9
      } else { alpha <- optim.param$quantile }

      new.noise.var <- noise.var/(n.ite+1-i.time.steps)
      pred <- predict(object=model, newdata=model@X, type="UK", checkNames = FALSE)
      mk <- pred$mean
      sk <- pred$sd
      sk <- sk*sqrt(model@n/(model@n-1))
      qk <- mk + qnorm(alpha)*sk
      q.min <- min(qk)

      # Choose new observation based on EQI
      oEGO <- max_EQI(model=model, new.noise.var=new.noise.var, 
                              alpha=alpha, q.min=q.min, lower=lower, upper=upper, parinit=parinit, control=control)
      x.new <- oEGO$par
      EQI.global.search <- oEGO$val

      # Compare with values at DoE
      EQI.doe <- rep(0, 1, model@n)
      for (i in 1:model@n)
      {  EQI.doe[i] <- EQI(x=model@X[i,], model=model, new.noise.var=new.noise.var, q.min=q.min, alpha=alpha)}
      i.best <- which.max(EQI.doe)

      # Local gradient search near the best DoE point
      controlLocalSearch <- list()
      controlLocalSearch$pop.size <- 1
      controlLocalSearch$max.generations <- 1
      controlLocalSearch$wait.generations <- 0
      controlLocalSearch$BFGSburnin <- 0
      controlLocalSearch$BFGSmaxit <- 20
      oEGO <- max_EQI(model=model, new.noise.var=new.noise.var, 
                              alpha=alpha, q.min=q.min, lower=lower, upper=upper, parinit=model@X[i.best,], control=controlLocalSearch)
      x.new2 <- oEGO$par
      EQI.local.search <- oEGO$val

      # Choose best DoE if difference is not significant with respect to local or global search
      if ( (max(EQI.local.search, EQI.global.search) - EQI.doe[i.best])/EQI.doe[i.best] < 0.001 )
      {
         x.new <- t(as.numeric(model@X[i.best,]))
         add.obs <- FALSE
      } else
      {
         if (EQI.local.search > EQI.global.search ) { x.new <- x.new2 }
      }
    }
    #-------------------------------------------------------------------------------------------------------
    #-- Find the new observation based on the MINIMAL QUANTILE criterion -----------------------------------
    #-------------------------------------------------------------------------------------------------------
    else if (optim.crit =="min.quantile")
    {
      if (is.null(optim.param))
      { alpha <- 0.1 
      } else { alpha <- optim.param$quantile }

      # Choose new observation based on minimum quantile
      oEGO <- min_quantile(model=model, alpha=alpha, lower=lower, upper=upper, parinit=parinit, control=control)
      x.new <- oEGO$par
      q <- oEGO$val

      # Compare with values at DoE
      pred <- predict(object=model, newdata=model@X, type="UK", checkNames = FALSE)
      mk <- pred$mean
      sk <- pred$sd
      q.doe <- mk + qnorm(alpha)*sk
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
      { alpha <- 0.75
      } else { alpha <- optim.param$quantile }

      # Find current best and y.min
      pred <- predict(object=model, newdata=model@X, type="UK", checkNames = FALSE)
      mk <- pred$mean
      sk <- pred$sd
      qk <- mk + qnorm(alpha)*sk
      y.min <- mk[which.min(qk)]

      # Choose new observation based on AEI
      oEGO <- max_AEI(model=model, y.min=y.min, new.noise.var=noise.var, lower=lower, upper=upper, parinit=parinit, control=control)
      x.new <- oEGO$par
      AEI <- oEGO$val

      # Compare with values at DoE
      AEI.doe <- rep(0, 1, model@n)
      for (i in 1:model@n)
      {  AEI.doe[i] <- AEI(x=model@X[i,], model=model, y.min=y.min, new.noise.var=noise.var)}
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
      oEGO <- max_AKG(model=model, new.noise.var=noise.var, lower=lower, upper=upper, parinit=parinit, control=control)
      x.new <- oEGO$par
      AKG <- oEGO$val

      # Compare with values at DoE
      AKG.doe <- rep(0, 1, model@n)
      for (i in 1:model@n)
      {  AKG.doe[i] <- AKG(x=model@X[i,], model=model, new.noise.var=noise.var)}
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

    #-------------------------------------------------------------------------------------------------------
    #-- Case 1: unknown noise variance ---------------------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------
    if (noiseReEstimate==TRUE)
    {
      #-- Update estim.model and noise.var -----------------------------------------------------------------
      model@control$upper.alpha <- model@covariance@sd2 / (model@covariance@sd2 + nugget.LB)
      old.noise.var <- noise.var

      new.estim.model <- try(km(formula=estim.model@trend.formula, covtype=estim.model@covariance@name,
                       design=rbind(estim.model@X, as.matrix(x.new)),
                       response=rbind(estim.model@y, as.matrix(y.new)),
                       lower=estim.model@lower, upper=estim.model@upper, parinit=covparam2vect(estim.model@covariance),
                       nugget=noise.var, nugget.estim=TRUE, control=model@control))

      # if km crashes: update model with new observation and old hyperparameters
      if (typeof(new.estim.model)=="character")
      { print("Error in hyperparameter estimation of estim.model - old hyperparameter values used instead")
        new.estim.model <- try(km(formula=model@trend.formula, covtype=model@covariance@name,
                        design=rbind(estim.model@X, as.matrix(x.new)),
                        response=rbind(estim.model@y, as.matrix(y.new)),
                        lower=model@lower, upper=model@upper, coef.cov=covparam2vect(estim.model@covariance),
                        nugget=old.noise.var, control=model@control))
      }
      estim.model <- new.estim.model
      noise.var <- estim.model@covariance@nugget

      #-- Update model -------------------------------------------------------------------------------------
      # (1) If a new point is added to the DOE
      if (add.obs==TRUE)
      {
        model@covariance@sd2 <- estim.model@covariance@sd2
        model@covariance@range.val  <- covparam2vect(estim.model@covariance)
        model@trend.coef <- estim.model@trend.coef
        model@noise.var <- noise.var*c(1/obs.n.rep)

        newmodel <- try(computeAuxVariables(model))
        newmodel <- try(update_km(model=newmodel, NewX=as.matrix(x.new), NewY=as.matrix(y.new), CovReEstimate=FALSE,new.noise.var=noise.var))

        # If km crashes: update model with new observation and old hyperparameters
        if (typeof(newmodel)=="character")
        { print("Error occured during hyperparameter reestimation")
          newmodel <- update_km(model=model, NewX=as.matrix(x.new), NewY=as.matrix(y.new), CovReEstimate=FALSE,new.noise.var=old.noise.var)
        }

      } else
      # (2) If an already existing observation is repeated
      {
        model@covariance@sd2 <- estim.model@covariance@sd2
        model@covariance@range.val  <- covparam2vect(estim.model@covariance)
        model@trend.coef <- estim.model@trend.coef
        newmodel <- model
        newmodel@noise.var <- noise.var*c(1/obs.n.rep)
        newmodel@y[i.best] <- (1/noise.var + 1/model@noise.var[i.best])^-1 * (y.new/noise.var + model@y[i.best]/model@noise.var[i.best])
        newmodel@noise.var[i.best] <- noise.var*model@noise.var[i.best]/(noise.var+model@noise.var[i.best])
        newmodel <- try(computeAuxVariables(newmodel))
        
        # If km crashes: try to build it again with a noise variance multiplied by 2
        n.try <- 0
        while (typeof(newmodel)=="character" & n.try < 10)
        {  print("Error in model building - noise variance multiplied by 2")
           noise.var <- 2*noise.var
           newmodel <- model
           newmodel@noise.var <- noise.var*c(1/obs.n.rep)
           newmodel@y[i.best] <- (1/noise.var + 1/model@noise.var[i.best])^-1 * (y.new/noise.var + model@y[i.best]/model@noise.var[i.best])
           newmodel@noise.var[i.best] <- noise.var*model@noise.var[i.best]/(noise.var+model@noise.var[i.best])
           newmodel <- try(computeAuxVariables(newmodel))
           n.try <- n.try + 1
        }
      }
      model <- newmodel
      print("new model nugget")
      print(model@noise.var)

    } else

    #-------------------------------------------------------------------------------------------------------
    #-- Case 2: known noise variance (no estim.model) ------------------------------------------------------
    #-------------------------------------------------------------------------------------------------------
    {
      print("Update model")

      #-- (1) If a new point is added to the DOE -----------------------------------------------------------
      if (add.obs==TRUE)
      {
        print("Adding obs")
        newmodel <- try(update_km(model=model, NewX=as.matrix(x.new), NewY=as.matrix(y.new), CovReEstimate=CovReEstimate,new.noise.var=noise.var))

        print(newmodel)
        # If km crashes: update model with new observation and old hyperparameters
        if (typeof(newmodel)=="character")
        { print("Error occured during hyperparameter reestimation")
          newmodel <- update_km(model=model, NewX=as.matrix(x.new), NewY=as.matrix(y.new), CovReEstimate=FALSE,new.noise.var=noise.var)
        }
      } else

      #-- (2) If an already existing observation is repeated -----------------------------------------------
      {
        # Update observation value and noise variance
        model@y[i.best] <- (1/noise.var + 1/model@noise.var[i.best])^-1 * (y.new/noise.var + model@y[i.best]/model@noise.var[i.best])
        model@noise.var[i.best] <- noise.var*model@noise.var[i.best]/(noise.var+model@noise.var[i.best])

        # Reestimate hyperparameters if required
        if (CovReEstimate==TRUE)
        {
          newmodel <- try(km(formula=model@trend.formula, design=model@X, response=model@y, covtype=model@covariance@name, noise.var=model@noise.var,
                             parinit=covparam2vect(model@covariance),lower=model@lower, upper=model@upper))

          # If km crashes: try to build it again with a noise variance multiplied by 2
          if (typeof(newmodel)=="character")
          { print("Error occured during hyperparameter reestimation")
            newmodel <- km(formula=model@trend.formula, design=model@X, response=model@y, covtype=model@covariance@name, noise.var=model@noise.var,
                           coef.trend=model@trend.coef, coef.cov=covparam2vect(model@covariance), coef.var=model@covariance@sd2, lower=model@lower, upper=model@upper)
          }
        # If no hyperparameter reestimation
        } else
        {  newmodel <- computeAuxVariables(model)}
      }
      model <- newmodel
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
}

############################################################################################################
######## RETURN OPTIMIZATION RESULTS #######################################################################
############################################################################################################
optim.result$model <- model
optim.result$theta <- all.thetas
optim.result$X <- all.X
optim.result$y <- all.y
optim.result$noise.var <- all.noise.var
if (noiseReEstimate==TRUE)
{ optim.result$estim.model <- estim.model}

return(optim.result)
}
