#' Update of one or two Kriging models when adding new observation
#' 
#' Update of a noisy Kriging model when adding new observation, with or
#' without covariance parameter re-estimation.  When the noise level is
#' unkown, a twin model "estim.model" is also updated.
#' 
#' 
#' @param model a Kriging model of "km" class
#' @param x.new a matrix containing the new points of experiments
#' @param y.new a matrix containing the function values on the points NewX
#' @param noise.var scalar: noise variance
#' @param type kriging type: "SK" or "UK"
#' @param add.obs boolean: if TRUE, the new point does not exist already in
#' the design of experiment model@@X
#' @param index.in.DOE optional integer: if add.obs=TRUE, it specifies the
#' index of the observation in model@@X corresponding to x.new
#' @param CovReEstimate optional boolean specfiying if the covariance
#' parameters should be re-estimated (default value = TRUE)
#' @param NoiseReEstimate optional boolean specfiying if the noise variance
#' should be re-estimated (default value = TRUE)
#' @param estim.model optional input of "km" class. Required if
#' NoiseReEstimate=TRUE, in order to deal with repetitions.
#' @param nugget.LB optional scalar: is used to define a lower bound on the
#' noise variance.
#' @return A list containing: \item{model}{ The updated Kriging model }
#' \item{estim.model}{If NoiseReEstimate=TRUE, the updated estim.model}
#' \item{noise.var}{If NoiseReEstimate=TRUE, the re-estimated noise variance}
#' @author Victor Picheny 
#' 
#' @references V. Picheny and D. Ginsbourger (2013), Noisy kriging-based optimization
#' methods: A unified implementation within the DiceOptim package,
#' \emph{Computational Statistics & Data Analysis}
#' 
#' @export update_km_noisyEGO
update_km_noisyEGO <- function(model, x.new, y.new, noise.var=0, type="UK", add.obs=TRUE, index.in.DOE=NULL, 
                               CovReEstimate=TRUE, NoiseReEstimate=FALSE, estim.model=NULL, nugget.LB=1e-5)
{
  i.new <- index.in.DOE
  cov.reestim <- CovReEstimate
  trend.reestim <- FALSE
  if (type=="UK") {trend.reestim <- TRUE}
  
  #-------------------------------------------------------------------------------------------------------
  #-- Case 1: unknown noise variance ---------------------------------------------------------------------
  #-------------------------------------------------------------------------------------------------------
  if (NoiseReEstimate==TRUE)
  {
    obs.n.rep <- pmax( 1, round(noise.var / as.numeric(model@noise.var)) )
    
    #-- Update estim.model and noise.var -----------------------------------------------------------------
    model@control$upper.alpha <- model@covariance@sd2 / (model@covariance@sd2 + nugget.LB)
    old.noise.var <- noise.var 
    new.estim.model <- try(update(object=estim.model, newX=as.matrix(x.new),newy=as.matrix(y.new),newX.alreadyExist=FALSE,
                                     cov.reestim=cov.reestim, newnoise.var=NULL, kmcontrol=model@control, nugget.reestim =TRUE, trend.reestim=trend.reestim))
    
    # if km crashes: update model with new observation and old hyperparameters
    if (typeof(new.estim.model)=="character")
    { message("Error in hyperparameter estimation of estim.model - old hyperparameter values used instead",quote=FALSE)
      new.estim.model <- try(update(object=estim.model, newX=as.matrix(x.new),newy=as.matrix(y.new),newX.alreadyExist=FALSE,
                                       cov.reestim=FALSE, newnoise.var=NULL, kmcontrol=model@control, trend.reestim=trend.reestim))
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
      newmodel <- try(update(object=newmodel, newX=as.matrix(x.new), newy=as.matrix(y.new), cov.reestim=FALSE, newnoise.var=noise.var, kmcontrol=model@control, trend.reestim=trend.reestim))
      
      # If km crashes: update model with new observation and old hyperparameters
      if (typeof(newmodel)=="character")
      { message("Error occured during hyperparameter reestimation - old hyperparameter values used instead",quote=FALSE)
        newmodel <- update(object=model, newX=as.matrix(x.new), newy=as.matrix(y.new), cov.reestim=FALSE, newnoise.var=old.noise.var, kmcontrol=model@control, trend.reestim=trend.reestim)
      }
    } else
      # (2) If an already existing observation is repeated
    {
      model@covariance@sd2 <- estim.model@covariance@sd2
      model@covariance@range.val  <- covparam2vect(estim.model@covariance)
      model@trend.coef <- estim.model@trend.coef
      newmodel <- model
      newmodel@noise.var <- noise.var*c(1/obs.n.rep)
      newmodel@y[i.new] <- (1/noise.var + 1/model@noise.var[i.new])^-1 * (y.new/noise.var + model@y[i.new]/model@noise.var[i.new])
      newmodel@noise.var[i.new] <- noise.var*model@noise.var[i.new]/(noise.var+model@noise.var[i.new])
      newmodel <- try(computeAuxVariables(newmodel))
      
      # If km crashes: try to build it again with a noise variance multiplied by 2
      n.try <- 0
      while (typeof(newmodel)=="character" & n.try < 10)
      {  message("Error in model building - noise variance multiplied by 2",quote=FALSE)
         noise.var <- 2*noise.var
         newmodel <- model
         newmodel@noise.var <- noise.var*c(1/obs.n.rep)
         newmodel@y[i.new] <- (1/noise.var + 1/model@noise.var[i.new])^-1 * (y.new/noise.var + model@y[i.new]/model@noise.var[i.new])
         newmodel@noise.var[i.new] <- noise.var*model@noise.var[i.new]/(noise.var+model@noise.var[i.new])
         newmodel <- try(computeAuxVariables(newmodel))
         n.try <- n.try + 1
      }
    }
    model <- newmodel
    message(c("New model nugget: ", noise.var, "\n"))
    
  } else
  {  
  #-------------------------------------------------------------------------------------------------------
  #-- Case 2: known noise variance (no estim.model) ------------------------------------------------------
  #-------------------------------------------------------------------------------------------------------
  
  #-- (1) If a new point is added to the DOE -----------------------------------------------------------
  if (add.obs==TRUE)
  {  
    newmodel <- try(update(object=model, newX=as.matrix(x.new), newy=as.matrix(y.new), cov.reestim=cov.reestim, newnoise.var=noise.var,trend.reestim=trend.reestim))

    # If km crashes: update model with new observation and old hyperparameters
    if (typeof(newmodel)=="character")
    { message("Error occured during hyperparameter reestimation - old hyperparameter values used instead")
      newmodel <- update(object=model, newX=as.matrix(x.new), newy=as.matrix(y.new), cov.reestim=FALSE, newnoise.var=noise.var, trend.reestim=trend.reestim)
    }
  } else
    
  #-- (2) If an already existing observation is repeated -----------------------------------------------
  {
    # Update observation value and noise variance
    model@y[i.new] <- (1/noise.var + 1/model@noise.var[i.new])^-1 * (y.new/noise.var + model@y[i.new]/model@noise.var[i.new])
    model@noise.var[i.new] <- noise.var*model@noise.var[i.new]/(noise.var+model@noise.var[i.new])
    
    # Reestimate hyperparameters if required
    if (cov.reestim==TRUE)
    {
      if (type=="UK")
      { newmodel <- try(km(formula=model@trend.formula, design=model@X, response=model@y, covtype=model@covariance@name, noise.var=model@noise.var,
                         parinit=covparam2vect(model@covariance),lower=model@lower, upper=model@upper,control=model@control))
      } else
      { newmodel <- try(km(formula=model@trend.formula, design=model@X, response=model@y, covtype=model@covariance@name, noise.var=model@noise.var,
                           coef.trend=model@trend.coef,
                           parinit=covparam2vect(model@covariance),lower=model@lower, upper=model@upper,control=model@control))
      }
      
      # If km crashes: try to build it again with old hyperparameters
      if (typeof(newmodel)=="character")
      { message("Error occured during hyperparameter reestimation - old hyperparameter values used instead")
        newmodel <- km(formula=model@trend.formula, design=model@X, response=model@y, covtype=model@covariance@name, 
                       noise.var=model@noise.var, control=model@control,
                       coef.trend=model@trend.coef, coef.cov=covparam2vect(model@covariance), 
                       coef.var=model@covariance@sd2, lower=model@lower, upper=model@upper)
      }
    } else # If no hyperparameter reestimation
    {  newmodel <- computeAuxVariables(model)
    }
  }
  model <- newmodel
  }
  upmod <- list(model=model)
  if (NoiseReEstimate==TRUE)
  { upmod$estim.model <- estim.model
    upmod$noise.var   <- noise.var
  }
  return(upmod)
}
