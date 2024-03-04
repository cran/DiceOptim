##' Augmented Expected Improvement
##' 
##' Evaluation of the Augmented Expected Improvement (AEI) criterion, which is
##' a modification of the classical EI criterion for noisy functions.  The AEI
##' consists of the regular EI multiplied by a penalization function that
##' accounts for the disminishing payoff of observation replicates. The current
##' minimum y.min is chosen as the kriging predictor of the observation with
##' smallest kriging quantile.
##' 
##' 
##' @param x the input vector at which one wants to evaluate the criterion
##' @param model a Kriging model of "km" class
##' @param new.noise.var the (scalar) noise variance of the future observation.
##' @param y.min The kriging predictor at the current best point (point with
##' smallest kriging quantile).  If not provided, this quantity is evaluated.
##' @param type Kriging type: "SK" or "UK"
##' @param envir environment for saving intermediate calculations and reusing
##' them within AEI.grad
##' @return Augmented Expected Improvement
##' @author Victor Picheny  
##' 
##' David Ginsbourger 
##' @references D. Huang, T.T. Allen, W.I. Notz, and N. Zeng (2006), Global Optimization of
##' Stochastic Black-Box Systems via Sequential Kriging Meta-Models,
##' \emph{Journal of Global Optimization}, 34, 441-466.
##' @keywords models
##' @examples
##' 
##' 
##' ##########################################################################
##' ###    AEI SURFACE ASSOCIATED WITH AN ORDINARY KRIGING MODEL        ####
##' ### OF THE BRANIN FUNCTION KNOWN AT A 12-POINT LATIN HYPERCUBE DESIGN ####
##' ##########################################################################
##' 
##' set.seed(421)
##' 
##' # Set test problem parameters
##' doe.size <- 12
##' dim <- 2
##' test.function <- get("branin2")
##' lower <- rep(0,1,dim)
##' upper <- rep(1,1,dim)
##' noise.var <- 0.2
##' 
##' # Generate DOE and response
##' doe <- as.data.frame(matrix(runif(doe.size*dim),doe.size))
##' y.tilde <- rep(0, 1, doe.size)
##' for (i in 1:doe.size)  {
##' y.tilde[i] <- test.function(doe[i,]) + sqrt(noise.var)*rnorm(n=1)
##' }
##' y.tilde <- as.numeric(y.tilde)
##' 
##' # Create kriging model
##' model <- km(y~1, design=doe, response=data.frame(y=y.tilde),
##'             covtype="gauss", noise.var=rep(noise.var,1,doe.size), 
##' 	    lower=rep(.1,dim), upper=rep(1,dim), control=list(trace=FALSE))
##' 
##' # Compute actual function and criterion on a grid
##' n.grid <- 12 # Change to 21 for a nicer picture
##' x.grid <- y.grid <- seq(0,1,length=n.grid)
##' design.grid <- expand.grid(x.grid, y.grid)
##' nt <- nrow(design.grid)
##' crit.grid <- rep(0,1,nt)
##' func.grid <- rep(0,1,nt)
##' 
##' crit.grid <- apply(design.grid, 1, AEI, model=model, new.noise.var=noise.var)
##' func.grid <- apply(design.grid, 1, test.function)
##' 
##' # Compute kriging mean and variance on a grid
##' names(design.grid) <- c("V1","V2")
##' pred <- predict.km(model, newdata=design.grid, type="UK")
##' mk.grid <- pred$m
##' sk.grid <- pred$sd
##' 
##' # Plot actual function
##' z.grid <- matrix(func.grid, n.grid, n.grid)
##' filled.contour(x.grid,y.grid, z.grid, nlevels=50, color = rainbow,
##' plot.axes = {title("Actual function");
##' points(model@@X[,1],model@@X[,2],pch=17,col="blue"); 
##' axis(1); axis(2)})
##' 
##' # Plot Kriging mean
##' z.grid <- matrix(mk.grid, n.grid, n.grid)
##' filled.contour(x.grid,y.grid, z.grid, nlevels=50, color = rainbow,
##' plot.axes = {title("Kriging mean");
##' points(model@@X[,1],model@@X[,2],pch=17,col="blue"); 
##' axis(1); axis(2)})
##' 
##' # Plot Kriging variance
##' z.grid <- matrix(sk.grid^2, n.grid, n.grid)
##' filled.contour(x.grid,y.grid, z.grid, nlevels=50, color = rainbow,
##' plot.axes = {title("Kriging variance");
##' points(model@@X[,1],model@@X[,2],pch=17,col="blue"); 
##' axis(1); axis(2)})
##' 
##' # Plot AEI criterion
##' z.grid <- matrix(crit.grid, n.grid, n.grid)
##' filled.contour(x.grid,y.grid, z.grid, nlevels=50, color = rainbow,
##' plot.axes = {title("AEI");
##' points(model@@X[,1],model@@X[,2],pch=17,col="blue"); 
##' axis(1); axis(2)})
##' 
##' @export AEI
AEI <- function(x, model, new.noise.var=0, y.min=NULL, type = "UK", envir=NULL)
{
  ########## Convert x in proper format(s) ###
  d <- length(x)
  if (d != model@d){ stop("x does not have the right size") }
  newdata.num <- as.numeric(x)
  newdata <- data.frame(t(newdata.num))
  colnames(newdata) = colnames(model@X)

  # Compute y.min if missing
  if (is.null(y.min))
  {
    pred <- predict(object=model, newdata=model@X, type=type, checkNames = FALSE)
    mk <- pred$mean
    sk <- pred$sd
    qk <- mk + qnorm(0.75)*sk
    y.min <- mk[which.min(qk)]
  }
  
  # Prediction en newdata en partant de X
  pred <- predict(object=model, newdata=newdata, type=type, checkNames = FALSE) 
  mk <- pred$mean
  sk <- pred$sd  
  
  xcr <- (y.min - mk)/sk 
  xcr.prob <- pnorm(xcr)
  xcr.dens <- dnorm(xcr)
  
  if (!is.null(envir)) 
  {	assign("xcr", xcr, envir=envir)
    assign("xcr.prob", xcr.prob, envir=envir)
    assign("xcr.dens", xcr.dens, envir=envir)
    assign("c", pred$c, envir=envir)
    assign("Tinv.c", pred$Tinv.c, envir=envir)
    assign("mk", mk, envir=envir)
    assign("sk", sk, envir=envir)
  }
  if (sk < sqrt(model@covariance@sd2)/1e6) { aei.val <- 0 }  else 
  { aei.val <- ((y.min - mk) * xcr.prob + sk * xcr.dens) * (1- sqrt(new.noise.var)/sqrt(new.noise.var + sk^2))}

  return(aei.val)
}
