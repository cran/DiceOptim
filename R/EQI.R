##' Expected Quantile Improvement
##' 
##' Evaluation of the Expected Quantile Improvement (EQI) criterion.
##' 
##' 
##' @param x the input vector at which one wants to evaluate the criterion
##' @param model a Kriging model of "km" class
##' @param new.noise.var (scalar) noise variance of the future observation.
##' Default value is 0 (noise-free observation).
##' @param beta Quantile level (default value is 0.9)
##' @param q.min Best kriging quantile. If not provided, this quantity is
##' evaluated.
##' @param type Kriging type: "SK" or "UK"
##' @param envir environment for saving intermediate calculations and reusing
##' them within EQI.grad
##' @return Expected Quantile Improvement
##' @author Victor Picheny 
##' 
##' David Ginsbourger 
##' 
##' @references 
##' Picheny, V., Ginsbourger, D., Richet, Y., Caplin, G. (2013). 
##' Quantile-based optimization of noisy computer experiments with tunable precision. 
##' \emph{Technometrics}, 55(1), 2-13.
##' @examples
##' 
##' 
##' ##########################################################################
##' ###    EQI SURFACE ASSOCIATED WITH AN ORDINARY KRIGING MODEL        ####
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
##' 
##' crit.grid <- apply(design.grid, 1, EQI, model=model, new.noise.var=noise.var, beta=.9)
##' func.grid <- apply(design.grid, 1, test.function)
##' 
##' # Compute kriging mean and variance on a grid
##' names(design.grid) <- c("V1","V2")
##' pred <- predict(model, newdata=design.grid, type="UK", checkNames = FALSE)
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
##' # Plot EQI criterion
##' z.grid <- matrix(crit.grid, n.grid, n.grid)
##' filled.contour(x.grid,y.grid, z.grid, nlevels=50, color = rainbow,
##' plot.axes = {title("EQI");
##' points(model@@X[,1],model@@X[,2],pch=17,col="blue"); 
##' axis(1); axis(2)})
##' 
##' @export EQI
EQI <- function(x, model, new.noise.var=0, beta=0.9, q.min=NULL, type = "UK", envir=NULL)
{
########## Convert x in proper format(s) ###
d <- length(x)
if (d != model@d){ stop("x does not have the right size") }
newdata.num <- as.numeric(x)
newdata <- data.frame(t(newdata.num))
colnames(newdata) = colnames(model@X)

######### Compute q.min if missing #########
if (is.null(q.min))
{ pred <- predict.km(model, newdata=model@X, type=type, checkNames = FALSE)
  q.min <- min(pred$mean + qnorm(beta)*pred$sd)  }

######### Compute prediction at x #########
predx <- predict.km(model, newdata=newdata, type=type, checkNames = FALSE)
mk.old <- predx$mean
sk.old <- predx$sd

######### Intermediate values ##########
v <- predx$Tinv.c
c <- predx$c

######### mq and sq ##########
if (sk.old < sqrt(model@covariance@sd2)/1e6)
{ mq <- mk.old
  sq <- 0
} else
{ mq <- mk.old + qnorm(beta) * sqrt((new.noise.var * sk.old^2)/(new.noise.var + sk.old^2))
  sq <- sk.old^2/sqrt(new.noise.var + sk.old^2)
}

######### EQI ##########
xcr <- (q.min - mq)/sq
xcr.prob <- pnorm(xcr)
xcr.dens <- dnorm(xcr)

if (!is.null(envir)) {
assign("c", predx$c, envir=envir)
assign("v", v, envir=envir)
assign("xcr.prob", xcr.prob, envir=envir)
assign("xcr.dens", xcr.dens, envir=envir)
assign("mq", mq, envir=envir)
assign("sq", sq, envir=envir)
assign("sk.old", sk.old, envir=envir)
}

if (sk.old < sqrt(model@covariance@sd2)/1e6){ eqi.val <- 0} else 
{ eqi.val <- sq*(xcr*xcr.prob + xcr.dens) }

return(eqi.val)
}
