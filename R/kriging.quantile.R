#source("kriging.quantile.R")

kriging.quantile <- function(x, model, alpha=0.1, type = "UK", envir=NULL)
{

d <- length(x)
newdata <- matrix(x, 1, d)

# Prediction en newdata en partant de X
predx <- predict.km(model, newdata=newdata, type="UK")
kriging.mean <- predx$mean
kriging.sd <- predx$sd  

qk <- kriging.mean + qnorm(alpha)*kriging.sd

if (!is.null(envir)) {
	assign("kriging.mean", kriging.mean, envir=envir)
	assign("kriging.sd", kriging.sd, envir=envir)
	assign("c", predx$c, envir=envir)
 	assign("Tinv.c", predx$Tinv.c, envir=envir)
 }

return(res <- qk)
}
