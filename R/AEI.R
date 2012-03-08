#source("AEI.R")

AEI <- function(x, model, new.noise.var=0, y.min=NULL, type = "UK", envir=NULL)
{
#newdata <- t(newdata)
d <- length(x)
newdata <- matrix(x, 1, d)

# Compute y.min if missing
if (is.null(y.min))
{
  pred <- predict.km(model, newdata=model@X, type="UK")
  mk <- pred$mean
  sk <- pred$sd
  qk <- mk + qnorm(0.75)*sk
  y.min <- mk[which.min(qk)]
}

# Prediction en newdata en partant de X
pred <- predict.km(model, newdata, type="UK") 
mk <- pred$mean
sk <- pred$sd  

xcr <- (y.min - mk)/sk 
xcr.prob <- pnorm(xcr)
xcr.dens <- dnorm(xcr)

if (!is.null(envir)) {
   	#assign("pred", pred, envir=envir)
	assign("c", pred$c, envir=envir)
	assign("Tinv.c", pred$Tinv.c, envir=envir)
	assign("mk", mk, envir=envir)
 	assign("sk", sk, envir=envir)
   	assign("xcr", xcr, envir=envir)
	assign("xcr.prob", xcr.prob, envir=envir)
 	assign("xcr.dens", xcr.dens, envir=envir)
   }

res <- ((y.min - mk) * xcr.prob + sk * xcr.dens) * (1- sqrt(new.noise.var)/sqrt(new.noise.var + sk^2))

return(res)
}
