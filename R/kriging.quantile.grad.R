#source("kriging.quantile.R")

kriging.quantile.grad <- function(x, model, alpha=0.1, type="UK")
{


	d <- length(x)
	newdata <- matrix(x, 1, d)

	predx <- predict(object=model, newdata=newdata, type=type, checkNames = FALSE)
	kriging.mean <- predx$mean
	kriging.sd <- predx$sd

	T <- model@T
	X <- model@X
	z <- model@z
	covStruct <- model@covariance
	c <- predx$c 
	Tinv.c <- predx$Tinv.c
	#dc <- covVector.dx(as.numeric(x), X, covStruct, c)	
	#dc <- covVector.dx(x=as.numeric(x), X=X, object=covStruct, c=c)
	dc <- covVector.dx(x=as.numeric(newdata), X=model@X, object=model@covariance, c=c)
	W <- backsolve(t(T), dc, upper.tri=FALSE)
	kriging.mean.grad <- t(W)%*%z	

	u <- model@M
	v <- Tinv.c
	aux <- t(W)%*%v + (t(W)%*%u)* as.numeric((1-t(v)%*%u)/(t(u)%*%u))
     	kriging.sd.grad <- - aux / kriging.sd 

return(quantile.grad <- kriging.mean.grad + qnorm(alpha)*kriging.sd.grad)
}
