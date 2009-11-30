EI.grad <- function(x, model, envir) {
	
   get("kriging.sd", envir=envir)

	if (kriging.sd < 1e-06) {
		d <- length(x)
	   	rep(0,d)
   	} else {
		toget <- matrix(c("xcr", "xcr.prob", "xcr.dens", "c", "Tinv.c"),1,5)
    	apply(toget, 2, get, envir=envir)
		
		T <- model@T
		X <- model@X
		z <- model@z
		covStruct <- model@covariance
		c <- envir$c # Rajout David 15 Avril 2009
		Tinv.c<-envir$Tinv.c # Rajout David 15 Avril 2009

		dc <- covVector.dx(as.numeric(x), X, covStruct, c)	
		W <- backsolve(t(T), dc, upper.tri=FALSE)
		kriging.mean.grad <- t(W)%*%z
	
		u <- model@M
		v <- Tinv.c
		aux <- t(W)%*%v + (t(W)%*%u)* as.numeric((1-t(v)%*%u)/(t(u)%*%u))
     	kriging.sd.grad <- - aux / envir$kriging.sd 
# Modif (rajout "envir$") David 15 Avril 2009

		return(- kriging.mean.grad * envir$xcr.prob + kriging.sd.grad * envir$xcr.dens) 
# Modif (rajout "envir$") David 15 Avril 2009
	}
}
