EI <- function (x, model, type="UK", envir=NULL) {

	m <- min(model@y)
	d <- length(x)
   x <- matrix(x, 1, d)
   
   predx <- predict.km(object=model, newdata=x, type=type)
   kriging.mean <- predx$mean
   kriging.sd <- predx$sd
	
   xcr <- (m - kriging.mean)/kriging.sd
    
	if (kriging.sd < 1e-06) {
 	    	res <- 0
    	}
   else {
		xcr.prob <- pnorm(xcr)
      	xcr.dens <- dnorm(xcr)	        
	   	res <- (m - kriging.mean) * xcr.prob + kriging.sd * xcr.dens
		if (!is.null(envir)) {
			assign("xcr", xcr, envir=envir)
			assign("xcr.prob", xcr.prob, envir=envir)
	   		assign("xcr.dens", xcr.dens, envir=envir)
	   		}
	}
    
   if (!is.null(envir)) {
   		assign("kriging.sd", kriging.sd, envir=envir)
		assign("c", predx$c, envir=envir)
 		assign("Tinv.c", predx$Tinv.c, envir=envir)
 	}
   	
	return(res)
}
