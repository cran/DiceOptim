`max_qEI.CL` <- function(model, npoints, L, lower, upper, parinit=NULL, control=NULL) {

	n1 <- nrow(model@X)
		
   for (s in 1:npoints) {
		oEGO <- max_EI(model=model, lower=lower, upper=upper, parinit=parinit, control=control)
  		model@X <- rbind(model@X, oEGO$par)
      	model@y <- rbind(model@y, L, deparse.level=0)   		
		model@F <- trendMatrix.update(model, Xnew=data.frame(oEGO$par))
		
		if (model@noise.flag) {
				# heterogenous case : use 0 nugget for new points
				model@noise.var = c(model@covariance@nugget, 0)
		}
		
		model <- computeAuxVariables(model)
	
	}

	return(list(par = model@X[(n1+1):(n1+npoints),, drop=FALSE], value = model@y[(n1+1):(n1+npoints),, drop=FALSE]))
}
