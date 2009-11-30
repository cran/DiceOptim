EGO.nsteps <-function(model, fun, nsteps, lower, upper, parinit=NULL, control=NULL, kmcontrol=NULL) {

	n <- nrow(model@X)
	
	if (is.null(kmcontrol$penalty)) kmcontrol$penalty <- model@penalty
	if (length(model@penalty==0)) kmcontrol$penalty <- NULL 

	if (is.null(kmcontrol$optim.method)) kmcontrol$optim.method <- model@optim.method 
	if (is.null(kmcontrol$parinit)) kmcontrol$parinit <- model@parinit
	if (is.null(kmcontrol$control)) kmcontrol$control <- model@control
		

	for (i in 1:nsteps) {
		oEGO<-max_EI(model, lower, upper, parinit, control)
		
		model@X<-rbind(model@X, oEGO$par)
		model@y<-rbind(model@y, fun(t(oEGO$par)))
		
		kmcontrol$parinit <- covparam2vect(model@covariance)
		
		if (model@param.estim) {
			model <- km(formula=model@trend.formula, design=model@X, response=model@y, 
		   		 covtype=model@covariance@name, lower=model@lower, upper=model@upper, 
                 nugget=NULL, penalty=kmcontrol$penalty, optim.method=kmcontrol$optim.method, 
		    	parinit=kmcontrol$parinit, control=kmcontrol$control, gr=model@gr)
		} else {
			coef.cov <- covparam2vect(model@covariance)
			model <- km(formula=model@trend.formula, design=model@X, response=model@y, 
		   		 covtype=model@covariance@name, coef.trend=model@trend.coef, coef.cov=coef.cov, coef.var=model@covariance@sd2, nugget=NULL)
		}
	}
	
	return(list(par=model@X[(n+1):(n+nsteps),, drop=FALSE], 
               value=model@y[(n+1):(n+nsteps),, drop=FALSE], npoints=1, nsteps=nsteps, lastmodel=model))

}
