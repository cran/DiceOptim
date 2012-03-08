max_EI <-function(model, lower, upper, parinit=NULL, control=NULL) {
        
	EI.envir <- new.env()	
	environment(EI) <- environment(EI.grad) <- EI.envir 
	#EI.optim(model, lower, upper, parinit, control)


	trend.names <- colnames(model@F)
	if (length(trend.names)==1) {
		if (trend.names=="(Intercept)") gr = EI.grad
	}
	else {
		gr = NULL
	}     

	d <- ncol(model@X)
    
	if (is.null(control$pop.size))  control$pop.size <- floor(4 + 3 * log(d))
    	if (is.null(control$max.generations))  control$max.generations <- 5
    	if (is.null(control$wait.generations))  control$wait.generations <- 2
    	if (is.null(control$BFGSburnin)) control$BFGSburnin <- 0
	if (is.null(parinit))  parinit <- lower + runif(d) * (upper - lower)
     
	domaine <- cbind(lower, upper)


	o <- genoud(EI, nvars=d, max=TRUE, 
		pop.size=control$pop.size,
		max.generations=control$max.generations, 
		wait.generations=control$wait.generations,
           	hard.generation.limit=TRUE, starting.values=parinit, MemoryMatrix=TRUE, 
           	Domains=domaine, default.domains=10, solution.tolerance=0.01,
          	gr=gr, boundary.enforcement=2, lexical=FALSE, gradient.check=FALSE, BFGS=TRUE,
          	data.type.int=FALSE, hessian=TRUE, unif.seed=floor(runif(1,max=10000)), int.seed=floor(runif(1,max=10000)),
          	print.level=0, share.type=0, instance.number=0,
          	output.path="stdout", output.append=FALSE, project.path=NULL,
         	P1=50, P2=50, P3=50, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0,
           	P9mix=NULL, 
		BFGSburnin=control$BFGSburnin,
          	BFGSfn=NULL, BFGShelp=NULL,
         	cluster=FALSE, balance=FALSE, debug=FALSE, model=model, envir=EI.envir)
                            
   o$par <- t(as.matrix(o$par))
	colnames(o$par) <- colnames(model@X)
	o$value <- as.matrix(o$value)
	colnames(o$value) <- "EI"
	return(list(par=o$par, value=o$value)) 

}
