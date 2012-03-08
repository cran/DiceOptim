#source("max_EQI.R")
max_EQI <-function(model, new.noise.var=0, alpha=0.9, q.min=NULL, type = "UK", lower, upper, parinit=NULL, control=NULL) {

	EQI.envir <- new.env()	
	environment(EQI) <- environment(EQI.grad_optim) <- EQI.envir 

	trend.names <- colnames(model@F)
	if (length(trend.names)==1) {
		if (trend.names=="(Intercept)") gr = EQI.grad_optim
	}
	else {
		gr = NULL
	}  

	d <- ncol(model@X)
	if(d<=6) N <- 3*2^d else N <- 32*d 
   	if (is.null(control$BFGSmaxit)) control$BFGSmaxit <- N
	if (is.null(control$solution.tolerance))  control$solution.tolerance <- 0
	if (is.null(control$pop.size))  control$pop.size <- N
    	if (is.null(control$max.generations))  control$max.generations <- 10
    	if (is.null(control$wait.generations))  control$wait.generations <- 2
    	if (is.null(control$BFGSburnin)) control$BFGSburnin <- 0
	if (is.null(parinit))  parinit <- improvedLHS(N, d)
#lower + runif(d) * (upper - lower)
     
	domaine <- cbind(lower, upper)

	o <- genoud(EQI, nvars=d, max=TRUE,
		pop.size=control$pop.size,
		max.generations=control$max.generations, 
		wait.generations=control$wait.generations,
           	hard.generation.limit=TRUE, starting.values=parinit, MemoryMatrix=TRUE, 
           	Domains=domaine, default.domains=10, solution.tolerance=control$solution.tolerance,
          	gr=gr, boundary.enforcement=2, lexical=FALSE, gradient.check=FALSE, BFGS=TRUE,
          	data.type.int=FALSE, hessian=TRUE, unif.seed=floor(runif(1,max=10000)), int.seed=floor(runif(1,max=10000)), print.level=1, 
		share.type=0, instance.number=0, output.path="stdout", output.append=FALSE, project.path=NULL,
         	P1=50, P2=50, P3=50, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0, P9mix=NULL, 
		BFGSburnin=control$BFGSburnin, BFGSfn=NULL, BFGShelp=NULL, control=list("maxit"=control$BFGSmaxit), 
		cluster=FALSE, balance=FALSE, debug=FALSE, model=model, 
                new.noise.var=new.noise.var, alpha=alpha, q.min=q.min, type=type, envir=EQI.envir
		)
                            
    o$par <- t(as.matrix(o$par))
	colnames(o$par) <- colnames(model@X)
	o$value <- as.matrix(o$value)
	colnames(o$value) <- "EQI"  
	return(list(par=o$par, value=o$value)) 
}
