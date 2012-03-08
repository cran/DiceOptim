#source("max_EI.plugin.R")
max_EI.plugin <-function(model, plugin, type = "UK", lower, upper, parinit=NULL, control=NULL) {

	EI.envir <- new.env()	
	environment(EI.plugin) <- environment(EI.plugin.grad) <- EI.envir 

	trend.names <- colnames(model@F)
	if (length(trend.names)==1) {
		if (trend.names=="(Intercept)") gr = EI.plugin.grad
	}
	else {
		gr = NULL
	}     
          
	d <- ncol(model@X)
    
	if (is.null(control$pop.size))  control$pop.size <- floor(10 + 5 * log(d))# floor(4 + 3 * log(d)) # 200
    	if (is.null(control$max.generations))  control$max.generations <- 10 #50
    	if (is.null(control$wait.generations))  control$wait.generations <- 2 #2
    	if (is.null(control$BFGSburnin)) control$BFGSburnin <- 2 #0
	if (is.null(parinit))  parinit <- lower + runif(d) * (upper - lower)
     
	domaine <- cbind(lower, upper)

	o <- genoud(EI.plugin, nvars=d, max=TRUE, 
		pop.size=control$pop.size,
		max.generations=control$max.generations, 
		wait.generations=control$wait.generations,
           	hard.generation.limit=TRUE, starting.values=parinit, MemoryMatrix=TRUE, 
           	Domains=domaine, default.domains=10, solution.tolerance=0.00001, gr=gr, 
		boundary.enforcement=2, lexical=FALSE, gradient.check=TRUE, BFGS=TRUE,
          	data.type.int=FALSE, hessian=FALSE, unif.seed=floor(runif(1,max=10000)), int.seed=floor(runif(1,max=10000)),
          	print.level=1, share.type=0, instance.number=0,
          	output.path="stdout", output.append=FALSE, project.path=NULL,
         	P1=50, P2=50, P3=50, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0,
           	P9mix=NULL, BFGSburnin=control$BFGSburnin,BFGSfn=NULL, BFGShelp=NULL,
         	cluster=FALSE, balance=FALSE, debug=TRUE, 
         	model=model, plugin=plugin, type=type, envir=EI.envir 
		)
                            
    o$par <- t(as.matrix(o$par))
	colnames(o$par) <- colnames(model@X)
	o$value <- as.matrix(o$value)
	colnames(o$value) <- "EI_plugin"  
	return(list(par=o$par, value=o$value)) 
}
