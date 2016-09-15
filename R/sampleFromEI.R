##' Samples \code{n} points from a distribution proportional to the expected
##' improvement (EI) computed from a \code{km} object.
##' 
##' @title Sampling points according to the expected improvement criterion
##' @param model an object of class \code{\link[DiceKriging]{km}},
##' @param minimization logical specifying if EI is used in minimiziation or in
##' maximization,
##' @param n number of points to be sampled,
##' @param initdistrib matrix of candidate points. 
##' @param lower vector of lower bounds,
##' @param upper vector of upper bounds,
##' @param T optional scalar : if provided, it replaces the current minimum (or
##' maximum) of observations.
##' @return A \code{n*d} matrix containing the sampled points. If \code{NULL}, \code{1000*d} points 
##' are obtained by latin hypercube sampling,
##' @author Sebastien Marmin 
##' 
##' Clement Chevalier 
##' 
##' David Ginsbourger 
##' @seealso \code{\link{EI}}, \code{\link{km}}, \code{\link{qEI}}
##' @references
##' 
##' D.R. Jones, M. Schonlau, and W.J. Welch (1998), Efficient global
##' optimization of expensive black-box functions, \emph{Journal of Global
##' Optimization}, 13, 455-492.
##' @keywords optimization
##' @examples
##' 
##' 
##' \donttest{
##' 
##' set.seed(004)
##' 
##' # a 9-points factorial design, and the corresponding responses
##' d <- 2
##' n <- 9
##' design.fact <- expand.grid(seq(0,1,length=3), seq(0,1,length=3)) 
##' names(design.fact)<-c("x1", "x2")
##' design.fact <- data.frame(design.fact) 
##' names(design.fact)<-c("x1", "x2")
##' response.branin <- apply(design.fact, 1, branin)
##' response.branin <- data.frame(response.branin) 
##' lower <- c(0,0)
##' upper <- c(1,1)
##' names(response.branin) <- "y" 
##' 
##' 
##' # model identification
##' fitted.model <- km(~1, design=design.fact, response=response.branin, 
##'                    covtype="gauss", control=list(pop.size=50,trace=FALSE), parinit=c(0.5, 0.5))
##' 
##' # sample a 30 point batch
##' batchSize <- 30
##' x <- sampleFromEI(model = fitted.model, n = batchSize, lower = lower, upper = upper)
##' 
##' # graphics 
##' # displays the EI criterion, the design points in black and the EI-sampled points in red.
##' nGrid <- 15
##' gridAxe1 <- seq(lower[1],upper[1],length=nGrid)
##' gridAxe2 <- seq(lower[2],upper[2],length=nGrid)
##' grid <- expand.grid(gridAxe1,gridAxe2)
##' aa <- apply(grid,1,EI,model=fitted.model)
##' myMat <- matrix(aa,nrow=nGrid)
##' image(x = gridAxe1, y = gridAxe2, z = myMat, 
##'       col = colorRampPalette(c("darkgray","white"))(5*10), 
##'       ylab = names(design.fact)[1], xlab=names(design.fact)[2], 
##'       main = "Sampling from the expected improvement criterion", 
##'       axes = TRUE, zlim = c(min(myMat), max(myMat)))
##' contour(x = gridAxe1, y = gridAxe2, z = myMat, 
##'         add = TRUE, nlevels = 10)
##' points(x[,1],x[,2],pch=19,col='red')
##' points(fitted.model@@X[,1],fitted.model@@X[,2],pch=19)
##' }
##' 
##' @export sampleFromEI
sampleFromEI <- function(model,minimization=TRUE,n= 1,initdistrib=NULL,lower=rep(0,model@d),upper=rep(1,model@d),T=NULL){
	
  #Generic function to sample from a density proportional to EI
  
  #### Arguments:
  # model : a km object
  # minimization : do we perform a global minimization or maximization ?
  # n : the total size of the sample distributed from a density proportional to EI
  # initdistrib : how these points are generated (default: "LHS" sequence)
  # d : dimension of the input set
  # lower, upper : arrays of size d with the lower (resp. upper) bounds in each dimension
  # T: the threshold used to compute the EI
  
  if(is.null(model)){
    print("Error in samplefromEI. A km object needs to be given")
    return(NULL)
  }
  
  d <- model@d
  
  if (length(lower) != length(upper) ){
    print("Error in samplefromEI : lower and upper must have the same length")
    return(NULL)
  }
  
  if(is.null(initdistrib)) {
	initial.points <- t(lower+t(randLHS(n=d*1000,k=d))*(upper-lower))
	} else {
	initial.points <- initdistrib
	}
	if(d==1) initial.points <- matrix(initial.points,ncol=1)
	#prediction on these initial candidate points
	predictions <- predict(object=model,newdata=initial.points,type="UK",checkNames=FALSE,cov.compute = FALSE,se.compute = TRUE,light.return = TRUE)
	  
  mn <- predictions$mean
  sn <- predictions$sd
	
  if(is.null(T)){
    if(minimization) T <- min(model@y)
    if(!minimization) T <- max(model@y)
  }
  
  uu <- (mn-T)/sn
  if(minimization) uu <- -1*uu
  ei <- sn*( uu*pnorm(uu) + dnorm(uu) )
  
  if ( minimization)  ei[is.nan(ei)] <- (T-mn)
  if (!minimization)  ei[is.nan(ei)] <- (mn-T)  
  ei <- ei*(ei>0)
  if(sum(ei!=0)<n) {
    maxEI <- max(ei)
    if (maxEI > 0) {
      ei <- ei/maxEI+1/(1000*d)
    } else {
      ei <- rep(1,(1000*d))
    }
  }
  my.indices <- sample(1:(1000*d), n, replace = FALSE, prob = ei)
  my.points <- initial.points[my.indices,]
		
	if(d==1) my.points <- matrix(my.points,ncol=1)
	if(n==1) my.points <- matrix(my.points,ncol=d)
    
	return(my.points)
}

ranperm <- function(X, N) order(runif(N))
randLHS <- function(n,k) {
  P <- matrix(nrow = n, ncol = k)
  P <- apply(P, 2, ranperm, N = n)
  P <- P - 1 + matrix(runif(n * k), nrow = n, ncol = k)
  return(P/n)
}