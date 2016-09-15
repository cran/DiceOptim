##' Analytical expression of the multipoint expected improvement (qEI)
##' criterion
##' 
##' Computes the multipoint expected improvement criterion.
##' 
##' 
##' @param x a matrix representing the set of input points (one row corresponds
##' to one point) where to evaluate the qEI criterion,
##' @param model an object of class \code{km},
##' @param plugin optional scalar: if provided, it replaces the minimum of the
##' current observations,
##' @param type "SK" or "UK" (by default), depending whether uncertainty
##' related to trend estimation has to be taken into account,
##' @param minimization logical specifying if EI is used in minimiziation or in
##' maximization,
##' @param fastCompute if TRUE, a fast approximation method based on a
##' semi-analytic formula is used (see [Marmin 2014] for details),
##' @param eps the value of \emph{epsilon} of the fast computation trick.
##' Relevant only if \code{fastComputation} is TRUE,
##' @param envir an optional environment specifying where to get intermediate
##' values calculated in \code{\link{qEI}}.
##' @return The multipoint Expected Improvement, defined as \deqn{qEI(X_{new})
##' := E[ ( min(Y(X)) - min(Y( X_{new} )) )_{+} | Y(X)=y(X)],}
##' 
##' where \eqn{X} is the current design of experiments, \eqn{ X_{new} } is a
##' new candidate design, and \eqn{Y} is a random process assumed to have
##' generated the objective function \eqn{y}.
##' @author Sebastien Marmin 
##' 
##' Clement Chevalier 
##' 
##' David Ginsbourger 
##' @seealso \code{\link{EI}}
##' @references
##' 
##' C. Chevalier and D. Ginsbourger (2014) Learning and Intelligent
##' Optimization - 7th International Conference, Lion 7, Catania, Italy,
##' January 7-11, 2013, Revised Selected Papers, chapter Fast computation of
##' the multipoint Expected Improvement with applications in batch selection,
##' pages 59-69, Springer.
##' 
##' D. Ginsbourger, R. Le Riche, L. Carraro (2007), A Multipoint Criterion for
##' Deterministic Parallel Global Optimization based on Kriging. The
##' International Conference on Non Convex Programming, 2007.
##' 
##' S. Marmin. Developpements pour l'evaluation et la maximisation du critere
##' d'amelioration esperee multipoint en optimisation globale (2014). Master's
##' thesis, Mines Saint-Etienne (France) and University of Bern (Switzerland).
##' 
##' D. Ginsbourger, R. Le Riche, and L. Carraro. Kriging is well-suited to
##' parallelize optimization (2010), In Lim Meng Hiot, Yew Soon Ong, Yoel
##' Tenne, and Chi-Keong Goh, editors, \emph{Computational Intelligence in
##' Expensive Optimization Problems}, Adaptation Learning and Optimization,
##' pages 131-162. Springer Berlin Heidelberg.
##' 
##' J. Mockus (1988), \emph{Bayesian Approach to Global Optimization}. Kluwer
##' academic publishers.
##' 
##' M. Schonlau (1997), \emph{Computer experiments and global optimization},
##' Ph.D. thesis, University of Waterloo.
##' @keywords models, parallel optimization
##' 
##' @examples
##' 
##' \donttest{
##' set.seed(007)
##' 
##' # Monte-Carlo validation
##' 
##' # a 4-d, 81-points grid design, and the corresponding response
##' d <- 4; n <- 3^d
##' design <- do.call(expand.grid,rep(list(seq(0,1,length=3)),d))
##' names(design) <- paste("x",1:d,sep="")
##' y <- data.frame(apply(design, 1, hartman4))
##' names(y) <- "y"
##' 
##' # learning
##' model <- km(~1, design=design, response=y, control=list(trace=FALSE))
##' 
##' # pick up 10 points sampled from the 1-point expected improvement
##' q <- 10
##' X <- sampleFromEI(model,n=q)
##' # simulation of the minimum of the kriging random vector at X
##' t1 <- proc.time()
##' newdata <- as.data.frame(X)
##' colnames(newdata) <- colnames(model@@X)
##' 
##' krig  <- predict(object=model, newdata=newdata,type="UK",se.compute=TRUE, cov.compute=TRUE)
##' mk <- krig$mean
##' Sigma.q <- krig$cov
##' mychol <- chol(Sigma.q)
##' nsim <- 300000
##' white.noise <- rnorm(n=nsim*q)
##' minYsim <- apply(crossprod(mychol,matrix(white.noise,nrow=q)) + mk,2,min)
##' # simulation of the improvement (minimization)
##' qImprovement <- (min(model@@y)-minYsim)*((min(model@@y)-minYsim) > 0)
##' 
##' # empirical expectation of the improvement and confident interval (95%)
##' eiMC <- mean(qImprovement)
##' confInterv <- c(eiMC - 1.96*sd(qImprovement)*1/sqrt(nsim),eiMC + 1.96*sd(qImprovement)*1/sqrt(nsim))
##' 
##' # MC estimation of the qEI
##' print(eiMC) 
##' t2 <- proc.time()
##' # qEI with analytical formula
##' qEI(X,model,fastCompute= FALSE)
##' t3 <- proc.time()
##' # qEI with fast computation trick
##' qEI(X,model)
##' t4 <- proc.time()
##' t2-t1  # Time of MC computation
##' t3-t2  # Time of normal computation
##' t4-t3  # Time of fast computation
##' }
##' 
##' @importFrom mnormt dmnorm pmnorm
##' @export qEI
qEI <- function(x,model ,plugin = NULL, type = "UK", minimization = TRUE, fastCompute = TRUE, eps = 10^(-5), envir = NULL) {
  d <- model@d
  if (!is.matrix(x)) {
    x <- matrix(x,ncol = d)
  }
  xb <- rbind(model@X,x)
  nExp <- model@n
  x <- unique(round(xb,digits = 8))
  if ((nExp+1)>length(x[,1])) {return (0)}
  x <- matrix(x[(nExp+1):length(x[,1]),],ncol=d)
  if(is.null(plugin) &&  minimization) plugin <- min(model@y)
  if(is.null(plugin) && !minimization) plugin <- max(model@y)

  if (length(x[,1]) == 1) {
    return(EI(x=x,model=model,plugin=plugin,type=type,envir=envir))
  }

  krig  <- predict(object=model, newdata=x,type=type,se.compute=FALSE, cov.compute=TRUE,checkNames = FALSE)

  sigma <- krig$cov
  mu <- krig$mean
  
  q <- length(mu)
  
  pk <- first_term <- second_term <- rep(0,times=q)
  if (!fastCompute) {
    symetric_term <- matrix(0,q,q)
    non_symetric_term <- matrix(0,q,q)
  }
  for(k in 1:q){
    #covariance matrix of the vector Z^(k)
    Sigma_k <- covZk( sigma=sigma, index=k )
    
    #mean of the vector Z^(k)
    mu_k <- mu - mu[k] 
    mu_k[k] <- -mu[k]
    if(minimization) mu_k <- -mu_k
    
    b_k <- rep(0,times=q)
    b_k[k] <- -plugin
    if(minimization) b_k <- -b_k
    pk[k] <- pmnorm(x= b_k - mu_k  ,varcov=Sigma_k,maxpts=q*200)[1]
    
    first_term[k] <- (mu[k] - plugin)*pk[k]
    if(minimization) first_term[k] <- -first_term[k]
    
    if (fastCompute) {
	    second_term[k] <- 1/eps*(pmnorm(x= b_k + eps*Sigma_k[,k]- mu_k  ,varcov=Sigma_k,maxpts=q*200)[1] - pk[k])
    } else {
        for(i in 1:q){
        	non_symetric_term[k,i] <- Sigma_k[i,k]
	       if (i>=k) {
          mik <- mu_k[i]
	         sigma_ii_k <- Sigma_k[i,i]
	            bik <- b_k[i]
	            phi_ik <- dnorm(x=bik, mean=mik , sd= sqrt(sigma_ii_k))
	  
	            #need c.i^(k) and Sigma.i^(k)
	            cik <- get_cik(b = b_k , m = mu_k , sigma = Sigma_k , i=i  )
	            sigmaik <- get_sigmaik( sigma = Sigma_k , i=i )
	            Phi_ik <- pmnorm(x = cik  ,varcov=sigmaik,maxpts=(q-1)*200)[1]
	            symetric_term[k,i] <- phi_ik*Phi_ik
	          }
          }
        }
  }
  if (!fastCompute) {
    symetric_term <- symetric_term + t(symetric_term)
    diag(symetric_term) <- 0.5 * diag(symetric_term) 
    second_term <- sum(symetric_term*non_symetric_term)
  }  
  if (!is.null(envir)) {
    assign("pk",pk,envir=envir)
    if (fastCompute == FALSE) {
      assign("symetric_term",symetric_term,envir = envir)
    }
    assign("kriging.mean", mu, envir = envir)
    assign("kriging.cov", sigma, envir = envir)
    assign("Tinv.c", krig$Tinv.c, envir = envir)
    assign("c", krig$c, envir = envir)
  }
  somFinal <- sum(first_term,second_term)
  return( somFinal )
}



get_sigmaik <- function( sigma , i ){
  
  result <- sigma
  q <- nrow(sigma)
  
  for(u in 1:q){
    for(v in 1:q){
      if((u!=i) && (v!=i)){
        result[u,v] <- sigma[u,v] - sigma[u,i]*sigma[v,i]/sigma[i,i] 
      }else{
        result[u,v] <- 0
      }
    }
  }

  result <- result[-i,-i]

  return(result)
}


get_cik <- function( b , m , sigma , i  ){
  
  sigmai <- sigma[i,] / sigma[i,i]
  
  result <- (b - m) - ( b[i] - m[i] ) * sigmai
  result <- result[-i]
  
  return(result)
  
}



covZk <- function(sigma,index){
  
  result <- sigma
  q <- nrow(sigma)
  
  result[index,index] <- sigma[index,index]
  for(i in 1:q){
    if(i!=index){
      result[index,i] <- result[i,index] <- sigma[index,index] - sigma[index,i]
    }
  }
  for(i in 1:q){
    for(j in i:q){
      if((i!=index)&&(j!=index)){
        result[i,j] <- result[j,i] <- sigma[i,j] + sigma[index,index] - sigma[index,i] - sigma[index,j]
      }
    }
  }
  
  return(result)
  
}
