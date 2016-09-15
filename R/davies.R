#' @useDynLib DOlab davies
MYdavies <- function(q,lambda,h = rep(1,length(lambda)),
	delta = rep(0,length(lambda)),sigma=0,lim=10000,acc=0.0001) {

  r <- length(lambda)
  if (length(h) != r) stop("lambda and h should have the same length!")
  if (length(delta) != r) stop("lambda and delta should have the same length!")
  
  out <- .C("davies",
  				lambdas=as.double(lambda),
  				noncentral=as.double(delta),
  				df=as.integer(h),
  				r=as.integer(r),
  				sigma=as.double(sigma),
  				q=as.double(q),
  				lq=as.integer(length(q)),
  				lim=as.integer(lim),
  				acc=as.double(acc),
  				trace=double(7*length(q)),
  				ifault=integer(length(q)),
  				res=double(length(q)),
				  package = "DOlab")

  ## out$res <- 1 - out$res
  
  return(list(trace=matrix(out$trace, ncol=7, byrow=TRUE),ifault=out$ifault,Qq=out$res))
  
}

