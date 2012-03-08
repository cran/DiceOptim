EQI.grad_optim <- function(x, model, new.noise.var=0, alpha=0.9, q.min=NULL, type = "UK", envir){

d <- length(x)
newdata <- matrix(x, 1, d)

toget <- matrix(c("xcr.prob", "xcr.dens", "mk.old", "sk.old", "c", "z", "T", "U", "V", "w", "mq", "sq"),1,12)
    	apply(toget, 2, get, envir=envir)

######### Intermediate values ##########
mk.old <- envir$mk.old #predx$mean
sk.old <- envir$sk.old
T <- envir$T #model@T
z <- envir$z #model@z
U <- envir$U #model@M
V <- envir$V #predx$Tinv.c
c <- envir$c #predx$c
w <- as.numeric(envir$w)

dc <- covVector.dx(x=as.numeric(newdata), X=model@X, object=model@covariance, c=c)
W <- backsolve(t(T), dc, upper.tri=FALSE)

######### Gradient of m, s2 and lambda with the first model ##########
mk.old.grad <- t(W)%*%z
sk2.old.grad <-  -2*( t(W)%*%V + as.numeric(1 - t(V)%*%U)*(t(W)%*%U)/w )

######### mq and sq ##########
mq <- envir$mq
sq <- envir$sq

######### Gradients of mq and sq ##########
mq.grad  <- mk.old.grad + as.numeric(qnorm(alpha)*(new.noise.var^(3/2))/(2*sk.old*(sk.old^2 + new.noise.var)^(3/2)))*sk2.old.grad
sq2.grad <- as.numeric((sk.old^2)*(2*new.noise.var + sk.old^2)/(new.noise.var + sk.old^2)^2)*sk2.old.grad
sq.grad <- sq2.grad/as.numeric(2*sq)

######### Gradient of EQI ##########
# xcr <- envir$xcr #(q.min - mq)/sq
xcr.prob <- as.numeric(envir$xcr.prob) #as.numeric(pnorm(xcr))
xcr.dens <- as.numeric(envir$xcr.dens) #as.numeric(dnorm(xcr))

grad.EQI <- - xcr.prob * mq.grad + xcr.dens * sq.grad

return(grad.EQI)
}
