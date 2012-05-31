EQI.grad <- function(x, model, new.noise.var=0, alpha=0.9, q.min=NULL, type = "UK"){

d <- length(x)
newdata <- matrix(x, 1, d)

######### Compute q.min if missing #########
if (is.null(q.min))
{
  pred <- predict.km(model, newdata=model@X, type=type)
  mk <- pred$mean
  sk <- pred$sd
  sk <- sk*sqrt(model@n/(model@n-1))
  qk <- mk + qnorm(alpha)*sk
  q.min <- min(qk)
}

######### Compute prediction at x #########
predx <- predict(model, newdata=newdata, type=type, checkNames = FALSE)
mk.old <- predx$mean

######### Intermediate values ##########
T <- model@T
z <- model@z
U <- model@M
V <- predx$Tinv.c
c <- predx$c
dc <- covVector.dx(x=as.numeric(newdata), X=model@X, object=model@covariance, c=c)
W <- backsolve(t(T), dc, upper.tri=FALSE)
w <- as.numeric(t(U)%*%U)

sk.old <- sqrt(model@covariance@sd2 - t(V)%*%V + (1 - t(V)%*%U)^2/(t(U)%*%U))

######### mq and sq ##########
mq <- mk.old + qnorm(alpha) * sqrt((new.noise.var * sk.old^2)/(new.noise.var + sk.old^2))
sq <- sk.old^2/sqrt(new.noise.var + sk.old^2)

######### Gradient of m, s2 and lambda with the first model ##########
mk.old.grad <- t(W)%*%z
sk2.old.grad <-  -2*( t(W)%*%V + as.numeric(1 - t(V)%*%U)*(t(W)%*%U)/w )

######### Gradients of mq and sq ##########
mq.grad  <- mk.old.grad + as.numeric(qnorm(alpha)*(new.noise.var)^(3/2)/(2*sk.old*(sk.old^2 + new.noise.var)^(3/2)))*sk2.old.grad
sq2.grad <- as.numeric((sk.old^2)*(2*new.noise.var + sk.old^2)/(new.noise.var + sk.old^2)^2)*sk2.old.grad
sq.grad <- sq2.grad/as.numeric(2*sq)

######### Gradient of EQI ##########
xcr <- (q.min - mq)/sq
xcr.prob <- as.numeric(pnorm(xcr))
xcr.dens <- as.numeric(dnorm(xcr))

grad.EQI <- - xcr.prob * mq.grad + xcr.dens * sq.grad

return(grad.EQI)
}
