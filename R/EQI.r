EQI <- function(x, model, new.noise.var=0, alpha=0.9, q.min=NULL, type = "UK", envir=NULL)
{

d <- length(x)
newdata <- matrix(x, 1, d)

######### Compute q.min if missing #########
# Compute q.min if missing
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
predx <- predict.km(model, newdata=newdata, type=type)
mk.old <- predx$mean
#sk.old <- predx$sd

######### Intermediate values ##########
T <- model@T
z <- model@z
U <- model@M
V <- predx$Tinv.c
c <- predx$c
w <- as.numeric(t(U)%*%U)
sk.old <- sqrt(model@covariance@sd2 - t(V)%*%V + (1 - t(V)%*%U)^2/w)

######### mq and sq ##########
mq <- mk.old + qnorm(alpha) * sqrt((new.noise.var * sk.old^2)/(new.noise.var + sk.old^2))
sq <- sk.old^2/sqrt(new.noise.var + sk.old^2)

######### EQI ##########
xcr <- (q.min - mq)/sq
xcr.prob <- pnorm(xcr)
xcr.dens <- dnorm(xcr)

if (!is.null(envir)) {
assign("mk.old", mk.old, envir=envir)
assign("sk.old", sk.old, envir=envir)
assign("c", predx$c, envir=envir)
assign("V", V, envir=envir)
assign("z", z, envir=envir)
assign("T", T, envir=envir)
assign("U", U, envir=envir)
assign("w", w, envir=envir)
assign("xcr.prob", xcr.prob, envir=envir)
assign("xcr.dens", xcr.dens, envir=envir)
assign("mq", mq, envir=envir)
assign("sq", sq, envir=envir)
}

eqi.alt <- sq*(xcr*xcr.prob + xcr.dens)

return(eqi.alt)
}
