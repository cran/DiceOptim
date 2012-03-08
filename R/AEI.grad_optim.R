AEI.grad_optim <- function(x, model, new.noise.var=0, y.min=NULL, type = "UK", envir){

d <- length(x)
x <- matrix(x, 1, d)

# Compute y.min if missing
if (is.null(y.min))
{
  pred <- predict.km(model, newdata=model@X, type="UK")
  mk <- pred$mean
  sk <- pred$sd
  qk <- mk + qnorm(0.75)*sk
  y.min <- mk[which.min(qk)]
}

toget <- matrix(c("xcr", "xcr.prob", "xcr.dens", "c", "Tinv.c", "mk","sk"),1,7)
    	apply(toget, 2, get, envir=envir)

# Model parameters
#pred <- envir$pred #predict.km(model, x, type="UK")
mk <- envir$mk #pred$mean
sk <- envir$sk #pred$sd

T <- model@T
X <- model@X
z <- model@z
covStruct <- model@covariance
c <- envir$c #pred$c
v <- envir$Tinv.c #pred$Tinv.c
u <- model@M

# AEI
xcr <- envir$xcr #(y.min - mk)/sk
xcr.prob <- envir$xcr.prob #pnorm(xcr)
xcr.dens <- envir$xcr.dens #dnorm(xcr)

ei.val <- (y.min - mk) * xcr.prob + sk * xcr.dens
pen <- (1- sqrt(new.noise.var)/sqrt(new.noise.var + sk^2))

aei.val <- ei.val * pen

# EI gradient
#dc <- covVector.dx(as.numeric(x), X, covStruct, c)
dc <- covVector.dx(x=as.numeric(x), X=X, object=covStruct, c=c)
W <- backsolve(t(T), dc, upper.tri=FALSE)

mk.grad <- t(W)%*%z
sk2.grad <-  -2*( t(W)%*%v + (t(W)%*%u)* as.numeric(1 - t(v)%*%u)/as.numeric(t(u)%*%u) )
sk.grad <- sk2.grad /(2*sk)

#kriging.sd.grad <- - t(W)%*%v + (t(W)%*%u)* as.numeric((1-t(v)%*%u)/(t(u)%*%u)) / sk

ei.grad <- - mk.grad * xcr.prob + sk.grad * xcr.dens

# Penalization gradient
pen.grad <- - ei.val * (1 - sqrt(new.noise.var))*sk*sk.grad*(new.noise.var + sk^2)^(-3/2)

# AEI gradient
aei.grad.val <- ei.grad*pen + ei.val*pen.grad
return(aei.grad.val)
}
