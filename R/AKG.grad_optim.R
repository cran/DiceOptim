AKG.grad_optim <- function(x, model, new.noise.var=0, type = "UK", envir){

d <- length(x)
newdata <- matrix(x, 1, d)
type="UK"
tau2.new <- new.noise.var

toget <- matrix(c("mk.x", "c.x", "V.x", "sk.x", "mk.X", 
"V.X", "mu.x", "cn", "sQ", "T", "z", "U", "Isort", "Iremove", "A1", "at", "bt", "ct"),1,18)
    	apply(toget, 2, get, envir=envir)

######### Compute prediction at x #########
mk.x <- envir$mk.x  #predx$mean
c.x  <- envir$c.x   #predx$c
V.x  <- envir$V.x   #predx$Tinv.c
sk.x <- envir$sk.x  #sqrt(model@covariance@sd2 - t(V.x)%*%V.x + (1 - t(V.x)%*%U)^2/(t(U)%*%U))
mk.X <- envir$mk.X  #predX$mean
V.X  <- envir$V.X   #predX$Tinv.c
mu.x <- envir$mu.x  #(1 - t(V.x)%*%U)/(t(U)%*%U)
cn   <- envir$cn
sQ   <- envir$sQ    #b[model@n+1]

T <- envir$T #model@T
z <- envir$z #model@z
U <- envir$U #model@M

dc.x <- covVector.dx(object=model@covariance, x=as.numeric(newdata), X=model@X, c=c.x)

######### Gradient of a ##########
a.grad <- matrix(0,model@d,model@n+1)
W <- backsolve(t(T), dc.x, upper.tri=FALSE)
a.grad[,model@n+1] <- t(W)%*%z

######### Gradient of b ##########
b.grad <- matrix(0,model@d, model@n+1)
cn.grad <- matrix(0,model@d, model@n)

sk2.x.grad <- -2*(t(W)%*%V.x + (t(W)%*%U)*as.numeric(mu.x))

for (i in 1:model@n)
{
  V.i <- V.X[,i]
  mu.i <- (1 - t(V.i)%*%U)/(t(U)%*%U)
  cn.grad[,i] <- t(dc.x[i,] - t(V.i + as.numeric(mu.i)*U)%*%W)
  b.grad[,i] <- cn.grad[,i]/sqrt(tau2.new+sk.x^2) - sk2.x.grad*as.numeric(cn[i]/2/(tau2.new+sk.x^2)^(3/2))
}
b.grad[,model@n+1] <- as.numeric(sk.x^2*(2*tau2.new+sk.x^2)/(tau2.new+sk.x^2)^2)*sk2.x.grad/as.numeric(2*sQ)

######### Careful: the AKG is written for MAXIMIZATION #########
######### Minus signs have been added where necessary ##########
Isort <- envir$Isort
Iremove <- envir$Iremove
A1 <- envir$A1

at <- envir$at
bt <- envir$bt
ct <- envir$ct

# print(at)
# print(bt)
# print(ct)
# print(W)

# Compute m_min
m_min <- min(mk.X)
m_min <- min(m_min,mk.x)
I.min <- which.min(c(m_min,mk.x))
m_min.grad <- matrix(0,model@d,1)
if (I.min == 2)
{  m_min.grad <- a.grad[,model@n+1]}

######## Sort and reduce A and B ####################
b.grad <- b.grad[,Isort]	
a.grad <- -a.grad[,Isort]

if (length(Iremove) > 0)
{  b.grad <- b.grad[,-Iremove]
   a.grad <- a.grad[,-Iremove]
}

# Build 'tilde' data
at.grad <- a.grad[,A1+1]
bt.grad <- b.grad[,A1+1]
nt <- length(at)
ct.grad <- matrix(0,model@d,nt+1)

for (i in 1:(nt-1))
{
  ct.grad[,i+1] <- (at.grad[,i]-at.grad[,i+1]) - (at[i]-at[i+1])*(bt.grad[,i+1]-bt.grad[,i])/(bt[i+1]-bt[i])
}

######### AGK ##########
grad.AKG <- matrix(0,model@d,1)     
for (k in 1:nt)
{  
  grad.AKG <- grad.AKG + at.grad[,k]*(pnorm(ct[k+1])-pnorm(ct[k])) + bt.grad[,k]*(dnorm(ct[k])-dnorm(ct[k+1])) + ct.grad[,k+1]*dnorm(ct[k+1])*(at[k]+bt[k]*ct[k+1]) - ct.grad[,k]*dnorm(ct[k])*(at[k]+bt[k]*ct[k])
}

grad.AKG <- grad.AKG - (-m_min.grad)

return(grad.AKG)
}
