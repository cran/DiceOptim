AKG.grad <- function(x, model, new.noise.var=0, type = "UK"){

d <- length(x)
newdata <- matrix(x, 1, d)
type="UK"
tau2.new <- new.noise.var

######### Compute prediction at x #########
predx <- predict(model, newdata=newdata, type=type, checkNames = FALSE)
mk.x <- predx$mean
c.x  <- predx$c
V.x <- predx$Tinv.c
T <- model@T
z <- model@z
U <- model@M
sk.x <- sqrt(model@covariance@sd2 - t(V.x)%*%V.x + (1 - t(V.x)%*%U)^2/(t(U)%*%U))

dc.x <- covVector.dx(model@covariance, x=as.numeric(newdata), X=model@X, c=c.x)

######### Compute prediction at X #########
predX <- predict.km(model, newdata=model@X, type=type, checkNames = FALSE)
mk.X <- predX$mean
V.X <- predX$Tinv.c

######### Compute cn #########
mu.x <- (1 - t(V.x)%*%U)/(t(U)%*%U)
cn <- rep(0, model@n+1)
for (i in 1:model@n)
{
  V.i <- V.X[,i]
  cn[i] <- c.x[i] - t(V.i)%*%V.x - mu.x*t(V.i)%*%U + mu.x
}
cn[model@n+1] <- sk.x^2

######### Compute a and b #########
a <- c(mk.X, mk.x)
b <- cn / sqrt(tau2.new + sk.x^2)
sQ <- b[model@n+1]

######### Gradient of a ##########
a.grad <- matrix(0,model@d,model@n+1)

W <- backsolve(t(T), dc.x, upper.tri=FALSE)
a.grad[,model@n+1] <- t(W)%*%z

######### Gradient of b ##########
b.grad <- matrix(0,model@d, model@n+1)

sk2.x.grad <- -2*(t(W)%*%V.x + (t(W)%*%U)*as.numeric(mu.x))

cn.grad <- matrix(0,model@d, model@n)
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
A <- -a
B <- b
B.grad <- b.grad
A.grad <- -a.grad

######## Sort and reduce A and B ####################
# Sort by increasing B values
Isort <- order(x=B,y=A)
b <- B[Isort]	
b.grad <- B.grad[,Isort]
a <- A[Isort]	
a.grad <- A.grad[,Isort]

# Remove rows when b[i+1] == b[i]
Iremove <- numeric()
for (i in 1:(model@n))
{ if (b[i+1] == b[i])
  {  Iremove <- c(Iremove, i)}
}
if (length(Iremove) > 0)
{  b <- b[-Iremove]   
   a <- a[-Iremove]
   b.grad <- b.grad[,-Iremove]
   a.grad <- a.grad[,-Iremove]
}

# Initialize loop
nobs <- length(a)-1
C <- rep(0, nobs+2)
C[1] <- -1e36
C[length(C)] <- 1e36
A1 <- 0

# Loop: build array of indices A1
for (k in 2:(nobs+1))
{
  nondom <- 1

  if (k == nobs+1)
  {  nondom <- 1
  } else if ( (a[k+1] >= a[k]) && (b[k] == b[k+1]) )
  {  nondom <- 0}
            
  if (nondom == 1)       
  {
    loopdone <- 0
    count <- 0
    while ( loopdone == 0 && count < 1e3 )
    {
      count <- count + 1
      u <- A1[length(A1)] + 1
      C[u+1] <- (a[u]-a[k]) / (b[k] - b[u])
      if ((length(A1) > 1) && (C[u+1] <= C[A1[length(A1)-1]+2]))
      { A1 <- A1[-length(A1)]
      } else
      { A1 <- c(A1, k-1)
        loopdone <- 1
      }
    }
  }
}
# Build 'tilde' data
at <- a[A1+1]
bt <- b[A1+1]
ct <- C[c(1, A1+2)]
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

# Compute m_min
m_min <- min(mk.X)
m_min <- min(m_min,mk.x)
I.min <- which.min(c(m_min,mk.x))

m_min.grad <- matrix(0,model@d,1)
if (I.min == 2)
{  m_min.grad <- A.grad[,model@n+1]}

grad.AKG <- grad.AKG - (m_min.grad)

return(grad.AKG)
}
