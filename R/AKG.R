##' Approximate Knowledge Gradient (AKG)
##' 
##' Evaluation of the Approximate Knowledge Gradient (AKG) criterion.
##' 
##' 
##' @param x the input vector at which one wants to evaluate the criterion
##' @param model a Kriging model of "km" class
##' @param new.noise.var (scalar) noise variance of the future observation.
##' Default value is 0 (noise-free observation).
##' @param type Kriging type: "SK" or "UK"
##' @param envir environment for saving intermediate calculations and reusing
##' them within AKG.grad
##' @return Approximate Knowledge Gradient
##' @author Victor Picheny  
##' 
##' David Ginsbourger 
##' 
##' @references Scott, W., Frazier, P., Powell, W. (2011). 
##' The correlated knowledge gradient for simulation optimization of continuous parameters using gaussian process regression. 
##' \emph{SIAM Journal on Optimization}, 21(3), 996-1026.
##' @examples
##' 
##' ##########################################################################
##' ###    AKG SURFACE ASSOCIATED WITH AN ORDINARY KRIGING MODEL          ####
##' ### OF THE BRANIN FUNCTION KNOWN AT A 12-POINT LATIN HYPERCUBE DESIGN ####
##' ##########################################################################
##' set.seed(421)
##' # Set test problem parameters
##' doe.size <- 12
##' dim <- 2
##' test.function <- get("branin2")
##' lower <- rep(0,1,dim)
##' upper <- rep(1,1,dim)
##' noise.var <- 0.2
##' 
##' # Generate DOE and response
##' doe <- as.data.frame(matrix(runif(doe.size*dim),doe.size))
##' y.tilde <- rep(0, 1, doe.size)
##' for (i in 1:doe.size)  {
##'   y.tilde[i] <- test.function(doe[i,]) + sqrt(noise.var)*rnorm(n=1)
##' }
##' y.tilde <- as.numeric(y.tilde)
##' 
##' # Create kriging model
##' model <- km(y~1, design=doe, response=data.frame(y=y.tilde),
##'             covtype="gauss", noise.var=rep(noise.var,1,doe.size), 
##'             lower=rep(.1,dim), upper=rep(1,dim), control=list(trace=FALSE))
##' 
##' # Compute actual function and criterion on a grid
##' n.grid <- 12 # Change to 21 for a nicer picture
##' x.grid <- y.grid <- seq(0,1,length=n.grid)
##' design.grid <- expand.grid(x.grid, y.grid)
##' nt <- nrow(design.grid)
##' 
##' crit.grid <- apply(design.grid, 1, AKG, model=model, new.noise.var=noise.var)
##' func.grid <- apply(design.grid, 1, test.function)
##' 
##' # Compute kriging mean and variance on a grid
##' names(design.grid) <- c("V1","V2")
##' pred <- predict.km(model, newdata=design.grid, type="UK")
##' mk.grid <- pred$m
##' sk.grid <- pred$sd
##' 
##' # Plot actual function
##' z.grid <- matrix(func.grid, n.grid, n.grid)
##' filled.contour(x.grid,y.grid, z.grid, nlevels=50, color = topo.colors,
##'                plot.axes = {title("Actual function");
##'                             points(model@@X[,1],model@@X[,2],pch=17,col="blue"); 
##'                             axis(1); axis(2)})
##' 
##' # Plot Kriging mean
##' z.grid <- matrix(mk.grid, n.grid, n.grid)
##' filled.contour(x.grid,y.grid, z.grid, nlevels=50, color = topo.colors,
##'                plot.axes = {title("Kriging mean");
##'                             points(model@@X[,1],model@@X[,2],pch=17,col="blue"); 
##'                             axis(1); axis(2)})
##' 
##' # Plot Kriging variance
##' z.grid <- matrix(sk.grid^2, n.grid, n.grid)
##' filled.contour(x.grid,y.grid, z.grid, nlevels=50, color = topo.colors,
##'                plot.axes = {title("Kriging variance");
##'                             points(model@@X[,1],model@@X[,2],pch=17,col="blue"); 
##'                             axis(1); axis(2)})
##' 
##' # Plot AKG criterion
##' z.grid <- matrix(crit.grid, n.grid, n.grid)
##' filled.contour(x.grid,y.grid, z.grid, nlevels=50, color = topo.colors,
##'                plot.axes = {title("AKG");
##'                             points(model@@X[,1],model@@X[,2],pch=17,col="blue"); 
##'                             axis(1); axis(2)})
##' 
##' @export
AKG <- function(x, model, new.noise.var=0, type = "UK", envir=NULL)
{

  ######### Convert x in proper format(s) ####
  d <- length(x)
  if (d != model@d){ stop("x does not have the right size") }
  newdata.num <- as.numeric(x)
  newdata <- data.frame(t(newdata.num))
  colnames(newdata) = colnames(model@X)
  tau2.new <- new.noise.var

  ######### Compute prediction at x #########
  predx <- predict.km(model, newdata=newdata, type=type, checkNames = FALSE)
  mk.x <- predx$mean
  sk.x <- predx$sd
  c.x  <- predx$c
  V.x <- predx$Tinv.c
  T <- model@T
  z <- model@z
  U <- model@M
  F.x <- model.matrix(model@trend.formula, data=newdata)

  if (sk.x < sqrt(model@covariance@sd2)/1e6 || model@covariance@sd2 < 1e-20)
  { AKG <- 0
    tuuinv<- mk.X<- V.X<- mu.x<- cn<- sQ<- Isort<- Iremove<- A1<- at<- bt<- ct <- NULL
  } else
  {
    ######### Compute prediction at model@X #########
    predX <- predict.km(model, newdata=model@X, type=type, checkNames = FALSE)
    mk.X <- predX$mean
    V.X <- predX$Tinv.c
    F.X <- model@F
    
    ######### Compute m_min #########
    m_min <- min(c(mk.X,mk.x))
    
    ######### Compute cn #########
    if (type=="UK")
    { tuuinv <- solve(t(U)%*%U)
      mu.x <- (F.X - t(V.X)%*%U)%*% tuuinv %*% t(F.x - t(V.x)%*%U)
    } else
    { tuuinv <- mu.x <- 0}
    
    cn <- c.x - t(V.X)%*%V.x + mu.x
    cn <- c(cn, sk.x^2)
    
    ######### Compute a and b #########
    A <- c(mk.X, mk.x)
    B <- cn / sqrt(tau2.new + sk.x^2)
    sQ <- B[length(B)]
    ######### Careful: the AKG is written for MAXIMIZATION #########
    ######### Minus signs have been added where necessary ##########
    A <- -A
    
    ######## Sort and reduce A and B ####################
    nobs <- model@n
    Isort <- order(x=B,y=A)
    b <- B[Isort]
    a <- A[Isort]
    
    Iremove <- numeric()
    for (i in 1:(nobs))
    {
      if (b[i+1] == b[i])
      {  Iremove <- c(Iremove, i)}
    }
    
    if (length(Iremove) > 0)
    {  b <- b[-Iremove]
       a <- a[-Iremove]}
    nobs <- length(a)-1
    
    C <- rep(0, nobs+2)
    C[1] <- -1e36
    C[length(C)] <- 1e36
    A1 <- 0
    
    for (k in 2:(nobs+1))
    {
      nondom <- 1
      if (k == nobs+1)
      {  nondom <- 1
      } else if ( (a[k+1] >= a[k]) && (b[k] == b[k+1]) ) #((b[k] - b[k+1])^2 < 1e-18)
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
    at <- a[A1+1]
    bt <- b[A1+1]
    ct <- C[c(1, A1+2)]
    
    ######### AGK ##########
    maxNew <- 0       
    for (k in 1:length(at))
    {  maxNew <- maxNew + at[k]*(pnorm(ct[k+1])-pnorm(ct[k])) + bt[k]*(dnorm(ct[k]) - dnorm(ct[k+1]))}
    
    AKG <- maxNew - (-m_min)
  }

  if (!is.null(envir)) {
    assign("mk.x", c.x, envir=envir) 
    assign("c.x", c.x, envir=envir)
    assign("V.x", V.x, envir=envir)
    assign("sk.x", sk.x, envir=envir)
    assign("F.x", F.x, envir=envir)
    assign("tuuinv", tuuinv, envir=envir)
    assign("mk.X", mk.X, envir=envir)
    assign("V.X", V.X, envir=envir)
    assign("mu.x", mu.x, envir=envir)
    assign("cn", cn, envir=envir)
    assign("sQ", sQ, envir=envir)
    assign("Isort", Isort, envir=envir)
    assign("Iremove", Iremove, envir=envir)
    assign("A1", A1, envir=envir)
    assign("at", at, envir=envir)
    assign("bt", bt, envir=envir)
    assign("ct", ct, envir=envir)
  }
  
  return(AKG)
}
