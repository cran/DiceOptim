prob.of.feasible.improvement <- function(integration.points, model.fun, model.constraint, equality = FALSE, critcontrol = NULL, type = "UK")
{
  n.cst <- length(model.constraint)
  if (is.null(dim(x))) x <- matrix(x, nrow=1)
  if(n.cst == 1 && class(model.constraint) != "list")  model.constraint <- list(model.constraint)
  
  if (is.null(nrow(integration.points))) integration.points <- as.matrix(integration.points)
  
  n.integration.points <- nrow(integration.points)
  mn.X.cst <- sn.X.cst <- matrix(0, n.cst, nrow(integration.points))
  
  for (i in 1:n.cst){
    p.tst <- predict(model.constraint[[i]], newdata=integration.points, type=type, checkNames=FALSE)
    mn.X.cst[i,] <- p.tst$mean
    sn.X.cst[i,] <- p.tst$sd
  }
  
  p.tst <- predict(model.fun, newdata=integration.points, type=type, checkNames=FALSE)
  mn.X.obj <- p.tst$mean
  sn.X.obj <- p.tst$sd
  
  #--- Check equality setting ---------------------------------------------------
  if(is.null(critcontrol$tolConstraints)){
    tolvec <- rep(0, n.cst)
    if(any(equality != FALSE)) tolvec[equality] <- 0.05
    else                       tolvec <- critcontrol$tolConstraints
  }
  
  #--- Compute current best ----------------------------------------------------
  obs.cst <- c()
  for (i in 1:n.cst) obs.cst <- cbind(obs.cst, model.constraint[[i]]@y)
  feasibility <- test_feas_vec(cst=obs.cst, equality=equality, tolConstraints = critcontrol$tolConstraints)
  if (any(feasibility)) obj.min <- min(model.fun@y[feasibility])
  else                  obj.min <- Inf

  #--- Update model.fun ---------------------------------------------------
  krig.obj <- predict(object=model.fun, newdata=x, type=type, se.compute=TRUE, cov.compute=FALSE, checkNames=FALSE)   
  mk.obj <- krig.obj$mean
  sk.obj <- krig.obj$sd
  if (class(model.fun)=="km") {
    if (obj.min!=Inf) {
      obj.min.tilde <- (obj.min - mn.X.obj) / sn.X.obj
      p.obj.min.tilde <- pnorm(obj.min.tilde)
    } else {
      p.obj.min.tilde <- rep(1, n.integration.points)
    }
  } else {
    p.obj.min.tilde <- as.numeric(mn.X.obj <= obj.min)
  }
  
  #--- Update model.constraint ----------------------------------------------------
  Ttilde <- Ttildem <- Ttildep <- matrix(0, n.cst, n.integration.points)
  
  for (i in 1:n.cst) {
    if (any(equality != FALSE)) {
      # Equality setting
      if (equality[i]) {
        Ttildep[i,]  <- (tolvec[i] - mn.X.cst[i,] )/sn.X.cst[i,]
        Ttildem[i,]  <- (-tolvec[i] - mn.X.cst[i,] )/sn.X.cst[i,]
      } else {
        Ttilde[i,]  <- (0 - mn.X.cst[i,] )/sn.X.cst[i,]
      }
    } else {
      Ttilde[i,]  <- (0 - mn.X.cst[i,] )/sn.X.cst[i,] 
    }
  }
  
  #--- Compute criterion -------------------------------------------------
  if (any(equality != FALSE)) {
    p.f <- 1
    for (i in 1:n.cst) {
      if (equality[i]) {
        p.f <- p.f*(pnorm(Ttildep[i,]) - pnorm(Ttildem[i,]))
      } else {
        p.f <- p.f*pnorm(Ttilde[i,])
      }
    }
    oldcrit <- p.obj.min.tilde*p.f
  } else {
    oldcrit <- p.obj.min.tilde*apply(pnorm(Ttilde),2,prod)
  }
  return(oldcrit)
}