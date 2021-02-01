##' Test whether a set of constraints are violated or not, depending on their nature (equality or inequality) and tolerance parameters
##' 
##' @title Test constraints violation (vectorized)
##' @param cst matrix of constraints (one column for each constraint function)
##' @param equality either FALSE or a Boolean vector defining which constraints are treated as equalities
##' @param tolConstraints tolerance (vector) for all constraints. If not provided, set to zero for inequalities and 0.05 for equalities
##' @return A Boolean vector, TRUE if the point if feasible, FALSE if at least one constraint is violated
##' @author Mickael Binois 
##' 
##' Victor Picheny 
test_feas_vec <- function(cst, equality=FALSE, tolConstraints = NULL){
  if(is.null(dim(cst))){
    cst <- matrix(cst, nrow = 1) 
  }
  
  if(is.null(tolConstraints)){
    epsilon <- rep(0, ncol(cst))
    epsilon[equality] <- 0.05 ## Message?
  }else{
    epsilon <- tolConstraints
  }
  
  tmp <- rep(TRUE, nrow(cst))
  
  eps_mat <- matrix(epsilon, nrow = nrow(cst), ncol = ncol(cst), byrow = T)
  
  # other side for equality constraints
  if(any(equality)){
    tmp[apply(cst[,equality,drop=FALSE] + eps_mat[,equality, drop=FALSE], 1, max)<0] <- FALSE
  }
  
  # side for all inequalities
  tmp[apply(cst - eps_mat, 1, max)>0] <- FALSE
  
  return(tmp)
  
}

test_feas <- function(cst, equality=FALSE, tolConstraints = NULL){
  if(is.null(tolConstraints)){
    epsilon <- rep(0, length(cst))
    epsilon[equality] <- 0.05
  }else{
    epsilon <- tolConstraints
  }
  
  # other side for equality constraints
  if(any(equality)){
    if(any(cst[equality] < -epsilon[equality])){
      return(FALSE)
    }
  }
  
  # side for all inequalities
  if (any(cst > epsilon)){
    return(FALSE)
  } else {
    return(TRUE)
  }
}