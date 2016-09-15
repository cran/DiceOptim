##' Modification of the function  \code{\link[KrigInv]{integration_design}} from the package \code{\link[KrigInv]{KrigInv}} to 
##' be usable for SUR-based optimization with constraints.

##' @title Generic function to build integration points (for the SUR criterion)
##' 
##' @param integcontrol Optional list specifying the procedure to build the integration points and weights. 
##'        Many options are possible.\cr 
##'        A) If nothing is specified, 100*d points are chosen using the Sobol sequence. \cr
##'        B) One can directly set the field \code{integration.points} (p * d matrix) for prespecified integration points. 
##'        In this case these integration points and the corresponding vector \code{integration.weights} will be used 
##'        for all the iterations of the algorithm. \cr
##'        C) If the field \code{integration.points} is not set then the integration points are renewed at each iteration. 
##'        In that case one can control the number of integration points \code{n.points} (default: 100*d) and a specific 
##'        distribution \code{distrib}. Possible values for distrib are: "\code{sobol}", "\code{MC}" and "\code{SUR}"
##'         (default: "\code{sobol}"). \cr
##'        C.1) The choice "\code{sobol}" corresponds to integration points chosen with the Sobol sequence in dimension d (uniform weight). \cr
##'        C.2) The choice "\code{MC}" corresponds to points chosen randomly, uniformly on the domain. \cr
##'        C.3) The choice "\code{SUR}" corresponds to importance sampling distributions (unequal weights). \cr
##'         When important sampling procedures are chosen, \code{n.points} points are chosen using importance sampling among a discrete 
##'         set of \code{n.candidates} points (default: \code{n.points}*10) which are distributed according to a distribution \code{init.distrib} 
##'         (default: "\code{sobol}"). Possible values for \code{init.distrib} are the space filling distributions "\code{sobol}" and "\code{MC}" 
##'         or an user defined distribution "\code{spec}". The "\code{sobol}" and "\code{MC}" choices correspond to quasi random and random points 
##'         in the domain. If the "\code{spec}" value is chosen the user must fill in manually the field \code{init.distrib.spec} to specify 
##'         himself a n.candidates * d matrix of points in dimension d.
##' @param lower Vector containing the lower bounds of the design space.
##' @param upper Vector containing the upper bounds of the design space.
##' @param model.fun object of class \code{\link[DiceKriging]{km}} corresponding to the objective functions,
##' or, if the objective function is fast-to-evaluate, a \code{\link[DOlab]{fastfun}} object,
##' @param model.constraint either one or a list of objects of class \code{\link[DiceKriging]{km}}, one for each constraint function,
##' @param equality either \code{FALSE} if all constraints are for inequalities, else a vector of boolean indicating which are equalities 
##' @param critcontrol optional list of parameters (see \code{\link[DOlab]{crit_SUR_cst}}); here only the component \code{tolConstraints} is used.
##' @param min.prob This argument applies only when importance sampling distributions are chosen. 
##'       For numerical reasons we give a minimum probability for a point to
##'       belong to the importance sample. This avoids probabilities equal to zero and importance sampling
##'       weights equal to infinity. In an importance sample of M points, the maximum weight becomes 
##'       \code{1/min.prob * 1/M}.
##'       
##' @author
##' Victor Picheny 
##' 
##' Mickael Binois 
##' 
##' @references 
##' Chevalier C., Picheny V., Ginsbourger D. (2012), 
##' The KrigInv package: An efficient and user-friendly R implementation of Kriging-based inversion algorithms, 
##' \emph{Computational Statistics and Data Analysis}, 71, 1021-1034.
##' 
##' Chevalier C., Bect J., Ginsbourger D., Vazquez E., Picheny V., Richet Y. (2011), 
##' Fast parallel kriging-based stepwise uncertainty reduction with application to the identification of an excursion set,
##' \emph{Technometrics}, 56(4), 455-465.
##' 
##' V. Picheny (2014),
##' A stepwise uncertainty reduction approach to constrained global optimization,
##' \emph{Proceedings of the 17th International Conference on Artificial Intelligence and Statistics}, JMLR W&CP 33, 787-795.
##' 
##' @seealso \code{\link[DOlab]{crit_SUR_cst}} \code{\link[KrigInv]{integration_design}}
##' @return 
##' A list with components:
##' \itemize{
##' \item{\code{integration.points}}{ p x d matrix of p points used for the numerical calculation of integrals}
##' \item{\code{integration.weights}}{ a vector of size p corresponding to the weight of each point. If all the points are equally 
##' weighted, integration.weights is set to NULL}
##' }
##' @importFrom randtoolbox sobol
##' @export
integration_design_cst <- function(integcontrol = NULL, lower, upper, model.fun=NULL,  model.constraint=NULL, equality=FALSE, critcontrol=NULL, min.prob=0.001){
  result <- NULL
  d <- length(lower)
  if (length(lower) != length(upper)) {
    print("Error in integration_Parameters: 'lower' and 'upper' must have the same length")
    return(NULL)
  }
  
  ## Trivial case 1
  if(is.null(integcontrol$integration.points)) {
    n.int.points <- d * 100
    #     integration.points <- t(lower + t(sobol(n = n.int.points, dim = d)) * (upper - lower))
    integration.points <- matrix(rep(lower,n.int.points)+rep(upper-lower,n.int.points)*runif(n.int.points*d),byrow=TRUE,ncol=d)
    if (d == 1) 
      integration.points <- matrix(integration.points, ncol = 1)
    if (!is.null(model.fun)) colnames(integration.points) <- colnames(model.fun@X)
    result$integration.points <- integration.points
    result$integration.weights <- NULL
    return(result)
  }
  #Trivial case 2
  if(!is.null(integcontrol$integration.points)){
    #integration points pre-specified
    #     if(!is.null(model.fun) && d>1) colnames(integcontrol$integration.points) <- colnames(model.fun@X)
    if (class(integcontrol$integration.points)=="data.frame") integcontrol$integration.points <- as.matrix(integcontrol$integration.points)
    result$integration.points  <- integcontrol$integration.points
    result$integration.weights <- integcontrol$integration.weights
    return(result)
  }
  
  #--------------------------------------------------------
  #non trivial cases:
  if(is.null(integcontrol$n.points)) integcontrol$n.points <- d*100
  if(is.null(integcontrol$distrib)) integcontrol$distrib <- "MC"
  
  if(integcontrol$distrib=="sobol"){
    integration.points <- lower+sobol(n=integcontrol$n.points,dim=d)*(upper-lower)
    if(d==1) integration.points <- matrix(integration.points,ncol=1)
    if(!is.null(model.fun)) colnames(integration.points)<- colnames(model.fun@X)
    result$integration.points <- integration.points
    result$integration.weights<-NULL
    return(result)
  }
  
  if(integcontrol$distrib=="MC"){
    integration.points <- lower+matrix(runif(d*integcontrol$n.points),ncol=d)*(upper-lower)
    if(d==1) integration.points <- matrix(integration.points,ncol=1)
    if(!is.null(model.fun)) colnames(integration.points)<- colnames(model.fun@X)
    result$integration.points <- integration.points
    result$integration.weights<-NULL
    return(result)
  }
  
  if(integcontrol$distrib=="SUR"){
    if(is.null(integcontrol$n.candidates)) integcontrol$n.candidates <- integcontrol$n.points*10
    if(is.null(integcontrol$init.distrib)) integcontrol$init.distrib <- "MC"
    
    #generation of the initial candidates points:
    if(integcontrol$init.distrib=="sobol") initial.integration.points <- t(lower+t(sobol(n=integcontrol$n.candidates,dim=d))*(upper-lower))
    if(integcontrol$init.distrib=="MC") initial.integration.points    <- t(lower+t(matrix(runif(d*integcontrol$n.candidates),ncol=d))*(upper-lower))
    if(integcontrol$init.distrib=="spec") initial.integration.points  <- integcontrol$init.distrib.spec
    
    if (d==1) initial.integration.points <- matrix(initial.integration.points,ncol=1)
    
    #prediction on these initial candidate points
    if(is.null(model.fun || model.constraint)){
      print("Error in integration_Parameters: for the 'SUR' importance sampling distribution, 
            you must set the arguments 'model.fun' and 'model.constraint'")
      return(NULL)
    }
    
    #--------------------------------------------------------------
    if(integcontrol$distrib=="SUR"){
      Tau.n <- prob.of.feasible.improvement(integration.points=initial.integration.points, model.fun=model.fun, 
                                            model.constraint=model.constraint, equality=equality, critcontrol=critcontrol, type = "UK")
    }
    #--------------------------------------------------------------
    Tau.n.sum <- sum(Tau.n)
    if(Tau.n.sum==0) Tau.n.sum <- 1
    prob.n <- pmax(Tau.n/Tau.n.sum,min.prob/integcontrol$n.candidates)
    prob.n <- prob.n/sum(prob.n)
    weight.n <- 1/(prob.n*integcontrol$n.candidates*integcontrol$n.points)
    
    prob.n.copy <- c(0,prob.n)
    prob.n.cum  <- cumsum(prob.n.copy)
    
    my.indices <- findInterval(runif(integcontrol$n.points),prob.n.cum,all.inside=TRUE)
    integration.points <- initial.integration.points[my.indices,]
    integration.weights <- weight.n[my.indices]
    
    if(d==1) integration.points <- matrix(integration.points,ncol=1)
    if(integcontrol$n.points==1) integration.points <- matrix(integration.points,ncol=d)
    
    if(!is.null(model.fun)){
      if (length(model.fun) > 1){
        colnames(integration.points) <- colnames(model.fun[[1]]@X)
      } else {
        colnames(integration.points)<- colnames(model.fun@X) 
      }
    }
    result$integration.points  <- integration.points
    result$integration.weights <- integration.weights
    return(result)
  }
}