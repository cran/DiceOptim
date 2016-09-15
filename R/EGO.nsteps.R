##' Sequential EI maximization and model re-estimation, with a number of
##' iterations fixed in advance by the user
##' 
##' Executes \emph{nsteps} iterations of the EGO method to an object of class
##' \code{\link[DiceKriging]{km}}.  At each step, a kriging model is
##' re-estimated (including covariance parameters re-estimation) based on the
##' initial design points plus the points visited during all previous
##' iterations; then a new point is obtained by maximizing the Expected
##' Improvement criterion (\code{\link{EI}}).
##' 
##' 
##' @param model an object of class \code{\link[DiceKriging]{km}} ,
##' @param fun the objective function to be minimized,
##' @param nsteps an integer representing the desired number of iterations,
##' @param lower vector of lower bounds for the variables to be optimized over,
##' @param upper vector of upper bounds for the variables to be optimized over,
##' @param parinit optional vector of initial values for the variables to be
##' optimized over,
##' @param control an optional list of control parameters for optimization. One
##' can control
##' 
##' \code{"pop.size"} (default : [4+3*log(nb of variables)]),
##' 
##' \code{"max.generations"} (default :5),
##' 
##' \code{"wait.generations"} (default :2),
##' 
##' \code{"BFGSburnin"} (default :0),
##' 
##' of the function \code{\link[rgenoud]{genoud}}.
##' @param kmcontrol an optional list representing the control variables for
##' the re-estimation of the kriging model.  The items are the same as in
##' \code{\link[DiceKriging]{km}} :
##' 
##' \code{penalty}, \code{optim.method}, \code{parinit}, \code{control}.
##' 
##' The default values are those contained in \code{model}, typically
##' corresponding to the variables used in \code{\link[DiceKriging]{km}} to
##' estimate a kriging model from the initial design points.
##' @return A list with components:
##' 
##' \item{par}{a data frame representing the additional points visited during
##' the algorithm,}
##' 
##' \item{value}{a data frame representing the response values at the points
##' given in \code{par},}
##' 
##' \item{npoints}{an integer representing the number of parallel computations
##' (=1 here),}
##' 
##' \item{nsteps}{an integer representing the desired number of iterations
##' (given in argument),}
##' 
##' \item{lastmodel}{an object of class \code{\link[DiceKriging]{km}}
##' corresponding to the last kriging model fitted.}
##' @note Most EGO-like methods (EI algorithms) usually work with Ordinary
##' Kriging (constant trend), by maximization of the expected improvement.
##' Here, the EI maximization is also possible with any linear trend. However,
##' note that the optimization may perform much faster and better when the
##' trend is a constant since it is the only case where the analytical gradient
##' is available.
##' 
##' For more details on \code{kmcontrol}, see the documentation of
##' \code{\link[DiceKriging]{km}}.
##' @author David Ginsbourger 
##' 
##' Olivier Roustant 
##' 
##' @seealso \code{\link{EI}}, \code{\link{max_EI}}, \code{\link{EI.grad}}
##' @references
##' 
##' D.R. Jones, M. Schonlau, and W.J. Welch (1998), Efficient global
##' optimization of expensive black-box functions, \emph{Journal of Global
##' Optimization}, 13, 455-492.
##' 
##' J. Mockus (1988), \emph{Bayesian Approach to Global Optimization}. Kluwer
##' academic publishers.
##' 
##' T.J. Santner, B.J. Williams, and W.J. Notz (2003), \emph{The design and
##' analysis of computer experiments}, Springer.
##' 
##' M. Schonlau (1997), \emph{Computer experiments and global optimization},
##' Ph.D. thesis, University of Waterloo.
##' @keywords optimize
##' @examples
##' 
##' 
##' set.seed(123)
##' ###############################################################
##' ### 	10 ITERATIONS OF EGO ON THE BRANIN FUNCTION, 	   ####
##' ###	 STARTING FROM A 9-POINTS FACTORIAL DESIGN         ####
##' ###############################################################
##' 
##' # a 9-points factorial design, and the corresponding response
##' d <- 2
##' n <- 9
##' design.fact <- expand.grid(seq(0,1,length=3), seq(0,1,length=3)) 
##' names(design.fact)<-c("x1", "x2")
##' design.fact <- data.frame(design.fact) 
##' names(design.fact)<-c("x1", "x2")
##' response.branin <- apply(design.fact, 1, branin)
##' response.branin <- data.frame(response.branin) 
##' names(response.branin) <- "y" 
##' 
##' # model identification
##' fitted.model1 <- km(~1, design=design.fact, response=response.branin, 
##' covtype="gauss", control=list(pop.size=50,trace=FALSE), parinit=c(0.5, 0.5))
##' 
##' # EGO n steps
##' library(rgenoud)
##' nsteps <- 5 # Was 10, reduced to 5 for speeding up compilation
##' lower <- rep(0,d) 
##' upper <- rep(1,d)     
##' oEGO <- EGO.nsteps(model=fitted.model1, fun=branin, nsteps=nsteps, 
##' lower=lower, upper=upper, control=list(pop.size=20, BFGSburnin=2))
##' print(oEGO$par)
##' print(oEGO$value)
##' 
##' # graphics
##' n.grid <- 15 # Was 20, reduced to 15 for speeding up compilation
##' x.grid <- y.grid <- seq(0,1,length=n.grid)
##' design.grid <- expand.grid(x.grid, y.grid)
##' response.grid <- apply(design.grid, 1, branin)
##' z.grid <- matrix(response.grid, n.grid, n.grid)
##' contour(x.grid, y.grid, z.grid, 40)
##' title("Branin function")
##' points(design.fact[,1], design.fact[,2], pch=17, col="blue")
##' points(oEGO$par, pch=19, col="red")
##' text(oEGO$par[,1], oEGO$par[,2], labels=1:nsteps, pos=3)
##' 
##' 
##' ###############################################################
##' ### 	20 ITERATIONS OF EGO ON THE GOLDSTEIN-PRICE,  	   ####
##' ###	 STARTING FROM A 9-POINTS FACTORIAL DESIGN   	   ####
##' ###############################################################
##' \dontrun{
##' # a 9-points factorial design, and the corresponding response
##' d <- 2 
##' n <- 9
##' design.fact <- expand.grid(seq(0,1,length=3), seq(0,1,length=3)) 
##' names(design.fact)<-c("x1", "x2")
##' design.fact <- data.frame(design.fact) 
##' names(design.fact)<-c("x1", "x2")
##' response.goldsteinPrice <- apply(design.fact, 1, goldsteinPrice)
##' response.goldsteinPrice <- data.frame(response.goldsteinPrice) 
##' names(response.goldsteinPrice) <- "y" 
##' 
##' # model identification
##' fitted.model1 <- km(~1, design=design.fact, response=response.goldsteinPrice, 
##' covtype="gauss", control=list(pop.size=50, max.generations=50, 
##' wait.generations=5, BFGSburnin=10,trace=FALSE), parinit=c(0.5, 0.5), optim.method="BFGS")
##' 
##' # EGO n steps
##' library(rgenoud)
##' nsteps <- 10 # Was 20, reduced to 10 for speeding up compilation
##' lower <- rep(0,d) 
##' upper <- rep(1,d)     
##' oEGO <- EGO.nsteps(model=fitted.model1, fun=goldsteinPrice, nsteps=nsteps, 
##' lower, upper, control=list(pop.size=20, BFGSburnin=2))
##' print(oEGO$par)
##' print(oEGO$value)
##' 
##' # graphics
##' n.grid <- 15 # Was 20, reduced to 15 for speeding up compilation
##' x.grid <- y.grid <- seq(0,1,length=n.grid)
##' design.grid <- expand.grid(x.grid, y.grid)
##' response.grid <- apply(design.grid, 1, goldsteinPrice)
##' z.grid <- matrix(response.grid, n.grid, n.grid)
##' contour(x.grid, y.grid, z.grid, 40)
##' title("Goldstein-Price Function")
##' points(design.fact[,1], design.fact[,2], pch=17, col="blue")
##' points(oEGO$par, pch=19, col="red")
##' text(oEGO$par[,1], oEGO$par[,2], labels=1:nsteps, pos=3)
##' }
##' 
##' #######################################################################
##' ### 	nsteps ITERATIONS OF EGO ON THE HARTMAN6 FUNCTION,  	   ####
##' ###	  STARTING FROM A 10-POINTS UNIFORM DESIGN      	   ####
##' #######################################################################
##' 
##' \dontrun{
##' fonction<-hartman6
##' data(mydata)
##' a <- mydata
##' nb<-10
##' nsteps <- 3 # Maybe be changed to a larger value 
##' x1<-a[[1]][1:nb];x2<-a[[2]][1:nb];x3<-a[[3]][1:nb]
##' x4<-a[[4]][1:nb];x5<-a[[5]][1:nb];x6<-a[[6]][1:nb]
##' design <- data.frame(cbind(x1,x2,x3,x4,x5,x6)) 
##' names(design)<-c("x1", "x2","x3","x4","x5","x6")
##' n <- nrow(design)
##' 
##' response <- data.frame(q=apply(design,1,fonction)) 
##' names(response) <- "y" 
##' fitted.model1 <- km(~1, design=design, response=response, covtype="gauss", 
##' control=list(pop.size=50, max.generations=20, wait.generations=5, BFGSburnin=5,
##' trace=FALSE), optim.method="gen", parinit=rep(0.8,6))
##' 
##' res.nsteps <- EGO.nsteps(model=fitted.model1, fun=fonction, nsteps=nsteps, 
##' lower=rep(0,6), upper=rep(1,6), parinit=rep(0.5,6), control=list(pop.size=50, 
##' max.generations=20, wait.generations=5, BFGSburnin=5), kmcontrol=NULL)
##' print(res.nsteps)
##' plot(res.nsteps$value,type="l")
##' }
##' 
##' @export EGO.nsteps
EGO.nsteps <-function(model, fun, nsteps, lower, upper, parinit=NULL, control=NULL, kmcontrol=NULL) {

	n <- nrow(model@X)
	
	if (is.null(kmcontrol$penalty)) kmcontrol$penalty <- model@penalty
	if (length(model@penalty==0)) kmcontrol$penalty <- NULL 

	if (is.null(kmcontrol$optim.method)) kmcontrol$optim.method <- model@optim.method 
	if (is.null(kmcontrol$parinit)) kmcontrol$parinit <- model@parinit
	if (is.null(kmcontrol$control)) kmcontrol$control <- model@control
		

	for (i in 1:nsteps) {
		oEGO<-max_EI(model=model, lower=lower, upper=upper, parinit=parinit, control=control)
		
		model@X<-rbind(model@X, oEGO$par)
		model@y<-rbind(model@y, fun(t(oEGO$par)))
		
		kmcontrol$parinit <- covparam2vect(model@covariance)
    kmcontrol$control$trace=FALSE
		
		if (model@param.estim) {
			model <- km(formula=model@trend.formula, design=model@X, response=model@y, 
		   		 covtype=model@covariance@name, lower=model@lower, upper=model@upper, 
                 nugget=NULL, penalty=kmcontrol$penalty, optim.method=kmcontrol$optim.method, 
		    	parinit=kmcontrol$parinit, control=kmcontrol$control, gr=model@gr, iso=is(model@covariance,"covIso"))
		} else {
			coef.cov <- covparam2vect(model@covariance)
			model <- km(formula=model@trend.formula, design=model@X, response=model@y, 
		   		 covtype=model@covariance@name, coef.trend=model@trend.coef, coef.cov=coef.cov, coef.var=model@covariance@sd2, nugget=NULL, iso=is(model@covariance,"covIso"))
		}
	}
	
	return(list(par=model@X[(n+1):(n+nsteps),, drop=FALSE], 
               value=model@y[(n+1):(n+nsteps),, drop=FALSE], npoints=1, nsteps=nsteps, lastmodel=model))

}
