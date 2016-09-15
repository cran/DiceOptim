##' Sequential multipoint Expected improvement (qEI) maximizations and model
##' re-estimation
##' 
##' Executes \code{nsteps} iterations of the multipoint EGO method to an object
##' of class \code{\link[DiceKriging]{km}}.  At each step, a kriging model
##' (including covariance parameters) is re-estimated based on the initial
##' design points plus the points visited during all previous iterations; then
##' a new batch of points is obtained by maximizing the multipoint Expected
##' Improvement criterion (\code{\link{qEI}}).
##' 
##' The parameters of list \code{optimcontrol} are :
##' 
##' - \code{optimcontrol$method} : "BFGS" (default), "genoud" ; a string
##' specifying the method used to maximize the criterion (irrelevant when
##' \code{crit} is "CL" because this method always uses genoud),
##' 
##' - when \code{crit="CL"} :
##' 
##' + \code{optimcontrol$parinit} : optional matrix of initial values (must
##' have model@@d columns, the number of rows is not constrained),
##' 
##' + \code{optimcontrol$L} : "max", "min", "mean" or a scalar value specifying
##' the liar ; "min" takes \code{model@@min}, "max" takes \code{model@@max},
##' "mean" takes the prediction of the model ; When L is \code{NULL}, "min" is
##' taken if \code{minimization==TRUE}, else it is "max".
##' 
##' + The parameters of function \code{\link{genoud}}. Main parameters are :
##' \code{"pop.size"} (default : [N=3*2^model@@d for dim<6 and N=32*model@@d
##' otherwise]), \code{"max.generations"} (default : 12),
##' \code{"wait.generations"} (default : 2) and \code{"BFGSburnin"} (default :
##' 2).
##' 
##' - when \code{optimcontrol$method = "BFGS"} :
##' 
##' + \code{optimcontrol$nStarts} (default : 4),
##' 
##' + \code{optimcontrol$fastCompute} : if TRUE (default), a fast approximation
##' method based on a semi-analytic formula is used, see [Marmin 2014] for
##' details,
##' 
##' + \code{optimcontrol$samplingFun} : a function which sample a batch of
##' starting point (default : \code{\link{sampleFromEI}}),
##' 
##' + \code{optimcontrol$parinit} : optional 3d-array of initial (or candidate)
##' batches (for all \code{k}, parinit[,,k] is a matrix of size
##' \code{npoints*model@@d} representing one batch). The number of initial
##' batches (length(parinit[1,1,])) is not contrained and does not have to be
##' equal to \code{nStarts}. If there is too few initial batches for
##' \code{nStarts}, missing batches are drawn with \code{samplingFun} (default
##' : \code{NULL}),
##' 
##' - when \code{optimcontrol$method = "genoud"} :
##' 
##' + \code{optimcontrol$fastCompute} : if TRUE (default), a fast approximation
##' method based on a semi-analytic formula is used, see [Marmin 2014] for
##' details,
##' 
##' + \code{optimcontrol$parinit} : optional matrix of candidate starting
##' points (one row corresponds to one point),
##' 
##' + The parameters of the \code{\link{genoud}} function. Main parameters are
##' \code{"pop.size"} (default : \code{[50*(model@@d)*(npoints)]}),
##' \code{"max.generations"} (default : 5), \code{"wait.generations"} (default
##' : 2), \code{"BFGSburnin"} (default : 2).
##' 
##' 
##' @param fun the objective function to be optimized,
##' @param model an object of class \code{\link[DiceKriging]{km}} ,
##' @param npoints an integer repesenting the desired batchsize,
##' @param nsteps an integer representing the desired number of iterations,
##' @param lower vector of lower bounds for the variables to be optimized over,
##' @param upper vector of upper bounds for the variables to be optimized over,
##' @param crit "exact", "CL" : a string specifying the criterion used. "exact"
##' triggers the maximization of the multipoint expected improvement at each
##' iteration (see \code{\link{max_qEI}}), "CL" applies the Constant Liar
##' heuristic,
##' @param minimization logical specifying if we want to minimize or maximize
##' \code{fun},
##' @param optimcontrol an optional list of control parameters for the qEI
##' optimization (see details or \code{\link{max_qEI}}),
##' @param cov.reestim optional boolean specifying if the kriging
##' hyperparameters should be re-estimated at each iteration,
##' @param ...  optional arguments for \code{fun}.
##' @return A list with components:
##' 
##' \item{par}{a data frame representing the additional points visited during
##' the algorithm,}
##' 
##' \item{value}{a data frame representing the response values at the points
##' given in \code{par},}
##' 
##' \item{npoints}{an integer representing the number of parallel
##' computations,}
##' 
##' \item{nsteps}{an integer representing the desired number of iterations
##' (given in argument),}
##' 
##' \item{lastmodel}{an object of class \code{\link[DiceKriging]{km}}
##' corresponding to the last kriging model fitted,}
##' 
##' \item{history}{a vector of size \code{nsteps} representing the current
##' known optimum at each step.}
##' @author Sebastien Marmin 
##' 
##' Clement Chevalier 
##' 
##' David Ginsbourger 
##' @seealso \code{\link{qEI}}, \code{\link{max_qEI}}, \code{\link{qEI.grad}}
##' @references
##' 
##' C. Chevalier and D. Ginsbourger (2014) Learning and Intelligent
##' Optimization - 7th International Conference, Lion 7, Catania, Italy,
##' January 7-11, 2013, Revised Selected Papers, chapter Fast computation of
##' the multipoint Expected Improvement with applications in batch selection,
##' pages 59-69, Springer.
##' 
##' D. Ginsbourger, R. Le Riche, L. Carraro (2007), A Multipoint Criterion for
##' Deterministic Parallel Global Optimization based on Kriging. The
##' International Conference on Non Convex Programming, 2007.
##' 
##' D. Ginsbourger, R. Le Riche, and L. Carraro. Kriging is well-suited to
##' parallelize optimization (2010), In Lim Meng Hiot, Yew Soon Ong, Yoel
##' Tenne, and Chi-Keong Goh, editors, \emph{Computational Intelligence in
##' Expensive Optimization Problems}, Adaptation Learning and Optimization,
##' pages 131-162. Springer Berlin Heidelberg.
##' 
##' S. Marmin. Developpements pour l'evaluation et la maximisation du critere
##' d'amelioration esperee multipoint en optimisation globale (2014). Master's
##' thesis, Mines Saint-Etienne (France) and University of Bern (Switzerland).
##' 
##' J. Mockus (1988), \emph{Bayesian Approach to Global Optimization}. Kluwer
##' academic publishers.
##' 
##' M. Schonlau (1997), \emph{Computer experiments and global optimization},
##' Ph.D. thesis, University of Waterloo.
##' @keywords optimize
##' @examples
##' 
##' 
##' \dontrun{
##' set.seed(123)
##' #####################################################
##' ### 4 ITERATIONS OF EGO ON THE BRANIN FUNCTION,   ###
##' ### STARTING FROM A 9-POINTS FACTORIAL DESIGN     ###
##' #####################################################
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
##' nsteps <- 4 # Was 10, reduced to 4 for speeding up execution
##' lower <- rep(0,d) 
##' upper <- rep(1,d)
##' npoints <- 3 # The batchsize
##' oEGO <- qEGO.nsteps(model = fitted.model1, branin, npoints = npoints, nsteps = nsteps,
##' crit="exact", lower, upper, optimcontrol = NULL)
##' print(oEGO$par)
##' print(oEGO$value)
##' plot(c(1:nsteps),oEGO$history,xlab='step',ylab='Current known minimum')
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
##' text(oEGO$par[,1], oEGO$par[,2], labels=c(tcrossprod(rep(1,npoints),1:nsteps)), pos=3)
##' }
##' 
##' @export qEGO.nsteps
qEGO.nsteps <- function (fun, model, npoints, nsteps, lower = rep(0, model@d), upper = rep(1, model@d), crit = "exact", minimization = TRUE,  optimcontrol = NULL, cov.reestim=TRUE, ...) 
{
    history <- rep(NaN,nsteps)
    d <- model@d
    n <- nrow(model@X)
    parinit <- optimcontrol$parinit
    for (i in 1:nsteps) {
	    print(paste("step :",i,"/",nsteps))
	    history[i] <- min(model@y)*minimization + max(model@y)*(!minimization)
	    res <- max_qEI(model=model, npoints = npoints, lower = lower, upper = upper,crit=crit, minimization = minimization, optimcontrol = optimcontrol)
	    model <- update(object = model, newX = res$par, newy=as.matrix(apply(res$par,1, fun, ...), npoints), newX.alreadyExist=FALSE,
                                  cov.reestim = cov.reestim, kmcontrol = list(control = list(trace = FALSE)))
    }
    return(list(par = model@X[(n + 1):(n + nsteps * npoints), 
        , drop = FALSE], value = model@y[(n + 1):(n + nsteps * 
        npoints), , drop = FALSE], npoints = npoints, nsteps = nsteps, 
        lastmodel = model,history=history))
}

