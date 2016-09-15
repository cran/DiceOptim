

##' Kriging-based optimization methods for computer experiments
##' 
##' Sequential and parallel Kriging-based optimization methods relying on
##' expected improvement criteria.
##' 
##' \tabular{ll}{ Package: \tab DiceOptim\cr Type: \tab Package\cr Version:
##' \tab 2.0 \cr Date: \tab July 2016\cr License: \tab GPL-2 | GPL-3\cr }
##' 
##' @name DiceOptim-package
##' @aliases DiceOptim-package DiceOptim
##' @docType package
##' @note This work is a follow-up of DiceOptim 1.0, which was produced within
##' the frame of the DICE (Deep Inside Computer Experiments) Consortium between
##' ARMINES, Renault, EDF, IRSN, ONERA and TOTAL S.A.
##' 
##' The authors would like to thank Yves Deville for his precious advice in R
##' programming and packaging, as well as the DICE members for useful
##' feedbacks, and especially Yann Richet (IRSN) for numerous discussions
##' concerning the user-friendliness of this package.
##' 
##' Package \code{rgenoud} >=5.3.3. is recommended.
##' 
##' Important functions or methods: \tabular{ll}{ \code{EGO.nsteps} \tab Standard Efficient Global Optimization algorithm with a fixed number of iterations (nsteps) \cr 
##' \tab---with model updates including re-estimation of covariance hyperparameters \cr
##' 
##' \code{EI} \tab Expected Improvement criterion (single infill point, noise-free, constraint free problems)\cr
##' 
##' \code{max_EI} \tab Maximization of the EI criterion. No need to specify any objective function \cr
##' 
##' \code{qEI.nsteps} \tab EGO algorithm with batch-sequential (parallel) infill strategy \cr
##' 
##' \code{noisy.optimizer} \tab EGO algorithm for noisy objective functions \cr
##' 
##' \code{EGO.cst} \tab EGO algorithm for (non-linear) constrained problems \cr
##' 
##'\code{easyEGO.cst} \tab User-friendly wrapper for \code{EGO.cst}}
##'
##' @author Victor Picheny (INRA, Castanet-Tolosan, France)
##' 
##' David Ginsbourger (Idiap Research Institute and University of Bern, Switzerland)
##' 
##' Olivier Roustant (Mines Saint-Etienne, France).
##' 
##' with contributions by M. Binois, C. Chevalier, S. Marmin and T. Wagner
##' @references N.A.C. Cressie (1993), \emph{Statistics for spatial data},
##' Wiley series in probability and mathematical statistics.
##' 
##' D. Ginsbourger (2009), \emph{Multiples metamodeles pour l'approximation et
##' l'optimisation de fonctions numeriques multivariables}, Ph.D. thesis, Ecole
##' Nationale Superieure des Mines de Saint-Etienne, 2009.
##' \url{https://tel.archives-ouvertes.fr/tel-00772384}
##' 
##' D. Ginsbourger, R. Le Riche, and L. Carraro (2010), chapter "Kriging is
##' well-suited to parallelize optimization", in \emph{Computational
##' Intelligence in Expensive Optimization Problems}, Studies in Evolutionary
##' Learning and Optimization, Springer.
##' 
##' D.R. Jones (2001), A taxonomy of global optimization methods based on
##' response surfaces, \emph{Journal of Global Optimization}, 21, 345-383.
##' 
##' D.R. Jones, M. Schonlau, and W.J. Welch (1998), Efficient global
##' optimization of expensive black-box functions, \emph{Journal of Global
##' Optimization}, 13, 455-492.
##' 
##' W.R. Jr. Mebane and J.S. Sekhon (2011), Genetic optimization
##' using derivatives: The rgenoud package for R, \emph{Journal of Statistical
##' Software}, \bold{51}(1), 1-55, \url{http://www.jstatsoft.org/v51/i01/}.
##' 
##' J. Mockus (1988), \emph{Bayesian Approach to Global Optimization}. Kluwer
##' academic publishers.
##' 
##' V. Picheny and D. Ginsbourger (2013), Noisy kriging-based optimization
##' methods: A unified implementation within the DiceOptim package,
##' \emph{Computational Statistics & Data Analysis}, 71, 1035-1053. 
##' 
##' C.E. Rasmussen and C.K.I. Williams (2006), \emph{Gaussian Processes for
##' Machine Learning}, the MIT Press, \url{http://www.GaussianProcess.org/gpml}
##' 
##' B.D. Ripley (1987), \emph{Stochastic Simulation}, Wiley.
##' 
##' O. Roustant, D. Ginsbourger and Yves Deville (2012), DiceKriging,
##' DiceOptim: Two R Packages for the Analysis of Computer Experiments by
##' Kriging-Based Metamodeling and Optimization, \emph{Journal of Statistical
##' Software}, \bold{42}(11), 1--26, \url{https://www.jstatsoft.org/article/view/v042i11}.
##' 
##' T.J. Santner, B.J. Williams, and W.J. Notz (2003), \emph{The design and
##' analysis of computer experiments}, Springer.
##' 
##' M. Schonlau (1997), \emph{Computer experiments and global optimization},
##' Ph.D. thesis, University of Waterloo.
##' @examples
##'  
##' set.seed(123)
##' \donttest{
##' ###############################################################
##' ###	2D optimization USING EGO.nsteps and qEGO.nsteps   ########
##' ###############################################################
##' 
##' # a 9-points factorial design, and the corresponding response
##' d <- 2
##' n <- 9
##' design.fact <- expand.grid(seq(0,1,length=3), seq(0,1,length=3)) 
##' names(design.fact)<-c("x1", "x2")
##' design.fact <- data.frame(design.fact) 
##' names(design.fact)<-c("x1", "x2")
##' response.branin <- data.frame(apply(design.fact, 1, branin))
##' names(response.branin) <- "y" 
##' 
##' # model identification
##' fitted.model1 <- km(~1, design=design.fact, response=response.branin, 
##' covtype="gauss", control=list(pop.size=50,trace=FALSE), parinit=c(0.5, 0.5))
##' 
##' ### EGO, 5 steps ##################
##' library(rgenoud)
##' nsteps <- 5
##' lower <- rep(0,d) 
##' upper <- rep(1,d)     
##' oEGO <- EGO.nsteps(model=fitted.model1, fun=branin, nsteps=nsteps, 
##' lower=lower, upper=upper, control=list(pop.size=20, BFGSburnin=2))
##' print(oEGO$par)
##' print(oEGO$value)
##' 
##' # graphics
##' n.grid <- 15
##' x.grid <- y.grid <- seq(0,1,length=n.grid)
##' design.grid <- expand.grid(x.grid, y.grid)
##' response.grid <- apply(design.grid, 1, branin)
##' z.grid <- matrix(response.grid, n.grid, n.grid)
##' contour(x.grid, y.grid, z.grid, 40)
##' title("EGO")
##' points(design.fact[,1], design.fact[,2], pch=17, col="blue")
##' points(oEGO$par, pch=19, col="red")
##' text(oEGO$par[,1], oEGO$par[,2], labels=1:nsteps, pos=3)
##' 
##' ### Parallel EGO, 3 steps with batches of 3 ##############
##' nsteps <- 3
##' lower <- rep(0,d) 
##' upper <- rep(1,d)
##' npoints <- 3 # The batchsize
##' oEGO <- qEGO.nsteps(model = fitted.model1, branin, npoints = npoints, nsteps = nsteps,
##' crit="exact", lower, upper, optimcontrol = NULL)
##' print(oEGO$par)
##' print(oEGO$value)
##' 
##' # graphics
##' contour(x.grid, y.grid, z.grid, 40)
##' title("qEGO")
##' points(design.fact[,1], design.fact[,2], pch=17, col="blue")
##' points(oEGO$par, pch=19, col="red")
##' text(oEGO$par[,1], oEGO$par[,2], labels=c(tcrossprod(rep(1,npoints),1:nsteps)), pos=3)
##' 
##' ##########################################################################
##' ### 2D OPTIMIZATION, NOISY OBJECTIVE                                   ###
##' ##########################################################################
##' 
##' set.seed(10)
##' library(DiceDesign)
##' # Set test problem parameters
##' doe.size <- 9
##' dim <- 2
##' test.function <- get("branin2")
##' lower <- rep(0,1,dim)
##' upper <- rep(1,1,dim)
##' noise.var <- 0.1
##' 
##' # Build noisy simulator
##' funnoise <- function(x)
##' {     f.new <- test.function(x) + sqrt(noise.var)*rnorm(n=1)
##'       return(f.new)}
##' 
##' # Generate DOE and response
##' doe <- as.data.frame(lhsDesign(doe.size, dim)$design)
##' y.tilde <- funnoise(doe)
##' 
##' # Create kriging model
##' model <- km(y~1, design=doe, response=data.frame(y=y.tilde),
##'      covtype="gauss", noise.var=rep(noise.var,1,doe.size), 
##'      lower=rep(.1,dim), upper=rep(1,dim), control=list(trace=FALSE))
##' 
##' # Optimisation with noisy.optimizer
##' optim.param <- list()
##' optim.param$quantile <- .7
##' optim.result <- noisy.optimizer(optim.crit="EQI", optim.param=optim.param, model=model,
##' 		n.ite=5, noise.var=noise.var, funnoise=funnoise, lower=lower, upper=upper,
##' 		NoiseReEstimate=FALSE, CovReEstimate=FALSE)
##' 
##' print(optim.result$best.x)
##' 
##' ##########################################################################
##' ### 2D OPTIMIZATION, 2 INEQUALITY CONSTRAINTS                          ###
##' ##########################################################################
##' set.seed(25468)
##' library(DiceDesign)
##' 
##' fun <- goldsteinprice
##' fun1.cst <- function(x){return(-branin(x) + 25)}
##' fun2.cst <- function(x){return(3/2 - x[1] - 2*x[2] - .5*sin(2*pi*(x[1]^2 - 2*x[2])))}
##' constraint <- function(x){return(c(fun1.cst(x), fun2.cst(x)))}
##' 
##' lower <- rep(0, 2)
##' upper <- rep(1, 2)
##' 
##' ## Optimization using the Expected Feasible Improvement criterion
##' res <- easyEGO.cst(fun=fun, constraint=constraint, n.cst=2, lower=lower, upper=upper, budget=10, 
##'                    control=list(method="EFI", inneroptim="genoud", maxit=20))
##' 
##' cat("best design found:", res$par, "\n")
##' cat("corresponding objective and constraints:", res$value, "\n")
##' 
##' # Objective function in colour, constraint boundaries in red
##' # Initial DoE: white circles, added points: blue crosses, best solution: red cross
##' 
##' n.grid <- 15
##' test.grid <- expand.grid(X1 = seq(0, 1, length.out = n.grid), X2 = seq(0, 1, length.out = n.grid))
##' obj.grid <- apply(test.grid, 1, fun)
##' cst1.grid <- apply(test.grid, 1, fun1.cst)
##' cst2.grid <- apply(test.grid, 1, fun2.cst)
##' filled.contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), nlevels = 50,
##'                matrix(obj.grid, n.grid), main = "Two inequality constraints",
##'                xlab = expression(x[1]), ylab = expression(x[2]), color = terrain.colors, 
##'                plot.axes = {axis(1); axis(2);
##'                             contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), 
##'                                     matrix(cst1.grid, n.grid), level = 0, add=TRUE,
##'                                     drawlabels=FALSE, lwd=1.5, col = "red")
##'                             contour(seq(0, 1, length.out = n.grid), seq(0, 1, length.out = n.grid), 
##'                                     matrix(cst2.grid, n.grid), level = 0, add=TRUE,drawlabels=FALSE,
##'                                     lwd=1.5, col = "red")
##'                             points(res$history$X, col = "blue", pch = 4, lwd = 2)       
##'                             points(res$par[1], res$par[2], col = "red", pch = 4, lwd = 2, cex=2) 
##'                }
##' )}

NULL



