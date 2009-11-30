qEI <- function 
(newdata, model, type="UK", MC.samples=10000, return.I=FALSE) 
{
   library(MASS) #Useful for simulation of Gaussian Vectors

   p <- predict.km(object=model, newdata=newdata, type=type,   	cov.compute=TRUE) #Prediction, along with cond cov

   # Conditional simulation of the process at newdata	
   cond.simu <- mvrnorm(n=MC.samples, mu=p$mean, Sigma=p$cov) 
	
   # Matrix with improvement values for each point and realization
   improvements.matrix <- 0.5*( (min(model@y)-cond.simu) + abs(min      
   (model@y)-cond.simu) ) #positive part
   
   # Computation of multipoint improvement values at newdata
   qI <- rep(0, MC.samples) 
   for(i in seq(1,MC.samples)){qI[i] <- max(improvements.matrix[i,])}

   # Computation of the multipoint improvement estimate, and its sd
   qEI.estimate <- mean(qI)
   qEI.sd.estimate <- sd(qI)/sqrt(MC.samples)
	
	# Case where improvements.matrix is needed amon the outputs
	if(return.I==TRUE){
   	return(
	res <- list(qEI=qEI.estimate, qEI.sd=qEI.sd.estimate,  	I=improvements.matrix))}
	
	# Case where improvements.matrix is not  needed
	else{
   	return(
	res <- list(qEI=qEI.estimate, qEI.sd=qEI.sd.estimate, I=NULL))}
}

