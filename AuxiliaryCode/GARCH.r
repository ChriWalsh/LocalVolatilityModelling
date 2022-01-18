estimation.GARCH <- function(est_results,x_unit,y_unit,data_returns)
{ 
# Inputs:
# ---------
# est_results  -  the estimation results from the smoothbackfitting procedure 
# x_unit:      -  matrix containing the covariate, which were used to estimate 
#                 the additive components via the smooth backfitting algorithm
# y_unit:      -  dependent variable used to estimate the additive model
#                 via the smooth backfitting algorithm
# data_returns -  Dependent variable in the multiplicative model, (for 
#                 observations used to estimate the nonparametric components).
# 
# Outputs:  
# ---------
# A list with 
# coef.Garch.m: The estimated GARCH coefficients when the residual is constructed from interpolating the multiplicative model
# coef.Garch.a: The estimated GARCH coefficients when the residual is constructed from interpolating the additive model
# fit.Garch.m:  The fitted conditional variance for residual construction from multiplicative model
# fit.Garch.a:  The fitted conditional variance for residual construction from additive model

  # Extraction of GARCH errors
  # --------------------------
  
  #1a. Extract GARCH errors by interpolating multiplicative functions (and set variance to one)
  #---------------------------------------------------------------------------------------------------
 
  m.ipol   <- lin_ipol(est_results$grid,exp(est_results$sbf),x_unit)
  
  product  <- matrix(1,ncol=1,nrow=dim(m.ipol)[1])
  for(i in 1:dim(m.ipol)[2])
      product <- product * m.ipol[,i]
  
  eps.m     <- data_returns/(exp(est_results$m0) * product)^(1/2)
  var.eps.m <- mean(eps.m^2)

  #1b. Extract GARCH errors by interpolating additive functions (and set variance to one)
  #--------------------------------------------------------------------------------------
  
  m.ipol <- lin_ipol(est_results$grid,est_results$sbf,x_unit)
  m.ipol <- sqrt(exp(m.ipol))							

  product <- matrix(1,ncol=1,nrow=dim(m.ipol)[1])
  for(i in 1:dim(m.ipol)[2])
      product <- product * m.ipol[,i]
  
  eps.a  <- data_returns / (sqrt(exp(est_results$m0))*product)
  var.eps.full <- mean(eps.a^2)

  #1c. Extract GARCH errors for the Feng model by interpolating multiplicative function 
  #------------------------------------------------------------------------------------  
     
  y.reg <- y_unit
  x.reg <- x_unit[,1]
  grid  <- seq(1/Ngrid,1,by=1/(Ngrid)) 
  # Set the bandwidth equal to the one chosen in time direction of the more general model
  bw.nw <- est_results$bw[1]

  # Epanechnikow Kernel 
  epan  <- function(x){return(0.75*(1-x^2)*((sign(1-x^2)+1)/2))}

  # Nadaraya-Watson estimator
  nw.fun <- function(x.reg,y.reg,bw.nw){
    nom     <- rep(0,Ngrid)
		denom   <- rep(0,Ngrid)
		weight0 <- rep(0,Ngrid)
		weight1 <- rep(0,Ngrid)
		nw      <- rep(0,Ngrid)

		for(j in 1:Ngrid){
		  weight0[j] <- sum(epan((x.reg - grid[j])/bw.nw))
 		  weight1[j] <- sum(epan((x.reg - grid[j])/bw.nw) * y.reg)
		}

		for(j in 1:Ngrid){   
		  if(weight0[j] < 0.00001)
 		    nw[j] <- 0
   	  else
  	    nw[j] <- weight1[j]/weight0[j]
		 }
		return(nw)
	}
  
  trend.only          <- exp(nw.fun(x.reg,y.reg,bw.nw))
  trend.only.values   <- lin_ipol(as.matrix(est_results$grid[,1]),as.matrix(trend.only),as.matrix(x_unit[,1]))
  
  eps.trend.only      <- data_returns/sqrt(trend.only.values)
  var.eps.trend.only  <- mean(eps.trend.only^2)

  #2. GARCH estimation using garchFit from the "fGarch" package.
  
  Fit.Garch.m          <- garchFit(formula = ~ garch(1, 1), data=eps.m,include.mean=FALSE,cond.dist="QMLE")
  Fit.Garch.a          <- garchFit(formula = ~ garch(1, 1), data=eps.a,include.mean=FALSE,cond.dist="QMLE")
  Fit.Garch.trend.only <- garchFit(formula = ~ garch(1, 1), data=eps.trend.only,include.mean=FALSE,cond.dist="QMLE")
  Fit.Garch.simple     <- garchFit(formula = ~ garch(1, 1), data=data_returns,include.mean=FALSE,cond.dist="QMLE")  

  result <- list(coef.Garch.m=Fit.Garch.m@fit$coef,
                 coef.Garch.a=Fit.Garch.a@fit$coef,
                 coef.Garch.trend.only=Fit.Garch.trend.only@fit$coef,
                 coef.Garch.simple=Fit.Garch.simple@fit$coef,
                 fit.Garch.m=Fit.Garch.m@h.t,
                 fit.Garch.a=Fit.Garch.a@h.t,
                 fit.Garch.trend.only=Fit.Garch.trend.only@h.t,
                 fit.Garch.simple=Fit.Garch.simple@h.t,
                 trend.only=trend.only,
                 var.eps.trend.only=var.eps.trend.only,
                 var.eps.full = var.eps.full,
                 eps.m=eps.m,
                 eps.a=eps.a,
                 eps.trend.only=eps.trend.only)

  return(result)

}

