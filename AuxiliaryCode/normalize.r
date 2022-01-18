normalize.tau.functions <- function(est_results,x_unit,y_unit,x.normalize.point,conf,data_returns){
  # Inputs:
  # ---------
  # est_results       -  the estimation results from the smoothbackfitting procedure 
  # x_unit:           -  matrix containing the covariate, which were used to estimate 
  #                      the additive components via the smooth backfitting algorithm
  # y_unit:           -  dependent variable used to estimate the additive model
  #                      via the smooth backfitting algorithm
  # x.normalize.point -  the normalization point for the plots of the estimates
  # conf:             -  (pointwise) confidence bands 
  # data_returns:     -  original dependent variable, i.e. the (log) returns 
  
  
  # Outputs:  
  # ---------
  # List containing:
  # tau.sqrd:       - Normalized estimated of multiplicative components.
  # tau.sqrd.upper: - Upper pointwise confidence band (95%-CI).
  # tau.sqrd.lower: - Corresponding lower pointwise confidence band.
  # tau0.sqrd :     - Unconditional volatility level 
  
  
   #---------------------------------------------------------------------------
   # (1) Estimate the GARCH errors by interpolating additive functions
   #    (and normalizing variance to one)
   #---------------------------------------------------------------------------
  
   m.ipol          <- lin_ipol(est_results$grid,est_results$sbf,x_unit)
   tau.sqrd.ipol.a <- exp(m.ipol)
   m.ipol          <- sqrt(exp(m.ipol))		   
   # m.ipol gives value of tau functions at data points				
      
   product <- matrix(1,ncol=1,nrow=dim(m.ipol)[1])
   for(i in 1:dim(m.ipol)[2])
       product <- product * m.ipol[,i]
  
   eps.a       <- data_returns / (sqrt(exp(est_results$m0))*product)
   var.eps.a0  <- mean(eps.a^2)
   var.eps.a   <- mean((eps.a-mean(eps.a))^2)     
   # Variance of GARCH error imposing the theoretical mean of 0

   eps.a  <- eps.a / sqrt(var.eps.a)
   # Normalizng the GARCH errors to have variance of 1

   #---------------------------------------------------------------------------
   # (2) Normalization of tau functions given the normalization point 
   #---------------------------------------------------------------------------
   tau.sqrd     <- exp(est_results$sbf)
   variance.eps <- var.eps.a 	
   eps.sqrd     <- eps.a^2

   conf.upper <- conf$conf.mult.upper
   conf.lower <- conf$conf.mult.lower  

   cl      <- dim(tau.sqrd)[2] - 1
   m.norm  <- lin_ipol(as.matrix(est_results$grid[,-1],ncol=cl),as.matrix(tau.sqrd[,-1],ncol=cl),x.normalize.point)

   for(k in 2:dim(est_results$sbf)[2])
     { tau.sqrd[,k]   <- tau.sqrd[,k] / m.norm[k-1]
       conf.upper[,k] <- conf$conf.mult.upper[,k]/m.norm[k-1]
       conf.lower[,k] <- conf$conf.mult.lower[,k]/m.norm[k-1]
     }	
   tau.sqrd[,1]   <- tau.sqrd[,1] * prod(m.norm) * variance.eps
   conf.upper[,1] <- conf$conf.mult.upper[,1]*prod(m.norm)*variance.eps
   conf.lower[,1] <- conf$conf.mult.lower[,1]*prod(m.norm)*variance.eps
   
   # --- Bands for (scaled) multiplicative model ---    

   conf.upper <- (exp(conf$conf.add.upper-est_results$sbf) - 1) * tau.sqrd
   conf.lower <- (exp(conf$conf.add.lower-est_results$sbf) - 1) * tau.sqrd

   tau.sqrd.upper <- tau.sqrd + conf.upper
   tau.sqrd.lower <- tau.sqrd + conf.lower

   #---------------------------------------------------------------------------
   # (3) Normalization such that trend component can be interpreted as 
   #     unconditional volatility level. 
   #---------------------------------------------------------------------------

   eps.tau0.sqrd   <- data_returns^2 / tau.sqrd.ipol.a[,1]

   variance.eps.tau0<- mean(eps.tau0.sqrd)
   tau0.sqrd        <- variance.eps.tau0*exp(est_results$sbf[,1])  
   

   result <- list(tau.sqrd=tau.sqrd, tau.sqrd.upper=tau.sqrd.upper, tau.sqrd.lower=tau.sqrd.lower, tau0.sqrd = tau0.sqrd)

   return(result)

}






