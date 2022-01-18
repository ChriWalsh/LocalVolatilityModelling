confidence.bands <- function(est_results,x_unit,y_unit){  

   # Inputs:
   # ---------
   # est_results  -  the estimation results from the smoothbackfitting procedure 
   # x_unit:      -  matrix containing the covariate, which were used to estimate 
   #                 the additive components via the smooth backfitting algorithm
   # y_unit:      -  dependent variable used to estimate the additive model
   #                 via the smooth backfitting algorithm
   # Outputs:  
   # ---------
   # List containing:
   # conf.add.upper: - Upper pointwise confidence band (95%-CI) for 
   #                   components in additive model.
   # conf.add.lower: - Corresponding lower pointwise confidence band.
   # conf.mult.upper: - Upper pointwise confidence band (95%-CI) for 
   #                    components in multiplicative model.
   # conf.mult.lower: - Corresponding lower pointwise confidence band.
  
  
  
   # Calculation of the Conditional Variance of the Error Term
   # ---------------------------------------------------------

   sbfvalues <- lin_ipol(est_results$grid,est_results$sbf,x_unit)        
   sbfsums   <- rowSums(sbfvalues)
   resid     <- (y_unit - est_results$m0 - sbfsums)
   epan      <- function(x){0.75*(1-x^2)*((sign(1-x^2)+1)/2)}
   var_resid <- matrix(0,ncol = dim(x_unit)[2], nrow = dim(est_results$grid)[1]) 
   reg_y     <- resid^2
   n0        <- dim(x_unit)[1]

   # Estimate the residual variance in time direction
   len <- length(acf(resid,type="correlation",plot=FALSE)$acf)
   
   Bartlett.weight <- rep(1,len)
  
   for( i in 2:len) 
   Bartlett.weight[i] <- 2*(Bartlett.weight[i] - i/(len+1))

   LR.var.resid <- sum(acf(resid,type="covariance",plot=FALSE)$acf*Bartlett.weight)

   var_resid[,1] <- (3/5)*n0^(-1/5)*est_results$bw[1]^(-1)*LR.var.resid

   # Estimate the residual variance in other directions
   for(j in 1:dim(est_results$grid)[1])
   {   for(k in 2:dim(x_unit)[2])
       {    reg_x1    <- x_unit[,k] - est_results$grid[j,k]
            weight    <- epan(reg_x1/est_results$bw[k])
            const_fit <- lm(reg_y ~ 1, weights=weight)
            var_resid[j,k] <- as.double(const_fit$coefficient[1])
            var_resid[j,k] <- (3/5) * n0^(-1/5) * est_results$bw[k]^(-1) * (var_resid[j,k] / est_results$kde[j,k])   
       }
   }

   # Calculation of Confidence Bands
   # -------------------------------


   # --- Bands for additive model --- #

   conf_upper_add <- n0^(-2/5) * qnorm(.975) * sqrt(var_resid)
   conf_lower_add <- n0^(-2/5) * qnorm(.025) * sqrt(var_resid)

   m_upper_add <- est_results$sbf + conf_upper_add
   m_lower_add <- est_results$sbf + conf_lower_add

   # --- Band for multiplicative model --- #
   
   conf_upper_mult <- n0^(-2/5) * qnorm(.975) * sqrt(var_resid) * exp(est_results$sbf)
   conf_lower_mult <- n0^(-2/5) * qnorm(.025) * sqrt(var_resid) * exp(est_results$sbf)

   result <- list(conf.add.upper=m_upper_add,conf.add.lower=m_lower_add,conf.mult.upper=conf_upper_mult,conf.mult.lower=conf_lower_mult)

   return(result)
}


