bw_plugin <- function(x,y,bw,data_returns,Max=100){
    # x        regressors 
    # y        response    
    # bw       bandwidth vector (starting values)
    # Max      maximal number of loops for bandwidth selection
    
    # Additional tuning parameter used to estimate derivative of density as in Park and Mammen (2005)
    c.deriv <- 2
  
    # Include possibility to keep track of bw iterates
    bw.iterates     <- matrix(0,ncol=dim(x)[2],nrow=(Max+1))
    bw.iterates[1,] <- bw
    
    # Store the different estimates for the bias and variance
    bias.iterates   <- matrix(0,ncol=dim(x)[2],nrow=(Max+1))
    
    var.hom.iterates <- matrix(0,ncol=dim(x)[2],nrow=(Max+1))
    
    
    library(stats)
    pos <- 1
    while(pos <= dim(x)[1])
    {     if(prod(x[pos,] <= 1 & x[pos,] >= 0) == 0)
          {  x <- x[-pos,]
             y <- y[-pos]         
          }   
          else
             pos <- pos+1
    }
    
    epan <- function(x){return(0.75*(1-x^2)*((sign(1-x^2)+1)/2))}
    epan.deriv <- function(x){return(-1.5*x*((sign(1-x^2)+1)/2))}
  

    bw_minus_2 <- 0
    bw_old <- 0.5 * bw
    loops  <- 0 
    

    
    # Loop for the bandwidth Selection
    #----------------------------------
    # (Less stringent stopping criterion than in Park and Mammen, which had relative change of 10^{-3}. 
    #  With addition that the change between the current and 2 iterations ago should be small.)     
    
    while(prod(abs((bw - bw_old)/bw_old) < 0.01) == 0 & prod(abs((bw - bw_minus_2)/bw_minus_2) < 0.0001) == 0  & loops < Max)
    {
         # Print the iteration of the bandwidth selection procedure
         loop.print <- paste("Loop", loops ,"out of max of", Max ,"in bandwidth selection procedure")
         print(loop.print)   

         # Run the smooth backfitting algorithm  
         add <- sbfnw(x,y,bw,N=Ngrid)
    
         # Extract the additive fit and interpolate to get values at data points
         m_add   <- add$sbf
         gridpts <- add$grid
         sbfvalues <- lin_ipol(gridpts,m_add,x)
    
         # Get residuals from the additive fit
         sbfsums   <- rowSums(sbfvalues)
         resid     <- (y - mean(y) - sbfsums)
   

         # Step 1: Calculate an estimator for the asymptotic variance
         #-------------------------------------------------------------
         #
         # (i) Asymptotic variance for the time component is based on long run variance formula. 
         # This uses asymptotic variance expression for log transform of the sample mean of 
         # squared GARCH(1,1) and the ARMA(1,1) representation of a squared GARCH process. For details see 
         # Section S.4 of supplementary material of the paper.  
         # (ii) Asymptotic variance for stochastic components is estimated by the empirical variance of the additive residuals.
         # (This simplifies from having to estimate the conditional variance for each direction.)
    
         
         # Initialize the asymptotic variance
         asy.variance  <- vector(mode="double",length=dim(x)[2])
      
         # (i) Fit GARCH model and use formula for long run variance
    
         # Extract GARCH errors by interpolating multiplicative functions
         m.ipol   <- lin_ipol(add$grid,exp(add$sbf),x_unit)
    
         product  <- matrix(1,ncol=1,nrow=dim(m.ipol)[1])
         for(i in 1:dim(m.ipol)[2])
            product <- product * m.ipol[,i]
         eps.m     <- data_returns/(exp(add$m0) * product)^(1/2)
    
         # Save the squared residual as a time series object
         eps.m.sq.ts  <- as.ts(eps.m*eps.m)
    
         # Fit an ARMA(1,1) model to the squared residual process.
         para.fit <- arima(eps.m.sq.ts,order=c(1,0,1),include.mean=TRUE)
    
         # Extract the parameter estimates
         coeff.sq.ARMA <- coef(para.fit)
         var.eta.ARMA <- para.fit$sigma2
    
         # Calculate the long-run variance based on asymptotic variance formula 
         # for mean of log of ARMA(1,1) obtained from delta method.
         lr.var <- var.eta.ARMA/coeff.sq.ARMA[3]^2*(1+coeff.sq.ARMA[2])^2
         asy.variance[1] <- lr.var
         
         # (ii) Estimate the conditional variance expressions by approximating it with the 
         # empirical variance of the residuals, i.e. using the homoscedastic case.
             
         var.hom <- var(resid) 
         for(k in 2:dim(x)[2])
            asy.variance[k] <- var.hom

         # Rescale by the factor from the Epanechnikov Kernel 
         asy.variance <- 3/5 * asy.variance
         
         # Step 2: Calculate an estimator for the asymptotic bias
         # (Plugin estimateds for terms in leading term of asymptotic bias 
         #  under independence of regressors)
         #-----------------------------------------------------------------------------------
         
         # Initialize bias 
         bias <-  matrix(0,nrow=dim(x)[1],ncol=dim(x)[2])    
         
         
         # Initialize density and derivative of density 
         dens        <- matrix(0,nrow=dim(x)[1],ncol=dim(x)[2]) 
         dens.deriv  <- matrix(0,nrow=dim(x)[1],ncol=dim(x)[2]) 
    
         # Get density estimates at observed covariate values
         # Loop over dimensions/columns of x
         for(k in 1:dim(x)[2])
            {   # Loop over observations/rows of x
                for(j in 1:dim(x)[1])
                {   weight     <- epan((x[,k] - x[j,k]) / (bw[k]))
                    dens[j,k]  <- 1/(dim(x)[1]*bw[k])*sum(weight)    
                }       
            }
    
         # Get estimated density derivative at observed covariate values
         # Loop over dimensions/columns of x
         for(k in 1:dim(x)[2])
         {   # Loop over observations/rows of x
             for(j in 1:dim(x)[1])
             {   weight     <- epan.deriv((x[,k] - x[j,k]) / (c.deriv*bw[k]))
                 dens.deriv[j,k]  <- 1/(dim(x)[1]*(c.deriv*bw[k])^2)*sum(weight)    
             }       
         }
    
         # Estimate first and second derivatives of univariate regression functions at observed covariate values
         # (Using local quadratic approximation to estimated univariate additive components.)
         m_1st.deriv <- matrix(0,nrow=dim(x)[1],ncol=dim(x)[2])  
         m_2nd.deriv <- matrix(0,nrow=dim(x)[1],ncol=dim(x)[2])
    
         # Loop over dimensions/columns of x
         for(k in 1:dim(x)[2])
         {   # Loop over observations/rows of x
             for(j in 1:dim(x)[1])
             {   reg_y    <- m_add[,k]
                 reg_x1   <- gridpts[,k] - x[j,k]
                 reg_x2   <- (gridpts[,k] - x[j,k]) * (gridpts[,k] - x[j,k])
                 weight   <- epan((gridpts[,k] - x[j,k]) / (c.deriv * bw[k]))
                 quad_fit <- lm(reg_y ~ reg_x1 + reg_x2, weights=weight)
                 m_1st.deriv[j,k] <- as.double(quad_fit$coefficients[2])
                 m_2nd.deriv[j,k] <- 2 * as.double(quad_fit$coefficients[3])       
              }       
         }
    
         # For time direction use simplified formula
         m_1st.deriv[,1] <- 0
         
         # Construct inverse of density at observed covariate values needed in bias formula
         dens.inv <- matrix(0,nrow=dim(x)[1],ncol=dim(x)[2]) 
    
         # Loop over dimensions/columns of x
         for(k in 1:dim(x)[2])
         {   # Loop over observations/rows of x
             for(j in 1:dim(x)[1])
             {   if(dens[j,k]<0.00001)
                 {    dens.inv[j,k] <- 0 
                 } else {
                 dens.inv[j,k]  <- 1/dens[j,k]
                 } 
             }
         }
    
         # Formula for leading term of bias from Mammen and Park 
         for(k in 1:dim(x)[2])
             bias[,k] <- 1/5*(m_1st.deriv[,k]*dens.deriv[,k]*dens.inv[,k]+1/2*m_2nd.deriv[,k])         
    
         # Normalize the bias components.   
         bias.norm <- rep(0,dim(x)[2])
    
         for(k in 1:dim(x)[2])
         {   bias.norm[k] <- mean(bias[,k])
         } 
    
         bias <- bias - bias.norm
    
         # Calculate the estimate of the squared asymptotic bias
         bias <- 1/dim(x)[1]*colSums(bias * bias)
    
         # --- Update bandwidth --- #
         # Update bandwith values from the last two iteration loops
         bw_minus_2 <- bw_old
         bw_old <- bw
    
         # Plug-in the formula for the estimate asymptotic bias and variance into the leading term
         # of the formula for the optimal bandwidth choice.
         bw     <- dim(x)[1]^(-1/5) * (asy.variance/4*bias)^(1/5)
    
         # Increase the loop counter and add the estimated bandwidth
         loops  <- loops + 1
         bw.iterates[loops+1,] <- bw    
    }

    # Return the selected bandwidth, the number of loops in the selection procedure and the
    # progressions of the selected bandwidth in the procedure.
    return(list(bw,loops,bw.iterates))
}
