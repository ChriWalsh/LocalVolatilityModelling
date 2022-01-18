scale.data <- function(YX,Bound){

# Inputs: 
# --------------
# YX:    matrix with 1st column containing id/date, 2nd column containing
#        dependent variable and further columns containing covariates
#        (in particular 4th column, which is the first covariate
#         should contain rescaled time.) 
# Bound: vector of length equal to dimension of covariates. Indicates
#        maximum distance between points in each dimension once
#        rescaled to unit interval  
#
# Outputs: 
# --------------
# List containing:
# x:           Covariates rescaled so that 0 and 1 are the cutoff points
# y_add:       Dependent variable squared and logged (Nothing is done to it here)
# x_unit:      Covariates rescaled to [0,1]^d with observations outside thrown away
# y_unit:      Values of squared and logged dependent variables for which rescaled covariates are not thrown away
# x_upper:     Value of covariate that is rescaled to 1
# x_lower:     Value of covariate that is rescaled to 0
# x.original:  Actual values of covariates (Nothing is done to it here)
# x.estimation:Actual covariates with values that would be outside [0,1]^d after rescaling replaced by NA
# id:          Identification number of data row for the variables in [0,1]^d
# y_mult_sqrt: Actual dependent variables for which covariates are in [0,1]^d


   id          <- YX[,1]			 
   y_mult_sqrt <- YX[,2]
   y_add       <- YX[,3]
   x           <- YX[,-(1:3)]										
   
   if(dim(x)[2] != length(Bound))
     stop('Dimension of Bound is incorrect!') 		
    
   # Scale covariate series to [0,1]^d

   x_lower      <- rep(0,dim(x)[2]) 
   x_upper      <- rep(1,dim(x)[2])

   x_normal     <- x
   x_scaled     <- x
   x.estimation <- x
   x.original   <- x
   id.original  <- id   

   y_unit       <- y_add

   empty_entry        <- rep(NA,dim(x.original)[2])

   eta           <- 1           # Initialize while loop
   deleted.count <- 0           # Count number of deleted observations

   # Iterate scaling over all dimensions of x

   while(eta > 0)
    {   dim.old <- dim(x_normal)[1]
        
        for(i in 1:dim(x_normal)[2])			
	{   temp          <- scaling(x_normal[,i],Bound[i])

	    # Scaling takes a vector and bound and returns another vector with the elements scaled to [0,1] such that 
            # no rescaled 'adjacent' points are further than bound away from each other.
            # x_lower and x_upper are the original value of the points now on 0 and 1. 		

    	    x_scaled[,i]  <- temp[[1]]
            x_lower[i]    <- temp[[2]]
    	    x_upper[i] 	  <- temp[[3]] 
	} 

	# Drop points whose covariates after rescaling are not in [0,1]^d 	

	pos              <- 1
        tot              <- 1        

	while(pos <= dim(x_scaled)[1])
	{   if(prod(x_scaled[pos,] <= 1 & x_scaled[pos,] >= 0) == 0)
   	 {  x_scaled           <- x_scaled[-pos,]
       	    x_normal           <- x_normal[-pos,]               # Not sure whether one needs x_normal
            y_mult_sqrt        <- y_mult_sqrt[-pos] 
            y_unit             <- y_unit[-pos]
            id                 <- id[-pos]
            x.estimation[tot,] <- empty_entry
            
            
            tot                <- tot + 1
            deleted.count      <- deleted.count + 1
        }   
  	    else
     	{    tot <- tot + 1
             pos <- pos + 1
	}
        }
	eta      <- dim.old - dim(x_normal)[1]
   }

   x_unit <- x_scaled

   for(i in 1:dim(x)[2])                                               
   {  x[,i] <- (x[,i] - x_lower[i])/(x_upper[i]-x_lower[i])
   }

   result <- list(id=id,x=x, y=y_add, x_unit=x_unit, y_unit=y_unit,x_upper=x_upper,x_lower=x_lower,x.original=x.original,x.estimation=x.estimation,id=id,y.mult.sqrt=y_mult_sqrt) 
   return(result)

}

