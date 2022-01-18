data.additive <- function(YX){

#   Inputs:
# --------------------------------------------
#    YX    - matrix with 1st column containing dependent variable
#            and other columns containing covariates
# --------------------------------------------

#   Outputs:
# --------------------------------------------
# A list containing:
# id: 	       Identification number of the "data row"
# y_mult_sqrt: Raw dependent variable 
# y_add:       Dependent variable squared and logged, i.e. y_add = (log(y_mult_sqrt))^(1/2)
# x:           List of potential covariates
# In particular all data rows with a missing variable have been deleted.   


#    Note: Important that order in YX is correct


   #Remove observations where one variable has "NA" 
   #------------------------------------------------
   
   YX <- delete_na(YX)     
  
   # Replace the date column with the row number
  # YX[,1] <- rownames(YX)
  # colnames(YX)[1] <- "id" 
   
   # Extract dependent variable and covariates
   id <- YX[,1]
   y_mult_sqrt <- YX[,2] 
   y_mult      <- y_mult_sqrt^2
   y_add 	   <- log(y_mult)
   x 	   <- as.matrix(YX[,-(1:2)],nrow=dim(YX)[1])
   colnames(x) <- colnames(YX)[-(1:2)]
   
   # Get rid of possible problems when defining y for zero returns
   if(min(y_add)== -Inf){ 
     YX <- cbind(id,y_mult_sqrt,y_mult,y_add,x)
     YX <- delete_na(YX)
   
     id 	   <- YX[,1]
     y_mult_sqrt   <- YX[,2] 
     y_mult        <- YX[,3]
     y_add         <- YX[,4]
     x 	           <- as.matrix(YX[,-(1:4)],nrow=dim(YX)[1])
     colnames(x)   <- colnames(YX)[-(1:4)]
   }   
   
   # Add a rescaled time trend variable 
   x.name <- colnames(x)
   x <- cbind(seq(0,1,length.out=dim(x)[1]),x)						
   colnames(x) <- c("Rescaled time",x.name)
   
   result <- cbind(id,y_mult_sqrt,y_add,x)
   result <- matrix(as.numeric(result),ncol=ncol(result))

   return(result)
}

