###########################
#                         #
#  Transformation of Data #
#                         #
###########################




lags <- function(x,lag.length,dates)

{    if(dates==TRUE)
        x <- x[,-1] 
        x <- as.matrix(x)

     cols <- dim(x)[2]
     rows <- dim(x)[1]
 
     result <- rbind(matrix(NA,ncol=cols,nrow=lag.length),as.matrix(x[1:(rows-lag.length),]))

     return(result)
}

leads <- function(x,lead.length,dates)

{    if(dates==TRUE)
        x <- x[,-1] 
        x <- as.matrix(x)
        
     cols <- dim(x)[2]
     rows <- dim(x)[1]
    
     result <- rbind(as.matrix(x[(lead.length+1):rows,]),matrix(NA,ncol=cols,nrow=lead.length))

     return(result)
}


dif <- function(x,d,dates)

{    if(dates==TRUE)
        x <- x[,-1] 

     x <- as.matrix(x)
         
     cols <- dim(x)[2]
     rows <- dim(x)[1]

     result <- rbind(matrix(NA,ncol=cols,nrow=d),diff(x,lag=d))

     return(result)
}     

growth_rates <- function(x,d,dates)

# Input: Matrix with time series in its columns. (No rescaled time in first column!)

{   if(dates==TRUE)
       x <- x[,-1]  

    x <- as.matrix(x)
    
    cols <- dim(x)[2]                             
    rows <- dim(x)[1] 
    
    for(i in 1:cols)
    {   x1          <- as.vector(x[,i])
        x1          <- (x1[d:rows] - x1[1:(rows-d+1)]) / x1[1:(rows-d+1)] 
        D           <- d-1                                       
        x[-(1:D),i] <- x1
        x[1:D,i]    <- NA
    }
               
    return(x)
}

av.diff <- function(x,d,dates)

{   if(dates==TRUE)
       x <- x[,-1]  

    x <- as.matrix(x)

    cols <- dim(x)[2]                             
    rows <- dim(x)[1]

    result <- matrix(NA,ncol=cols,nrow=rows)
    
    for(j in 1:cols)
    {   for(i in (2*d):rows)  
        {   temp <- mean(x[(i-(d-1)):i,j]) - mean(x[((i-d)-(d-1)):(i-d),j])
            #temp <- x[i,j] - mean(x[((i-d)-(d-1)):(i-d),j])  
            result[i,j] <- temp
        }
    }

    return(result)
}    


scaling <- function(x,bound)

{   x <- as.matrix(x)
    n <- dim(x)[1]
    d <- dim(x)[2]
             
    x_sorted <- matrix(0,ncol=d,nrow=n)
    distance <- matrix(0,ncol=d,nrow=n)                       
             
    LOWER <- rep(2,d)
    UPPER <- rep(n,d)
             
    for(i in 1:d)   
    {   x_sorted[,i] <- sort(x[,i])
              
        for(j in 2:n)
            distance[j,i] <- x_sorted[j,i] - x_sorted[(j-1),i]
               
        lower  <- 2
        lower0 <- 2
        upper  <- n
        upper0 <- n
        crit   <- 1
                                                                             
        while(crit > bound)
                 
        {     lower  <- lower - 1 + which.max(distance[lower:floor(n/2),i]) 
              upper  <- ceiling(n/2) - 1 + which.max(distance[ceiling(n/2):upper,i]) 
      
              if(lower > upper)
                 stop("Criterion too hard!") 
                 
              if((x_sorted[upper,i] - x_sorted[lower0-1,i]) == 0 | (x_sorted[upper0,i] - x_sorted[lower-1,i]) == 0)
                 stop("Attention: Division by zero!")
      
              if(distance[lower,i] < distance[upper,i])
              {  dist.max <- distance[upper,i]
                 crit     <- dist.max / (x_sorted[upper0,i] - x_sorted[lower0-1,i])                          
                 if(crit > bound)
                 {  upper    <- upper - 1
                    UPPER[i] <- upper   
                    upper0   <- upper
                 } 
              }
              else
              {  dist.max <- distance[lower,i] 
                 crit     <- dist.max / (x_sorted[upper0,i] - x_sorted[lower0-1,i])                          
                 if(crit > bound)
                 {  lower    <- lower + 1
                    LOWER[i] <- lower
                    lower0   <- lower
                 }
              }   
        }                 

        x[,i] <- x[,i] - x_sorted[LOWER[i]-1,i] 
        x[,i] <- x[,i] / (x_sorted[UPPER[i],i] - x_sorted[LOWER[i]-1,i]) 

        LOWER[i] <- x_sorted[LOWER[i]-1,i]
        UPPER[i] <- x_sorted[UPPER[i],i]
    }

    return(list(x,LOWER,UPPER))
}



####################
#                  #
# Cleaning of Data #
#                  #
####################



complete_date <- function(x,dates)

# Adds additional row for each missing date and
# gives the series the value NA at these dates 
# Input: matrix x with two columns, 
# the first of which contains dates.

                 {   x      <- as.matrix(x)                
		     dates1 <- match(dates,x[,1])
                     X      <- matrix(NA,nrow=length(dates),ncol=dim(x)[2])
                     X[,1]  <- dates             
                     pos    <- 1                     

                     for(i in 1:length(dates))
                     {   if(!is.na(dates1[i]))
                         {  X[i,-1] <- x[pos,-1]
                            pos <- pos+1
                         }
                     }
                
                     return(X) 
                 }


   
delete_na <- function(x)

# Deletes observations where one variable
# has a value "NA", "Inf" or "-Inf" for the return data

             {   m  <- 0
                    
                 for(j in 1:length(x[1,])) 
                 {   pos <- 1
                     for(i in 1:length(x[,1]))             
                     {   if(is.na(x[pos,j]) | x[pos,j] == Inf | x[pos,j] == -Inf)                    
                            x <- x[-pos,]
                         else
                            pos <- pos + 1   
                     }
                 } 
                  
                 return(x)
             }



###################################
#                                 #
# Generating rescaled time vector #
#                                 #
###################################



x_withx0 <- function(x)

# Generates rescaled time vector and sticks
# it as first column into the x-matrix

	    {   x     <- as.matrix(x)
                cols  <- dim(x)[2]                             
                rows  <- dim(x)[1] 
           
	        x0    <- as.matrix(seq(1/rows,1,by=1/rows)) 
                x     <- matrix(c(x0,x),ncol = cols + 1, byrow=FALSE)
             
                return(x)
	    }                          
		



