

# Interpolation of Functions
# --------------------------


lin_ipol <- function(gridpts,funcpts,xpts)


# gridpts: N times d matrix       
# funcpts: N times d matrix
# xpts:    n times d matrix 

# (gridpts[,j],funcpts[,j]) is the vector of (x,f_j(x)) values
# xpts[,j] are the points where the function f_j is to be interpolated


{     n <- length(xpts[,1])
      d <- length(xpts[1,])
      
      result <- matrix(0,nrow=n,ncol=d)
           
      for(i in 1:d)
      {   temp       <- approx(gridpts[,i],funcpts[,i],xout=xpts[,i],method="linear",rule=2)
                        # rule = 2 means that the function is extrapolated constantly if necessary
          result[,i] <- temp$y
      }

      return(result)
}




const_ipol_left <- function(gridpts,funcpts,xpts)


# gridpts: N times d matrix       
# funcpts: N times d matrix
# xpts:    n times d matrix 

# (gridpts[,j],funcpts[,j]) is the vector of (x,f_j(x)) values
# xpts[,j] are the points where the function f_j is to be interpolated


{     n <- length(xpts[,1])
      d <- length(xpts[1,])
      
      result <- matrix(0,nrow=n,ncol=d)
           
      for(i in 1:d)
      {   temp       <- approx(gridpts[,i],funcpts[,i],xout=xpts[,i],method="constant",f=0,rule=2)
                        # rule = 2 means that the function is extrapolated constantly if necessary
          result[,i] <- temp$y
      }

      return(result)
}



const_ipol_right <- function(gridpts,funcpts,xpts)


# gridpts: N times d matrix       
# funcpts: N times d matrix
# xpts:    n times d matrix 

# (gridpts[,j],funcpts[,j]) is the vector of (x,f_j(x)) values
# xpts[,j] are the points where the function f_j is to be interpolated


{     n <- length(xpts[,1])
      d <- length(xpts[1,])
      
      result <- matrix(0,nrow=n,ncol=d)
           
      for(i in 1:d)
      {   temp       <- approx(gridpts[,i],funcpts[,i],xout=xpts[,i],method="constant",f=1,rule=2)
                        # rule = 2 means that the function is extrapolated constantly if necessary
          result[,i] <- temp$y
      }

      return(result)
}
