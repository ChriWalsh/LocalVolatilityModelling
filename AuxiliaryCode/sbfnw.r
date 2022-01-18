sbfnw <- function(x,y,bw,grids=FALSE,N,iterate=50)
{

# x  n*d matrix of predictors
# y  vector of regressors
# bw vector of bandwidths
# N  number of gridpoints =  ng in ll


      n <- as.integer(dim(x)[1])
      d <- as.integer(dim(x)[2])
      
      if (!n == length(y))
             stop("Predictor and Response of different length")

      if (!length(bw) == d)
             stop("Wrong dimension of bandwidth vector")
      
      storage.mode(x)       <- "double"
      storage.mode(y)       <- "double"
      storage.mode(bw)      <- "double"
      storage.mode(N)       <- "integer"
      storage.mode(iterate) <- "integer"
  
      if(grids==TRUE)
         Grid <- Grid1
      else
         Grid <- rep(seq(0,1,1/(N-1)),d)
      
      L <- ceiling(0.5*d*(d-1)*N)
      L <- as.integer(L)

            
      xvec     <- as.vector(x)
      gridvec  <- as.vector(Grid)
      m0       <- vector(mode="double",length=1)
      kde      <- vector(mode="double",length=d*N)
      regbf    <- vector(mode="double",length=d*N)
      regnw    <- regbf       
      conv     <- vector(mode="double",length=d)
      fh_int   <- vector(mode="double",length=d)
      fh_int2a <- vector(mode="double",length=L)
      fh_int2b <- vector(mode="double",length=L)  
      reg_int  <- vector(mode="double",length=d) 
      

      result <- .C("sbfnw", xvec, y, n, d, gridvec, N, bw, kde, m0, regbf, regnw,
                   iterate, conv, fh_int, fh_int2a, fh_int2b, reg_int)
          
      list(grid=matrix(result[[5]], ncol=d),bw=result[[7]], kde=matrix(result[[8]], ncol=d),
           m0=result[[9]], sbf=matrix(result[[10]], ncol=d), nwr=matrix(result[[11]], ncol=d),
           iteration=result[[12]], conv=result[[13]], fh_int=matrix(result[[14]],ncol=d),
           fh_int2a=matrix(result[[15]],ncol=ceiling(0.5*d*(d-1))),
           fh_int2b=matrix(result[[16]],ncol=ceiling(0.5*d*(d-1))), 
           nwr_int=matrix(result[[17]],ncol=d)) 
}
     




