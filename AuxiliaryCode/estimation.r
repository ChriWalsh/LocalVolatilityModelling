
estimation <- function(x,y,bw.init,bw.method,Ngrid,data_returns){   

  #   Inputs:
  # --------------------------------------------
  #    x            - matrix with covariates that have been scaled to [0,1], 
  #                   the first covariate is rescaled time.  
  #    y            - dependent variable for the additive model. 
  #    bw.init      - (initial) bandwidth vector    
  #    bw.method    - Choose to: (i) run the estimation with bw.init ("adhoc") or
  #                   (ii) use bw.init to start the iterative procedure ("plugin")
  #    Ngrid        - number of equidistant grid points at which to estimate.
  #    data_returns - squared dependent variable of the multiplicative model, 
  #                   (i.e. squared returns.)  
  #
  # --------------------------------------------
  
  #   Outputs:
  # --------------------------------------------
  # A list containing:
  #
  # sbf        - matrix containing the estimated nonparametric function at the 
  #              the equidistant grid
  # grid       - the estimation points, i.e. the equidistant grid
  # m0         - the estimated constant in the additive model
  # bw         - bandwidth vector used to get final estimates.
  #
  # 
  # bw.iterates - matrix illustrating how iterative bandwidth selection proceeded
  # loops       - the number of loops the iterative bandwidth selection procedure
  #               ran
  #
  # kde       - (scaled) univariate kernel density estimates of the regressors
  #             at the grid points.
  # nwr       - (scaled) univariate Nadaraya-Watson (pilot) estimates 
  # iteration - number of iterations needed in the smooth backfitting algorithm  
  # conv      - "distance" between final smooth backfitting estimates and 
  #             last iteration run. 
  # fh_int, fh_int2a, fh_int2b - computed integrals in the smooth backfitting procedure.
  
  # Load the shared library and source the file to call the C code for the 
  # smooth backfitting
  dyn.load("AuxiliaryCode/sbfnw.so","sbfnw")
  source("AuxiliaryCode/sbfnw.r")
  
  # Code to setup bandwidth with user defined bandwidth choice   
  if(bw.method == "adhoc")
  {   bw  <- bw.init  
      loops <- 1
      bw.iterates <- list(bw.iterates=1)
  } 
  #  Code for bandwidth selection as proposed in the supplement of Walsh & Vogt (2022+)
  #  "Locally Stationary Multiplicative Voaltility" (JBES)
  if(bw.method == "plugin")
  {   # Load the bandwidth selection code
      source("AuxiliaryCode/bw_plugin_Loc_Vol.r")
      # Run the bandwidth selection code
      bw_temp  <- bw_plugin(x,y,bw.init,data_returns)
      # Store the bandwidth obtained from the bandwidth selection algorithm
      bw       <- bw_temp[[1]]
      # Store the number of iterations used in the bandwidth selection algorithm 
      loops    <- list(loops=bw_temp[[2]])
      # Store the progression of the selected bandwidth over the loops.
      bw.iterates <- list(bw.iterates=bw_temp[[3]])
  }

  # Run the actual smooth backfitting algorithm with the corresponding chosen or selected bandwidth
  add <- sbfnw(x,y,bw,N=Ngrid) 					
  # Add the progression of the bandwidth selection to the estimation results
  add <- c(add,loops,recursive=FALSE,bw.iterates)

  return(add)
}




