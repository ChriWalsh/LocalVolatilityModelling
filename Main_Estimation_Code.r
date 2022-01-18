rm(list=ls())

#-------------------------------
#
# Step 0: Call all necessary functions.
#
#-------------------------------

library(fGarch)

source("AuxiliaryCode/aux_functions1.r")
# Contains helper functions for interpolating functions

source("AuxiliaryCode/aux_functions2.r")
# Contains helper functions for preparing data



#--------------------------------------------------
#
# Step 1: Read in the data used for the application  
#
#--------------------------------------------------

mult_data <- read.csv("Data/Data_Application.csv")

# The model in the application uses 3 regressors obtained from FRED: 
# (Data starts on 1986-01-03 runs until 2019-03-21)
#  -T10yAaa.tlag: Lag of yield difference between 10 year treasury and Aaa bonds
#  -AaaBaa.lag:  Lag of yield difference between Aaa bonds and Baa Bonds  
#  -T10yT1y.lag: Lag of yield differenc between 10 year and 1 year treasuries

#==============================================================================
# Auxiliary step:
# --------------
# 
# Redo the date column in order to extract date vectors that are used in the  
# plots for the paper later.  Note: This is not essential for estimation. 
# Equally it is not essential for plotting in general and could be done 
# differently. (It is given here to reproduce the figures as give in the paper.)
temp <- strptime(as.character(mult_data[,1]), "%Y-%m-%d")
temp <- as.integer(format(temp, "%Y%m%d"))
mult_data[,1] <- temp

rm(temp)
#==============================================================================


#--------------------------------------------------
#
# Step 2: Transform the data to correspond with the accompanying additive model
#         (Includes the removal of observations with missing values (NAs))
#
#------------------------------------------------------------------------------

source("AuxiliaryCode/data.additive.r")     
# Contains: data.additive <- function(YX)
#
# Removes missing values and adds transform of dependent variable
# needed for additive model as well as rescaled time covariate. 
#
# More details as comments in the file.

add_data <- data.additive(mult_data)



#------------------------------------------------------------------------------
#
# Step 3: The regressors are rescaled to the unit hypercube for 
#         use in the smooth backfitting algorithm.
#         (Includes the possibility to ensure that no gaps larger 
#         than specified in Bound are present --  not used here)
#
#------------------------------------------------------------------------------

source("AuxiliaryCode/scale.data.r")        
# Contains: scale.data <- function(YX,Bound)
#
# Rescales the covariates to lie in the unit interval.
# (By setting BOUND one can ensure that data are not too
# thinly distributed anywhere on the support -- is not used in
# application.) 
#
# More details as comments in the file.

Bound     <- c(1,1,1,1)            
data_est_input <- scale.data(add_data,Bound)



#------------------------------------------------------------------------------
#
# Step 4: Estimate the nonparametric components
# 
#------------------------------------------------------------------------------

# Step 4a: Extract inputs for the smoothbackfitting porcedure.
#-------------------------------------------------------------------------------

y_unit <- data_est_input$y_unit
# Dependent variable transformed for additive model 
# (the logarithm of the  squared returns).

x_unit <- data_est_input$x_unit 
# Covariates normalized to [0,1]^d. 

data_returns <- data_est_input$y.mult.sqrt
# Original dependent variable, i.e. the (log) returns 
# (for the observations used in estimating the nonparametric components). 

# Step 4b: Estimate the nonparametric components.
#-------------------------------------------------------------------------------

source("AuxiliaryCode/estimation.r")
# Contains: estimation <- function(x,y,bw.init,bw.method,Ngrid) 
#
# This calls the code to estimate the nonparametric components.
# It sources: 
#
# (i)    The shared library file "AuxiliaryCode/sbfnw.so" which contains the function
#        "sbfnw". (The original C code is in the same folder. This will need 
#        creating first first.)
# (ii)   "AuxiliaryCode/sbfnw.r", which calls the C code. 
#[(iii)] If the iterative bandwidth selection in chosen, then the corresponding
#        code in "AuxiliaryCode/bw_plugin_Loc_Vol.r" is sourced.  
#
# More details as comments in the file.

bw.init   <- c(0.01,0.05,0.05,0.05)
# Initial bandwidth vector for the smoothing procedure

bw.method <- "plugin"
# Choose the bandwidth selection procedure. Options: 
# (i)  "adhoc"  - do the estimation with the initial bandwidth vector or
# (ii) "plugin" - use the bandwidth selection procedure outlined in the supplement

Ngrid <- 200
# Specify number of equi-distant grid points to use in each covariate direction for the estimation

est_results <- estimation(x_unit,y_unit,bw.init,bw.method,Ngrid,data_returns)



#------------------------------------------------------------------------------
#
# Step 5: Estimate the parameters of the GARCH(1,1)
# 
#------------------------------------------------------------------------------

source("AuxiliaryCode/GARCH.r")                        
# Contains: estimation.GARCH <- function(est_results,x_unit,y_unit,data_returns) 
#
# Using the 'fGarch' package the GARCH parameters in the model are estimated. 
# In addition, the estimators for the GARCH parameters in:
# (i)  the Feng model with rescaled time as the only regressor model and 
# (ii) the basic model that models returns as a GARCH(1,1).
#
# As of R4.0.0 a warning message appears after running garchFit.
# Warning message:
# Using formula(x) is deprecated when x is a character vector of length > 1.
# Consider formula(paste(x, collapse = " ")) instead.
#
# More details as comments in the file.

Garch.estimates     <- estimation.GARCH(est_results,x_unit,y_unit,data_returns)

# Summary:
# (1) The nonparametric component estimates are contained in est_results
# (2) The GARCH component estimates are contained in Garch.estimates

















