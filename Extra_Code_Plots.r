#------------------------------------------------------------------------------
#
# Step 0: (Re-)Do the estimation
# 
#------------------------------------------------------------------------------

# Either:
# (i)  run the estimation including the bandwidth selection procedure in 
#      "Main_Estimation_Code.r" 
# or  
# (ii) run the estimation without the bandwidth selection 
#      (set ' bw.method <- "adhoc" ' in ) and use the previously selected 
#      bandwidth vector as the initial bandwidth vector. Here this is
#      bw.init <- c(0.1682183, 0.1918573, 0.2298709, 0.1804567)
# 
#
#
# The rest of Step 0 here follows (ii) and re-estiamtes the model without the 
# bandwidth selection procedure. 

rm(list=ls())
library(fGarch)

source("AuxiliaryCode/aux_functions1.r")
source("AuxiliaryCode/aux_functions2.r")

# Step 0.1: Read in the data used for the application 
mult_data <- read.csv("Data/Data_Application.csv")

temp <- strptime(as.character(mult_data[,1]), "%Y-%m-%d")
temp <- as.integer(format(temp, "%Y%m%d"))
mult_data[,1] <- temp
rm(temp)


# Step 0.2: Transform the data to correspond with the accompanying additive model
source("AuxiliaryCode/data.additive.r")  
add_data <- data.additive(mult_data)

# Step 0.3: The regressors are rescaled to the unit hypercube for 
#           use in the smooth backfitting algorithm.
source("AuxiliaryCode/scale.data.r")
Bound     <- c(1,1,1,1)            
data_est_input <- scale.data(add_data,Bound)

# Step 0.4: Estimate the nonparametric components
y_unit <- data_est_input$y_unit
x_unit <- data_est_input$x_unit 
data_returns <- data_est_input$y.mult.sqrt
source("AuxiliaryCode/estimation.r")

# The next two lines are the only ones that deviate from "Main_estimation_Code.r"
bw.init   <- c(0.1682183, 0.1918573, 0.2298709, 0.1804567)
bw.method <- "adhoc"

Ngrid <- 200
est_results <- estimation(x_unit,y_unit,bw.init,bw.method,Ngrid,data_returns)

# Step 0.5: Estimate the parameters of the GARCH(1,1)
source("AuxiliaryCode/GARCH.r")      
Garch.estimates     <- estimation.GARCH(est_results,x_unit,y_unit,data_returns)







#------------------------------------------------------------------------------
#
# Step 1: Construct the graphs of the data used in Figure 1 of the paper
# 
#------------------------------------------------------------------------------

# Extract a date vector for the plots
date.plot  <- as.Date(as.character(data_est_input$id),format="%Y%m%d")
date.plot  <- seq(date.plot[1],date.plot[length(date.plot)],by=(date.plot[length(date.plot)] - date.plot[1])/(dim(data_est_input$x.original)[1]-1))

# Plot the log returns
plot(date.plot,data_returns,ylab="Log Returns",xlab="Time",type="l")
# Plot the covariates used in the application
plot(date.plot,data_est_input$x.original[,2],ylab="Yield Difference (10yr bill - Aaa)",xlab="Time",type="l")
plot(date.plot,data_est_input$x.original[,3],ylab="Yield Difference (Baa - Aaa)",xlab="Time",type="l")
plot(date.plot,data_est_input$x.original[,4],ylab="Yield Difference (10yr - 1 yr bill)",xlab="Time",type="l")







#------------------------------------------------------------------------------
#
# Step 2: Construct the graphs of the nonparametric estimates 
#         given in Figure 3 of the paper.
# 
#------------------------------------------------------------------------------

# Step 2a: Construct matric containing the support of estimation points
#------------------------------------------------------------------------------
d     <- dim(est_results$sbf)[2]                      
Ngrid <- dim(est_results$sbf)[1]                                   

x.plot <- matrix(0,ncol=dim(data_est_input$x_unit)[2],nrow=Ngrid)
for(k in 1:dim(data_est_input$x_unit)[2])
  x.plot[,k] <- seq(data_est_input$x_lower[k],data_est_input$x_upper[k],by=(data_est_input$x_upper[k] - data_est_input$x_lower[k])/(Ngrid-1))

# Step 2b: Construct the (pointwise) confidence bands for the 
# nonparametric estimates.
#------------------------------------------------------------------------------

source("AuxiliaryCode/confid.r")                       
# Contains: confidence.bands <- function(est_results,x_unit,y_unit)
#
# Returns (pointwise) confidence bands for the nonparametric components in the
# additive as well as the multiplicative model

conf <- confidence.bands(est_results,x_unit,y_unit)



# Step 2c: Construct the normalized estimates and confidence bands given in
#           Figure 3 of paper.
#------------------------------------------------------------------------------

# Set the normalization point to the median of the covariates
L_x               <- length(x_unit[1,]) - 1       
x.normalize.point <- matrix(0,ncol=L_x,nrow=1)

for(i in 1:L_x)
  x.normalize.point[1,i] <- median(x_unit[,i+1])  

source("AuxiliaryCode/normalize.r")   
# Contains: normalize.tau.functions <- function(est_results,x_unit,y_unit,x.normalize.point,conf,data_returns)
#
# Returns normalized estimates of the multiplicative components and their 
# (pointwise) confidence bands.

tau.sqrd <- normalize.tau.functions(est_results,x_unit,y_unit,x.normalize.point,conf,data_returns)

# Step 2d: Plot the graphs for Figure 3 of the paper.
#------------------------------------------------------------------------------

# Determine plotting range
min.plot <- vector(mode="double",length=d)
max.plot <- vector(mode="double",length=d)

for(k in 1:d)
{   min.plot[k] <- min(tau.sqrd$tau.sqrd.lower[,k]) 
max.plot[k] <- max(tau.sqrd$tau.sqrd.upper[,k]) 
}  

# Get expressions for the Y-axis labels
ylabel.mult <- c(expression(tilde(tau)[0]^2),expression(tilde(tau)[1]^2),expression(tilde(tau)[2]^2),expression(tilde(tau)[3]^2),
                 expression(tau[4]^2),expression(tau[5]^2),expression(tau[6]^2))

# Determine the position for the normalization point, i.e. undo the scaling
normal.line <- rep(0,d-1)
for(k in 1:d-1)
  normal.line[k] <- x.normalize.point[k]*(data_est_input$x_upper[k+1]-data_est_input$x_lower[k+1])+data_est_input$x_lower[k+1]

# Adjust the lower bound to fit in the marker for the normalization point 
# and adjust the upper bound so that the fits for the first and second regressor 
# are on the same scale
min.plot[-1] <- rep(0,3)
max.plot[-1] <- rep(max(max.plot[-1]),3)

# Construct the plot for the first regressor 
plot(x.plot[,2],tau.sqrd$tau.sqrd[,2],type="l",lwd=1.5,ylim=c(min.plot[2],max.plot[2]),xlab="Spread (in pp) between Aaa and 10yr treasuries",ylab=ylabel.mult[2])
lines(x.plot[,2],tau.sqrd$tau.sqrd.upper[,2],lty=2,lwd=1.5)
lines(x.plot[,2],tau.sqrd$tau.sqrd.lower[,2],lty=2,lwd=1.5) 
abline(h=1,lty=3,lwd=1.5)
points(x=normal.line[1],y=min.plot[2],pch=6)

# Construct the plot for the second regressor
plot(x.plot[,3],tau.sqrd$tau.sqrd[,3],type="l",lwd=1.5,ylim=c(min.plot[3],max.plot[3]),xlab="Spread (in pp) between Baa and Aaa bonds",ylab=ylabel.mult[3])
lines(x.plot[,3],tau.sqrd$tau.sqrd.upper[,3],lty=2,lwd=1.5)
lines(x.plot[,3],tau.sqrd$tau.sqrd.lower[,3],lty=2,lwd=1.5) 
abline(h=1,lty=3,lwd=1.5)
points(x=normal.line[2],y=min.plot[3],pch=6)

# Construct the plot for the third regressor. (Change the upper plot limit to 
# 2 as in paper.)
plot(x.plot[,4],tau.sqrd$tau.sqrd[,4],type="l",lwd=1.5,ylim=c(min.plot[4],2),xlab="Spread (in pp) between 10yr and 1yr treasuries",ylab=ylabel.mult[4])
lines(x.plot[,4],tau.sqrd$tau.sqrd.upper[,4],lty=2,lwd=1.5)
lines(x.plot[,4],tau.sqrd$tau.sqrd.lower[,4],lty=2,lwd=1.5) 
abline(h=1,lty=3,lwd=1.5)
points(x=normal.line[3],y=min.plot[4],pch=6)







#------------------------------------------------------------------------------
#
# Step 3: Plot the estimation results for the trend component given in Figure
#         2 of the paper.
# 
#------------------------------------------------------------------------------

# Construct a date vector for the estimation grid
date.plot2 <- as.Date(as.character(data_est_input$id),format="%Y%m%d")
date.plot2 <- seq(date.plot2[1],date.plot2[length(date.plot2)],by=(date.plot2[length(date.plot2)] - date.plot2[1])/(Ngrid-1))

# Construct the plot for the component that depends on (rescaled) time
plot(date.plot2,tau.sqrd$tau.sqrd[,1],type="l",lwd=1.5,ylim=c(min.plot[1],max.plot[1]),xlab="Time",ylab=ylabel.mult[1])
lines(date.plot2,tau.sqrd$tau.sqrd.upper[,1],lty=2,lwd=1.5)
lines(date.plot2,tau.sqrd$tau.sqrd.lower[,1],lty=2,lwd=1.5) 




#------------------------------------------------------------------------------
#
# Step 4: Collect GARCH parameter estimates as in Table 1 of the paper.
# 
#------------------------------------------------------------------------------

para <- list(Garch.estimates$coef.Garch.a,Garch.estimates$coef.Garch.trend.only,Garch.estimates$coef.Garch.simple)

GARCH.parameters           <- matrix(0,ncol=5,nrow=3)
colnames(GARCH.parameters) <- c("Model","omega","alpha","beta","Persistence")
GARCH.parameters[1,]       <- c("Simple",para[[3]],para[[3]][2]+para[[3]][3]) 
GARCH.parameters[2,]       <- c("Trend only",para[[2]],para[[2]][2]+para[[2]][3]) 
GARCH.parameters[3,]       <- c("Trend and Covariats",para[[1]],para[[1]][2]+para[[1]][3]) 

# Calculate the half life of volatility
HalfLife <- as.numeric(GARCH.parameters[,3])+as.numeric(GARCH.parameters[,4]) 
HalfLife <- 1 - log(2)/log(HalfLife)
GARCH.table <- cbind(GARCH.parameters,HalfLife) 

# Print the table containing the values to construct Table 1 in the paper.
GARCH.table







#------------------------------------------------------------------------------
#
# Step 5: Plot the unconditional volatility levels as in right panel of Figure 2.
# 
#------------------------------------------------------------------------------

# Estimate unconditional variance level model with trend only
trend.only   <- Garch.estimates$trend.only*Garch.estimates$var.eps.trend.only

# Set the limits for the plot
max.plot <- max(trend.only,tau.sqrd$tau0.sqrd)
min.plot <- min(trend.only,tau.sqrd$tau0.sqrd)


# Construct the plot of the unconditional volatility levels
plot(date.plot2,trend.only,type="l",lty=2,ylim=c(min.plot,max.plot+0.00010),xlab="Date",ylab="")
lines(date.plot2,tau.sqrd$tau0.sqrd)
legend("top",c("Model with covariates","Trend only model"),lty=c(1,2))

