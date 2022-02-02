# LocalVolatilityModelling

Date: 02.02.2022

This repository contains the data set and the R code to run the application in:

Walsh, C.P.  and Vogt, M. (2022), "Locally Stationary Multiplicative Volatility Modelling"
*Journal of Business Economics and Statistics*, Online at https://doi.org/10.1080/07350015.2022.2036612 

The code used for computing the smooth backfitting estimators is based on C code (originally written by Berthold Haag).

## (Very) brief overview of repository

The repository contains:

* the data in the folder **Data/**
* the main estimation code in the file **Main_Estimation_Code.r** to estimate the model components  
* additional code to reproduce the figures and the table from the application section of the paper is in the file **Extra_Code_Plots.r**
* auxiliary commented functions in the folder **AuxiliaryCode/**

### Data 

The data in the folder **Data/Data_Application.csv** consists of the following:

* The first column (*Date*) contains the date 
* The second column (*SP500*) contains the S&P 500 log returns
* The third column (*T10yAaa.lag*) contains the lagged yield difference between 10 year treasuries and Aaa bonds
* The fourth column (*Aaa.Baa.lag*) contains the lagged yield difference between Aaa bonds and Baa bonds
* The fifth column (*T10yT1y.lag*) contains the lagged yield difference between 10 year and 1 year treasuries

The raw data used to calculate the S&P 500 log returns were taken from Yahoo!Finance (https://finance.yahoo.com). The raw data used to calculate the lagged yield differences were taken from the FRED database of the Federal Reserve bank of St.Louis (https://fred.stlouisfed.org/).


## Preliminariy steps BEFORE running the code

1. In order to run the C code one needs to build a shared library 
which is then dynamically loaded:
+ To create the dynamic library: Start R; go to 
the **AuxiliaryCode/** folder and run `R CMD SHLIB sbfnw.c`.
+ The dynamic library is loaded on line 43 of 
**AuxiliaryCode/estimation.r**, which reads 
`dyn.load("AuxiliaryCode/sbfnw.so","sbfnw") ` 
For windows users this line will need adjusting to 
`dyn.load("AuxiliaryCode/sbfnw.dll","sbfnw")`. 

2. The estimation of the GARCH parameters is done with the *fGarch* package, so this will need to be installed, which currently requires R version 4.0.0 or above.


## More on the code

Once all the above preliminary preparations have been done then, the code can be run:

1. The data is loaded and the estimates are calculated when running 
**Main_Estimation_Code.r**. The code has been split into seperate steps with comments as to 
what each step does. The code for the corresponding functions are
all in the **AuxiliaryCode/** folder.  

2. Additional code to get (pointwise) confidence intervals and 
construct figures and table as in paper is provided in 
**Extra_Code_Plots.r**. 
(The code requires estimation of the model first as
mentioned in Step 0.) 
The code is then split into individual steps needed to 
reproduce the figures and the table from the paper again the 
functions used are in **AuxiliaryCode/** folder. 
