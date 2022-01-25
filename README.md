# LocalVolatilityModelling

Date: 21.01.2022

This repository contains the R code to run the application in the 
conditionally accepted article:

Walsh, C.P.  and Vogt, M. (2022+) "Locally Stationary Multiplicative Volatility Modelling"
*Journal of Business Economics and Statistics*

The code used for computing the smooth backfitting estimators is based on C code (originally written by Berthold Haag).

## (Very) brief overview of repository

The repository contains:

* the data in the folder **Data/**
* the main estimation code in the file **Main_Estimation_Code.r** to estimate the model components  
* additional code to reproduce the figures and the table from the application section of the paper is in the file **Extra_Code_Plots.r**
* auxiliary commented functions in the folder **AuxiliaryCode/**

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
