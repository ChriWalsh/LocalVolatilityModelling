Date: 17.01.2021

This repository contains the R code to run the application in the 
conditionally accepted article:

Walsh, C.P.  and Vogt, M. "Local Volatility Modelling"
Journal of Business Economics and Statistics

The code is based on a C code (originally written by Berthold Haag). 

BEFORE the code can be run one needs to build a shared library 
which is then dynamically loaded:

(1) To create the dynamic library go to 
the "AuxiliaryCode/" folder and run "R CMD SHLIB sbfnw.c".
(2) The dynamic library is loaded on line 43 of 
"AuxiliaryCode/estimation.r", which reads 
"dyn.load("AuxiliaryCode/sbfnw.so","sbfnw")"
(For windows users this line will need adjusting to 
"dyn.load("AuxiliaryCode/sbfnw.dll","sbfnw")".

BEFORE the code can be run it may be necessary to install the 
"fGarch" package. 
(This currently seems to require R version 4.0.0 or above.) 

Once all the above preparation has been done then:


(I) The data is loaded and the estimates are calculated using 
"Main_Estimation_Code.r".

The code has been split into seperate steps with comments as to 
what each step does. The code for the corresponding functions are
all in the "AuxiliaryCode/" folder. 

(II) Additional code to get (pointwise) confidence intervals and 
construct plots and Table as in paper is provided in 
"Extra_Code_Plots.r".

(The code requires estimation of the model first as
mentioned in Step 0.) 

The code is then split into individual steps needed to 
reproduce the figures and the table from the paper. 


