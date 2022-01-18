#include <math.h> 
#include <Rinternals.h> 
#include <R.h>
#include <stdlib.h> 
#include <Rmath.h>


  double awert(double x)
  { if(x < 0) 
      return(-x);
    else
      return(x);
  }


  
  double epanc(double x, double xi, double h)
  {       double lower, upper;
          if(xi > 1 || xi < 0)
	     return(0.0);
          else
	  {  if(awert(x-xi) < h) 
	     {  if(1-xi < h)  
                   upper = (1-xi) * (1 - (1-xi) * (1-xi)/(h*h*3.0))/h;
                else 
                   upper = 2.0/3.0;
                if(xi < h)
                   lower = -xi * (1 - xi * xi/(h*h*3.0))/h;
                else
                   lower = -2.0/3.0;
                return((1- (x-xi) * (x-xi)/(h*h))/(upper-lower)) ;
             }
             else
                return(0.0);
          }
  }



  double biweight(double x, double xi, double h)
  {       double lower, upper;
          if(xi > 1 || xi < 0)
	     return(0.0);
          else
	  {  if(awert(x-xi) < h) 
	     {  if(1-xi < h)  
		 upper = ((1-xi)/h) * (1 - (2.0/3.0) * pow((1-xi)/h,2) +  (1.0/5.0) * pow((1-xi)/h,4));
                else 
                   upper = 8.0/15.0;
                if(xi < h)
                   lower =  (-xi/h) * (1 - (2.0/3.0) * pow(-xi/h,2) +  (1.0/5.0) * pow(-xi/h,4));
                else
                   lower = -8.0/15.0;
                return((1 - pow((x-xi)/h,2)) * (1 - pow((x-xi)/h,2)) / (upper-lower));
             }
             else
                return(0.0);
          }
  }



  void sbfnw(double *x, double *y, int *n, int *dim, double *grid, int *ng, double *h, 
             double *dens, double *m0, double *mbf, double *reg, int *iterate, double *conv, 
             double *fh_int, double *fh_int2a, double *fh_int2b, double *reg_int)


             /* x        regressor x= (x_1^1,...,x_n^1,...,x_1^d,...,x_n^d)
                y        predictor y= (y_1,...,y_n)
                n        sample size
                dim      dimension of regressors d
                grid     grid vector grid= (g_1^1,...,g_N^1,...,g_1^d,...,g_N^d)
                ng       length of grid N
                h        vector of bandwidths (d-dim)
                dens     kernel density estimators
                m0       estimate of model constant
                mbf      backfitting estimators 
                         mbf= (m^1(g_1^1),...,m^1(g_N^1),...,m^d(g_1^d),...,m^d(g_N^d))          
                reg      vector of one-dimensional NW estimators
                iterate  max. number of iterations, 
                         retour: actual number of iterations
                conv     value of convergence criterion 
                fh_int   integrals of one-dimensional kernel density estimators
                fh_int2a two-dimensional kernel density estimators with the lower index
                         component integrated out (should be approximately equal to one-
                         dimensional kernel density estimator)       
                fh_int2b two-dimensional kernel density estimators with the higher index
                         component integrated out (should be approximately equal to one-
                         dimensional kernel density estimator)
		reg_int  integrals of the one-dimensional NW estimators */


  {     double kerx, dummy, dummy1, dummy2, dummy3, dummy4, critz, critn, delta, n0, ymean;
    double *weight, *fh, *fh_inv, *rh, *s00, *mbf0, *mean_mbf, *ind, *m0_j;
        int i, j, k, l, m, init1, init2, stop1, stop2, pos, posy, pos1, pos2, pos3, pos4, abbruch;

        
        
        /* Definition and initializing of some variables */
       
        weight = (double*) malloc(sizeof(double) * dim[0] * n[0] * ng[0]);
        rh = (double*) malloc(sizeof(double) * dim[0] * ng[0]);
        fh = (double*) malloc(sizeof(double) * dim[0] * ng[0]);
        fh_inv = (double*) malloc(sizeof(double) * dim[0] * ng[0]);
        mean_mbf = (double*) malloc(sizeof(double) * dim[0]);
        mbf0 = (double*) malloc(sizeof(double) * dim[0] * ng[0]);                   
        m0_j = (double*) malloc(sizeof(double) * dim[0]);       
        ind = (double*) malloc(sizeof(double) * n[0]);        


        for(i=0; i< dim[0]* ng[0]; i++)
        {   rh[i] = 0;
            fh[i] = 0;
        }

        for(i=0; i < *n; i++)
	    ind[i] = 1;
        
        n0 = 0;

        for(i=0; i < *n; i++)
        {   for(k=0; k < *dim; k++)
            {   if(x[i + k*n[0]] > 1 || x[i + k*n[0]] < 0)
                   ind[i] = 0;
	    }
	    n0 += ind[i];
        }

        ymean = 0;

        for(i=0; i < *n; i++)
	{   ymean += ind[i] * y[i];
        }

        ymean = ymean/n0;
        m0[0] = ymean;   

        
        
        /* Calculation of weights */

        init1 = 0;
        init2 = 0;
        stop1 = *n;
        stop2 = *ng;
        pos = 0;

        for(k=0; k< *dim; k++)
        {   posy = 0;
            for(i=init1; i< stop1; i++)
            {      for(j=init2; j< stop2; j++)
                   {    /* kerx = dnorm((grid[j] - x[i])/h[k],0,1,0); */
                         kerx = epanc(grid[j], x[i], h[k]); 
                         weight[pos] = kerx;                       
                         
                         /* weight = (w_11^1,...,w_1N^1,...,w_n1^1,...,w_nN^1,...,w_nN^d) */
                                                                  
                         rh[j] += ind[posy] * y[posy] * kerx;  /* Calculation of numerator of NW estimator */
                         fh[j] += ind[posy] * kerx;            /* Calculation of denominator of NW estimator */
                         pos++;
                   }
                   posy++;
            }
            init1 += *n;
            init2 += *ng;
            stop1 += *n;
            stop2 += *ng;
        }
         
   
        /* Calculation of one-dimensional estimators at the grid points */

        stop1 = ng[0] * dim[0];

        for(i=0; i< stop1; i++)
	{   if(fh[i] < 0.0001)
               fh_inv[i] = 0.0;
            else
               fh_inv[i] = 1.0/fh[i]; 
            reg[i] = rh[i] * fh_inv[i];
            mbf[i] = reg[i];
        }


        /* Calculation of kernel density estimators */

        for(k=0; k < *dim; k++)
	{   for(i=0; i < *ng; i++)
	       dens[i + k*ng[0]] = fh[i + k*ng[0]] / (n0 * h[k]);
               /* Note: have to rescale the bandwidth h[k] if the real
                  data do not come from the unit interval [0,1] */
        } 


        /* Normalization of starting value for the backfitting algorithm */  

        for(i=0; i < *dim; i++)
	{   dummy1 = 0; 
            for(j=0; j < *ng - 1; j++)
	    {   delta = grid[j + 1 + i*ng[0]] - grid[j + i*ng[0]];
                dummy1 += mbf[j + i*ng[0]] * dens[j + i*ng[0]] * delta;                           
            }
            for(j=0; j < *ng; j++)
 	        mbf[j + i*ng[0]] = mbf[j + i*ng[0]] - dummy1;
        }
       
              
        /* Calculation of two-dimensional kernel density estimators at the grid points 
           s00  = (s_12, s_13,s_14,...,s_1d,...s_d-1,d)
           s_12 = (s_12(g_1^1,g_1^2),...,s_12(g_1^1,g_N^2),...,s_12(g_N^1,g_N^2)) */
        
        s00 = (double*) malloc(sizeof(double) * dim[0] * (dim[0] - 1)* ng[0] * ng[0] /2);
        
        posy = 0;
        
        for(i=0; i< *dim-1; i++)
        {   init1 = i*n[0]*ng[0];
            for(j=i+1; j< *dim; j++)
            {   init2 = j* n[0] * ng[0];
                for(k=0; k< *ng; k++)
                {   for(l=0; l< *ng; l++)
                    {   s00[posy] = 0;
                        for(m=0; m< *n; m++)
                        {   s00[posy] += ind[m] * weight[init1 + m* ng[0] + k] * weight[init2 + m* ng[0] + l];}
                        posy++;
                    }
                }
            }
        }
        
      
        free(weight);
 
    

        /* Calculation of some integral expressions for the SBF algorithm 
           fh_int2a and fh_int2b are ordered in analogy to the vector s00 */


	for(i=0; i < *dim; i++)
	{   dummy1 = 0;
	    dummy2 = 0;
            dummy3 = 0;
            dummy4 = 0;
            for(j=0; j < *ng-1; j++)
	    {   delta = grid[j + 1 + i*ng[0]] - grid[j + i*ng[0]];
                dummy1 += fh[j + i*ng[0]] * delta;
	        dummy2 += fh[j + 1 + i*ng[0]] * delta;
                dummy3 += reg[j + i*ng[0]] * fh[j + i*ng[0]] * delta;
	        dummy4 += reg[j + 1 + i*ng[0]] * fh[j + 1 + i*ng[0]] * delta;
            }
            fh_int[i]  = (dummy1 + dummy2) / (2.0*n0*h[i]);
            reg_int[i] = (dummy3 + dummy4) / (2.0*n0*h[i]); 
        }  

        
        for(i=0; i< *dim; i++)
	{   if(fh_int[i] < 0.0001)
               fh_int[i] = 0.0;
            else
	       fh_int[i] = 1.0/fh_int[i]; 
            m0_j[i] = reg_int[i] * fh_int[i];
        }        

	       
        pos = 0;
        init1 = 0;
        stop1 = ng[0];
        

        for(l=0; l < *dim-1; l++)
	{  for(i=l+1; i < *dim; i++) 
	   {   for(j=init1; j < stop1; j++)
               {   dummy1 = 0;
	           dummy2 = 0;
                   for(k=0; k < *ng-1; k++)
		   {   delta = grid[k + 1 + l*ng[0]] - grid[k + l*ng[0]];
                       dummy1 += s00[j + k*ng[0]] * delta;
		       dummy2 += s00[j + (k+1)*ng[0]] * delta;                                
                   }
	           fh_int2a[pos] = ((dummy1 + dummy2) / (2.0*n0*h[i]*h[l])) * fh_int[l];    
                   pos++;
               }
	       init1 += ng[0]*ng[0];
               stop1 += ng[0]*ng[0];
           }
        }
       
        
	pos = 0;
        pos1 = 0;
        init1 = 0;
        stop1 = ng[0];


	for(i=0; i < *dim-1; i++)
	{  for(l=i+1; l < *dim; l++) 
	   {   for(j=0; j < *ng; j++)
               {   dummy1 = 0;
	           dummy2 = 0;
                   for(k=0; k < *ng-1; k++)
		   {   delta = grid[k + 1 + l*ng[0]] - grid[k + l*ng[0]];
                       dummy1 += s00[pos1*ng[0]*ng[0] + j*ng[0] + k] * delta;
	               dummy2 += s00[pos1*ng[0]*ng[0] + j*ng[0] + k + 1] * delta;
                   }
	           fh_int2b[pos] = ((dummy1 + dummy2) / (2.0*n0*h[i]*h[l])) * fh_int[l];    
                   pos++;
               }
               pos1++;
	       init1 += ng[0]*ng[0];
               stop1 += ng[0]*ng[0];
           }
        } 


   
	/* Smooth backfitting iteration */
        
        m=0;
        do
        {
          
          /* Initializing */
         
          pos1=0;
          pos2=0;
          posy=0;
          abbruch = 1; 
          
        
          m++;
          
          /* Calculation of the SBF estimator m^i at the grid points */
              
          for(i=0; i< *dim ; i++)            
          {   critz = 0;
              critn = 0;
              for(j=0; j< *ng; j++)              
              {   
                  /* Calculation of the integrals that are substracted */
                  
                  dummy=0;
                  pos=0;
                  
                  /* Integrals for k < i */
                  
                  for(k=0; k< i; k++)           
                  {   dummy1 =0;
		      dummy2 =0; 
                      pos1 = (k * (dim[0] - 1) - k*(k-1)/2 + (i - k) - 1) * ng[0] * ng[0] + j;
                      pos3 = (k * (dim[0] - 1) - k*(k-1)/2 + (i - k) - 1) * ng[0];
                      for(l=0; l< *ng-1; l++)
		      {   delta = grid[l + 1 + k*ng[0]] - grid[l + k*ng[0]];
                          dummy1 += mbf[pos] * ( s00[pos1] * fh_inv[j + i*ng[0]] / h[k] - fh_int2b[pos3] ) * delta;
			  dummy2 += mbf[pos+1] * ( s00[pos1 + *ng] * fh_inv[j + i*ng[0]] / h[k] - fh_int2b[pos3+1] ) * delta;
                          pos++;
                          pos1+= *ng;
                          pos3++;                        
                      }
                      pos++;
                      dummy += (dummy1 + dummy2) / 2.0;
                  } 
                 
                  pos+= *ng;                  /* k = i is omitted */
                  
                  /* Integrals for k > i */
                
                  for(k=i+1; k< *dim; k++)      
                  {   dummy1 =0;
		      dummy2 =0;
                      pos2 = (dim[0] * i - i * (i+1)/2 + k-i-1) * ng[0] * ng[0] + j * ng[0];
                      pos4 = (dim[0] * i - i * (i+1)/2 + k-i-1) * ng[0];
                      for(l=0; l< *ng-1; l++)
		      {   delta = grid[l + 1 + k*ng[0]] - grid[l + k*ng[0]];
                          dummy1 += mbf[pos] * ( s00[pos2] * fh_inv[j + i*ng[0]] / h[k] - fh_int2a[pos4] ) * delta;
			  dummy2 += mbf[pos+1] * ( s00[pos2+1] * fh_inv[j + i*ng[0]] / h[k]- fh_int2a[pos4+1] ) * delta;
                          pos++;
                          pos2++;
                          pos4++;
                      }
                      pos++;
                      dummy += (dummy1 + dummy2) / 2.0;
                  }
                                   
                  mbf0[posy] = mbf[posy];                    
                  mbf[posy] = reg[posy] - dummy - m0_j[i];          
                  
                  posy++; 
              }

              /* Normalization of the SBF estimator */

              dummy1 = 0;
              dummy2 = 0;

              for(j=0; j < *ng-1; j++)
	      {   delta = grid[j + 1 + i*ng[0]] - grid[j + i*ng[0]];
                  dummy1 += mbf[j + i*ng[0]] * fh[j + i*ng[0]] * delta / (n0*h[i]);
		  dummy2 += mbf[j + 1 + i*ng[0]] * fh[j + 1 + i*ng[0]] * delta / (n0*h[i]);      
              }

              mean_mbf[i] = (dummy1 + dummy2) / 2.0;
              
              for(j=0; j < *ng; j++)
              {   mbf[j + i*ng[0]] = mbf[j + i*ng[0]] - mean_mbf[i];
              }
  
              /* Convergence criterion */

              for(j=0; j < *ng; j++)
	      {   critz += (mbf[j + i*ng[0]] - mbf0[j + i*ng[0]]) * (mbf[j + i*ng[0]] - mbf0[j + i*ng[0]]);
		  critn += mbf0[j + i*ng[0]] * mbf0[j + i*ng[0]];
              }                     
              
              conv[i] = critz/(critn+0.00001);

              if( critz/(critn+0.00001) > 0.00001) 
                  abbruch = 0;
              else 
                  abbruch *=1;
          }
        
        }

        while(abbruch==0 && m < *iterate);
        

        iterate[0] = m;
        
        free(fh);
        free(s00);     
  }

         
   
