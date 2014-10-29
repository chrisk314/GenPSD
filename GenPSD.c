#include "GenPSD.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
void BetaDist(double beta_a, double beta_b, int NumElem, double *h)
{
  // Calculate prefactor for this choice of a and b.
  double B_factor;
  B_factor = tgamma(beta_a + beta_b) / ( tgamma(beta_a) * tgamma(beta_b) );
  
  // Perform integral to calculate grading curve values.
  h[0] = 0;
  for(int i=1; i<NumElem; i++){
    h[i] = i*2.0; 
  }
}


double ReducedDiam(double *d_min, double *d_max, double *d){
    double d_r = (d - d_min)/(d_max - d_min);
    return d_r;
}

