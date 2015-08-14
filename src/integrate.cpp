#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define VERY_BIG 1E200
//#define VERBOSE_INTEGRATE

#include "GenPSD.h"

int Integrate(double *result, double a, double b, const int KMAX,\
    const int N, const double Err, Func* MyFunc)
{

  double ss, ss_old = VERY_BIG; // Current and previous leading error estimate.
  double s[KMAX], temp[KMAX];   // First order estimates and temp storage.
  double del;                   // Integration subdivision width.
  
  del = (b-a)/N;
  s[0] = 0.0;
  for(int j=1; j<N; j++){
    s[0] += MyFunc->Eval(j*del);        // Add contributions for each integration.
  }
  s[0] = (del/2) * ( MyFunc->Eval(a) + 2*s[0] + MyFunc->Eval(b) );
  temp[0] = s[0];                       // Save result for truncation step.

  // Generate array of first order estimates.
  for(int k=1; k<=KMAX; k++){
    s[k] = 0.0;
    for(int j=0; j<N*pow(2,k-1); j++){
      s[k] += MyFunc->Eval(del/2 + j*del); // Add new contributions.
    }
    del = (b-a)/(N*pow(2,k));              // Halve interval for next step.
    s[k] = 0.5*s[k-1] + del*s[k];          // Combine with previous terms.
    temp[k] = s[k];                        // Save results for truncation step.

    // Truncate terms to improve estimate.
    for(int l=1; l<=k; l++){
      ss = (pow(4,l)*temp[k]-temp[l-1]) / (pow(4,l)-1); // Leading error term
      temp[l-1] = temp[k];                              // Remove unnecessary terms.
      temp[k] = ss;                                     // Save new best estimate.  
      
      if(l == k){
	    // Check convergence.
        if(fabs(ss-ss_old) <= Err){
#ifdef VERBOSE_INTEGRATE
          printf("Estimate converged to accuracy %+1.3e in %d function evaluations.\n",\
                Err, int(N*pow(2,k)+1));
#endif
	      *result = ss;
          return 0;
	    }
        ss_old = ss;
      }
    }
  }

#ifdef VERBOSE_INTEGRATE
  printf("Estimate failed to converge to accuracy %+1.3e\n", Err);
#endif

  return 1;

}
