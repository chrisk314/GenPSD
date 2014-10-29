#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define C4 0.498426033038178
#define VERY_BIG 1E200
#define VERBOSE_INTEGRATE

#include "GenPSD.h"

int Integrate(double *result, double a, double b, const int JMAX,
    const int N, const double Err, Func* MyFunc)
{

  double ss, ss_old = VERY_BIG; // Current and previous leading error estimate.
  double s[JMAX], temp[JMAX];   // First order estimates and temp storage.
  double del;                   // Integration subdivision width.
  
  del = (b-a)/N;
  s[0] = 0.0;
  for(int j=1; j<N; j++){
    s[0] += MyFunc->Eval(j*del);        // Add contributions for each integration.
  }
  s[0] = (del/2) * ( MyFunc->Eval(a) + 2*s[0] + MyFunc->Eval(b) );
  temp[0] = s[0];                       // Save result for truncation step.

  // Generate array of first order estimates.
  for(int k=2; k<=JMAX; k++){
    s[k-1] = 0.0;
    for(int j=0; j<N*pow(2,k-2); j++){
      s[k-1] += MyFunc->Eval(del/2 + j*del); // Add new contributions.
    }
    del = (b-a)/(N*pow(2,k-1));              // Halve interval for next step.
    s[k-1] = 0.5*s[k-2] + del*s[k-1];        // Combine with previous terms.
    temp[k-1] = s[k-1];                      // Save results for truncation step.

    // Truncate terms to improve estimate.
    for(int l=2; l<=k; l++){
      ss = (pow(4,l-1)*temp[k-1]-temp[l-2]) / (pow(4,l-1)-1); // Leading error term
      temp[l-2] = temp[k-1];                                  // Remove unnecessary terms.
      temp[k-1] = ss;                                         // Save new best estimate.  
      
      if(l == k){
	    // Check convergence.
        if(fabs(ss-ss_old) <= Err){
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
