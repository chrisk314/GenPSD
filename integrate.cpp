#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define C4 0.498426033038178

#include "GenPSD.h"

int Integrate(double *result, double a, double b, const int JMAX,
    const int N, const double Err, Func* MyFunc)
{

  double ss, s[JMAX], temp[JMAX], del;   //leading error term, first order estimates
  
  del = (b-a)/N;
  s[0] = 0.0;
  for(int j=1; j<N; j++){
    s[0] += MyFunc->Eval(j*del);        //add up contributions for each integration
  }
  s[0] = (del/2) * ( MyFunc->Eval(a) + 2*s[0] + MyFunc->Eval(b) );
  temp[0] = s[0]; //save results for truncation step
#ifdef VERBOSE_INTEGRATE
  printf("k: 1, I2^(0)N: %+1.3e\n", s[0]);
#endif

  //generate array of first order estimates
  for(int k=2; k<=JMAX; k++){
    s[k-1] = 0.0;
    for(int j=0; j<N*pow(2,k-2); j++){
      s[k-1] += MyFunc->Eval(del/2 + j*del);            //add up new contributions 
    }
    del = (b-a)/(N*pow(2,k-1));                //halve interval for next step
    s[k-1] = 0.5*s[k-2] + del*s[k-1];      //combine with terms from prev estimates 
    temp[k-1] = s[k-1]; //save results for truncation step
#ifdef VERBOSE_INTEGRATE
    printf("k: %d, I2^(%d)N: %+1.3e\n", k, k-1, s[k-1]);
#endif

    //truncate terms to improve estimate
    for(int l=2; l<=k; l++){
      ss = (pow(4,l-1)*temp[k-1]-temp[l-2]) / (pow(4,l-1)-1); //calculate leading error term
#ifdef VERBOSE_INTEGRATE
      printf("k: %d, l: %d, ss: %+1.3e\n", k, l, ss);
#endif    
      if(l == k){
	    // Check convergence.
	    if(fabs(ss - C4) <= Err){
#ifdef VERBOSE_INTEGRATE
	      printf("Estimate converged to desired accuracy %+1.3e after %d refinements "\
                "with %d function evaluations.\n", Err, k, (int)(N*pow(2,k-1)+1));
	      printf("Estimate: I2^%dN = %+1.3e\nActual: C(4)     = %+1.3e\n", k, ss, C4);
#endif
	      *result = ss;
          return 0;
	    }
      }
      temp[l-2] = temp[k-1];  // Remove old unnecessary terms.
      temp[k-1] = ss;         // Save new best error estimate.  
    }
  }

#ifdef VERBOSE_INTEGRATE
  printf("Estimate failed to converge to accuracy %+1.3e after %d MyFunc->Evaltion evaluations.\n",\
      Err, (int)(N*pow(2,k-2)));
#endif

  return 1;

}
