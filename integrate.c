#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define C4 0.498426033038178


int Integrate(double *result, double a, double b, const int JMAX,
    const int N, const double Err, double (*func)(double)){

  int j,k,l;                             //k: refinement level, loop indices
  double ss, s[JMAX], temp[JMAX], del;   //leading error term, first order estimates
  
  //generate array of first order estimates
  for(k=1;k<=JMAX;k++){
    if(k==1){                          //include integration end points on first run
      del = (b-a)/(N*pow(2,k-1));
      for(s[k-1]=0.0,j=1;j<N*pow(2,k-1);j++){
        s[k-1]+=func(j*del);        //add up contributions for each integration
      }
      s[k-1] = (del/2)*(func(a) + 2*s[k-1] + func(b));
    }
    else{
      for(s[k-1]=0.0,j=0;j<N*pow(2,k-2);j++){
        s[k-1]+=func(del/2 + j*del);            //add up new contributions 
      }
      del = (b-a)/(N*pow(2,k-1));                //halve interval for next step
      s[k-1] = 0.5*s[k-2] + del*s[k-1];      //combine with terms from prev estimates 
    }

    printf("k: %d, I2^(%d)N: %lf\n", k, k-1, s[k-1]);
    printf("C(4): %lf\n\n", C4);

    temp[k-1] = s[k-1]; //save results for truncation step

    //truncate terms to improve estimate
    for(l=2;l<=k;l++){
      ss=(pow(4,l-1)*temp[k-1]-temp[l-2]) / (pow(4,l-1)-1); //calculate leading error term
      
      printf("k: %d, l: %d, ss: %lf\n", k, l, ss);
    
      if(l==k){
	    //compare against result of C4 from mathematica
	    if(fabs(ss - C4) <= Err){
	      //reached desired accuracy for estimate
	      printf("Estimate converged to desired accuracy %lf after %d refinements "\
                "with %d function evaluations.\n", Err, k, (int)(N*pow(2,k-1)+1));
	      printf("Estimate: I2^%dN = %lf\nActual: C(4)     = %lf\n", k, ss, C4);

	      *result = ss;
          return 0;
	    }
      }
      temp[l-2]=temp[k-1];      //remove old unnecessary terms
      temp[k-1]=ss;         //save new best error estimate    
    }
  }

  //convergence failed if we're here
  printf("Estimate failed to converge to accuracy %lf after %d function evaluations.\n",\
      Err, (int)(N*pow(2,k-2)));
  return 1;
}
