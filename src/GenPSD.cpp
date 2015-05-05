#include "GenPSD.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//#include <boost/random.hpp>
//#include <boost/random/normal_distribution.hpp>
//#include <boost/generator_iterator.hpp>

#include "rand_gen.h"
#include "integrate.h"

#define VERY_BIG 1E200
#define VERY_BIG_INT 60000

Func::Func(double& _a, double& _b):
    a(_a), b(_b){};

double Func::Eval(double x){
    //return 3*x*x;
    //return cos(M_PI*x*x/2);
    return pow(x, a-1) * pow(1-x,b-1);
}


int BetaDist(double beta_a, double beta_b, int NumElem, double *dia, double *h){

  double B_factor;      // Prefactor for this choice of a and b.
  double IntStart=0.0, IntStop=1.0, x;
  int IntIter=10, IntSubDiv=26;
  double IntErr=1E-10;
  Func *MyFunc;

  // Set up integrand function and gamma fn. prefactor.
  MyFunc = new Func(beta_a, beta_b);
  B_factor = tgamma(beta_a + beta_b) / ( tgamma(beta_a) * tgamma(beta_b) );
  
  // Perform integral to calculate grading curve values.
  h[0] = 0;
  for(int i=1; i<NumElem; i++){
    x = (dia[i] - dia[0]) / (dia[NumElem-1] - dia[0]);
    if( Integrate(&h[i], IntStart, x, IntIter, IntSubDiv, IntErr,\
            MyFunc) != 0 )
    {
      printf("Integral did not converge in function Integrate called from "\
                "function BetaDist. Program will exit.\n");
      return 1;
    }
    h[i] *= B_factor; 
  }
  
  return 0;
}


double ReducedDiam(double *d_min, double *d_max, double *d){
  double d_r = (d - d_min)/(d_max - d_min);
  return d_r;
}


double GenDia(double d_lower, double d_upper, double y){
    return (d_lower*d_upper) / (d_upper + y*(d_lower - d_upper));
}


int GenPopulations(int *N_tot, int N_c, int N_pc_min, int N_p_min, int *N,\
                    double *h, double *dia){
    
  double MeanDia;
  double VolMeanDia;                    // Mean volume for given size class.
  double TotalVolEst;
  double SpherePrefactor;               
  double N_try[N_c];                    // Estimated particle no. at intermediate steps.
  double N_min=VERY_BIG;
  int SumParticles;

  SpherePrefactor = (1.0/6.0) * M_PI;

  // Estimate total volume.
  //MeanDia = GenDia(dia[0], dia[N_c], 0.5);// incorrect! Only correct when grading is uniform by volume fraction!
  for(int i=0; i<=N_c; i++){
	 if(h[i] > 0.5){
		 MeanDia = dia[i-1] + ((0.5 - h[i-1]) / (h[i] - h[i-1])) * (dia[i] - dia[i-1]);
		 break;
	 }
  }
  VolMeanDia = SpherePrefactor * MeanDia * MeanDia * MeanDia;
  TotalVolEst = VolMeanDia * N_p_min;

  // Get first estimate for size class populations from grading curve h(d).
  printf("\nGrading class trial populations:\n");
  SumParticles = 0;
  for(int i=0; i<N_c; i++){
    MeanDia = GenDia(dia[i], dia[i+1], 0.5);
    VolMeanDia = SpherePrefactor * MeanDia * MeanDia * MeanDia;
    N_try[i] = (h[i+1] - h[i]) * TotalVolEst / VolMeanDia ; // Using volume esimate.

    N[i] = int(N_try[i]);
    SumParticles += N[i];
    if(N_try[i] < N_min){
      N_min = N_try[i];
    }
    printf("N_try[%d] = %1.4lf\n", i, N_try[i]);
  }

  // Check that all size classes contain at least the minimum required number
  // of particles N_pc_min.
  if(N_min == 0){
    printf("Encountered empty size class. Program will exit.\n");
    exit(-1);
  }
  else if(N_min < N_pc_min){
    double RescaleFactor = double(N_pc_min)/N_min;
    SumParticles = 0;
    for(int i=0; i<N_c; i++){
      N[i] = int(RescaleFactor*N_try[i]);
      SumParticles += N[i];
    }
  }

  // Check total number of particles meets the required minimum N_p_min
  if(SumParticles < N_p_min){
    double RescaleFactor = double(N_p_min)/SumParticles;
    SumParticles = 0;
    for(int i=0; i<N_c; i++){
      N[i] = int(RescaleFactor*N[i]);
      SumParticles += N[i];
    }
  }

  // Pass back total particle number to calling function to allocate
  // memory for PSD.
  *N_tot = SumParticles;

  return 0;
}


void GenPSD(int N_tot, int N_c, int *N, double *y, double *dia, double seed){
  
  int offset;
  RandGenRealUni<double> Dice(seed, 0.0, 1.0);

  offset = 0;
  for(int i=0; i<N_c; i++){
    for(int j=0; j<N[i]; j++){
      y[offset+j] = GenDia(dia[i], dia[i+1], Dice());
    }
    SelectionSort(&y[offset], N[i]); // Sort diameters numerically.
    offset += N[i];
  }

}


void SelectionSort(double *a, int array_size){

  for(int i=0; i<array_size-1; ++i){
    int min;
    double temp;
    
    min = i;
    for (int j=i+1; j<array_size; ++j){
      if(a[j] < a[min]){
        min = j;
      }
    }
    temp = a[i];
    a[i] = a[min];
    a[min] = temp;
  }

}


double CheckPSD(double *psd, int N_tot, double *h, int N_c, double *dia){

	// Calculate accumulated volumes for each class
	double AccumVol[N_c], TotalVol, C_r;

	AccumVol[0] = 0.0;
	for(int i=0, j=0; i<N_tot; i++){
		if(psd[i] > dia[j+1]){
			j++;
			AccumVol[j] = AccumVol[j-1];
		}
		AccumVol[j] += M_PI * psd[i] * psd[i] * psd[i] / 6.0;
	}

	// Normalise accumulated volumes
	TotalVol = AccumVol[N_c-1];
	for(int i=0; i<N_c; i++){
		AccumVol[i] /= TotalVol;
	}

	// Calculate regression coefficient
	C_r = 0.0;
	for(int i=0; i<N_c; i++){
		C_r += (AccumVol[i] - h[i+1]) * (AccumVol[i] - h[i+1]);
	}
	C_r = sqrt(C_r);

	return C_r;
}
