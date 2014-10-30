#include "GenPSD.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/generator_iterator.hpp>

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


int BetaDist(double beta_a, double beta_b, int NumElem, double *h){

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
    x = i*IntStop/(NumElem-1);
    if( Integrate(&h[i], IntStart, x, IntIter, IntSubDiv, IntErr,\
            MyFunc) != 0 )
    {
      printf("Integral did not converge in function Integrate called from "\
                "function BetaDist. Program will exit.\n");
      return 1;
    }
    h[i] *= B_factor; 
  }
  
  MyFunc->~Func();

  return 0;
}


double GenDia(double d_lower, double d_upper, double y){
    return (d_lower*d_upper) / (d_upper + y*(d_lower - d_upper));
}

int GenPopulations(int *N_tot, int N_c, int N_pc_min, int N_p_min, int *N,\
                    double *h, double d_min, double d_max){
    
  double DiaDelta;                      // Width of each size class interval.
  double MeanDia;
  double VolMeanDia;                    // Mean volume for given size class.
  double TotalVolEst;
  double SpherePrefactor;               
  int N_min=VERY_BIG_INT, SumParticles;

  DiaDelta = (d_max - d_min) / N_c;
  SpherePrefactor = (1.0/6.0) * M_PI;

  // Estimate total volume.
  MeanDia = GenDia(d_min, d_max, 0.5);
  VolMeanDia = SpherePrefactor * MeanDia * MeanDia * MeanDia;  //
  TotalVolEst = VolMeanDia * N_p_min;

  // Get first estimate for size class populations from grading curve h(d).
  SumParticles = 0;
  for(int i=0; i<N_c; i++){
    MeanDia = GenDia(d_min + i*DiaDelta, d_min + (i+1)*DiaDelta, 0.5);
    VolMeanDia = SpherePrefactor * MeanDia * MeanDia * MeanDia;
    N[i] = (h[i+1] - h[i]) * ( TotalVolEst / VolMeanDia ); // Using volume esimate.
    //N[i] = (h[i+1] - h[i]) / VolMeanDia;  // Assuming unit volume.
    printf("N[%d] = %d\n", i, N[i]);
    SumParticles += N[i];
    if(N[i] < N_min){
      N_min = N[i];
    }
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
      N[i] = int(RescaleFactor*N[i]);
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


void GenPSD(int N_tot, int N_c, int *N, double *y, double d_min, double d_max,\
                double seed){
    
  boost::mt19937 MersenneTwister(seed);
  boost::uniform_real<float> uniform_range( 0.0, 1.0 );
  boost::variate_generator< boost::mt19937, boost::uniform_real<float> >\
    Dice(MersenneTwister, uniform_range);
  
  double DiaDelta, DiaStart, DiaStop;

  DiaDelta = (d_max - d_min) / N_c; // Class diameter width.
  
  for(int i=0; i<N_c; i++){
    DiaStart = d_min + i*DiaDelta;
    DiaStop = DiaStart + DiaDelta;
    int offset = 0; // Track first element in y[] for current size class.
    for(int j=0; j<N[i]; j++){
      y[offset+j] = GenDia(DiaStart, DiaStop, Dice()); 
    }
    SelectionSort(&y[offset], N[i]); // Sort diameters numerically.
    offset += N[i];
  }

}


double ReducedDiam(double *d_min, double *d_max, double *d){
  double d_r = (d - d_min)/(d_max - d_min);
  return d_r;
}


void SelectionSort(double *a, int array_size){

  for(int i=0; i<array_size-1; ++i){
     int min, temp;
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
