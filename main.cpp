#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include "GenPSD.h"
#include "integrate.h"
//#ifdef __cplusplus
//extern "C"
//{
//#endif
//#include "integrate.h"
//#ifdef __cplusplus
//}
//#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif
#define MAX_STR_LEN 50


int main(int argc, char **argv){
  
  char *InputFileName;
  FILE *InputFile, *OutputFile;

  // Check correct number of command line args and get file handle.
  if( !(argc == 2) ){
    printf("Invalid number of command line arguments. Program will exit.\n\
            Usage: GenPSD input_file\n");
    exit(-1);
  }
  else{
    InputFileName = argv[1];
    InputFile = fopen(InputFileName, "r");
    if(InputFile == NULL){
      printf("Input file %s could not be opened. Program will exit.\n", argv[1]);
      exit(-1);
    }
  }

  double beta_a, beta_b;
  double d_min, d_max;
  double C_rep, *h;
  int N_p_min, N_pc_min, N_tot;
  int N_c, *N;
  char line[MAX_STR_LEN];

  // Read parameters from file
  while( fgets(line, MAX_STR_LEN, InputFile) != NULL ){
    char *SubStr, *SubStrEnd;
    if( strncasecmp(line, "beta_a", 6) == 0 ){
      strtok(line, " ");
      SubStr = strtok(NULL, "\n");
      beta_a = strtod(SubStr, &SubStrEnd);
    }
    else if( strncasecmp(line, "beta_b", 6) == 0 ){
      strtok(line, " ");
      SubStr = strtok(NULL, "\n");
      beta_b = strtod(SubStr, &SubStrEnd);
    }
    else if( strncasecmp(line, "d_min", 5) == 0 ){
      strtok(line, " ");
      SubStr = strtok(NULL, "\n");
      d_min = strtod(SubStr, &SubStrEnd);
    }
    else if( strncasecmp(line, "d_max", 5) == 0 ){
      strtok(line, " ");
      SubStr = strtok(NULL, "\n");
      d_max = strtod(SubStr, &SubStrEnd);
    }
    else if( strncasecmp(line, "C_rep", 5) == 0 ){
      strtok(line, " ");
      SubStr = strtok(NULL, "\n");
      C_rep = strtod(SubStr, &SubStrEnd);
    }
    else if( strncasecmp(line, "N_p_min", 7) == 0 ){
      strtok(line, " ");
      SubStr = strtok(NULL, "\n");
      N_p_min = atoi(SubStr);
    }
    else if( strncasecmp(line, "N_pc_min", 8) == 0 ){
      strtok(line, " ");
      SubStr = strtok(NULL, "\n");
      N_pc_min = atoi(SubStr);
    }
    else if( strncasecmp(line, "N_c", 3) == 0 ){
      strtok(line, " ");
      SubStr = strtok(NULL, "\n");
      N_c = atoi(SubStr);
    }
  }

  // Display input parameters.
  printf("\n!!-----------------------------------------------------!!\n");
  printf("Input parameters read from file %s:\n"\
        "!!-----------------------------------------------------!!\n"\
        "beta_a = %lf\n"\
        "beta_b = %lf\n"\
        "d_min = %lf\n"\
        "d_max = %lf\n"\
        "C_rep = %lf\n"\
        "N_p_min = %d\n"\
        "N_pc_min = %d\n"\
        "N_c = %d\n\n",\
        InputFileName, beta_a, beta_b, d_min, d_max, C_rep, N_p_min, N_pc_min, N_c);
  
  // Allocate memory for analytical grading curve and PSD class populations.
  //h = malloc( (N_c+1) * sizeof(double) );
  //N = malloc( N_c * sizeof(int) );
  h = (double*)calloc(N_c+1, sizeof(double));
  N = (int*)calloc(N_c, sizeof(int));

  // Generate analytical grading curve with BetaDist.  
  BetaDist(beta_a, beta_b, N_c+1, h);
  for(int i=0; i<=N_c; i++){
    printf("h[%d] = %lf\n", i, h[i]);
  }

  // Estimate populations from grading curve with constraints.
  GenPopulations(&N_tot, N_c, N_pc_min, N_p_min, N, h, d_min, d_max);
  for(int i=0; i<N_c; i++){
    printf("N[%d] = %d\n", i, N[i]);
  }

  //double IntResult;
  //if( Integrate(&IntResult, -4, 0, 10, 26, 1E-10, MyFunc) == 1 ){
  //  printf("Function Integrate did not converge. Program will exit.\n");
  //  exit(-1);
  //}
  //printf("IntResult = %+1.15e\n", IntResult);

  return 0;

}
