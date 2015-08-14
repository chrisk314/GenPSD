#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <iostream>
#include <ctime>
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
//#define M_PI 3.14159265358979323846264338327
#define M_PI 3.141592653589793
#endif
#define MAX_STR_LEN 50


int main(int argc, char **argv){
  
  char *InputFileName, *OutputFileName;
  FILE *InputFile, *OutputFile;

  // Check correct number of command line args and get file handle.
  if( !(argc == 2) ){
    printf("Invalid number of command line arguments. Program will exit.\n\
            Usage: GenPSD <input file>\n");
    return 1;
  }
  else{
    InputFileName = argv[1];
    InputFile = fopen(InputFileName, "r");
    if(InputFile == NULL){
      printf("Input file %s could not be opened. Program will exit.\n", argv[1]);
      return 1;
    }
  }

  double beta_a, beta_b;
  double d_min, d_max;
  double *dia, *h;
  int N_p_min, N_pc_min, N_tot;
  int N_c, *N;
  double *y, C_rep_target, C_rep;
  bool class_uni;
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
      C_rep_target = strtod(SubStr, &SubStrEnd);
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
    else if( strncasecmp(line, "class_uni", 9) == 0 ){
      strtok(line, " ");
      SubStr = strtok(NULL, "\n");
      class_uni = atoi(SubStr);
    }
  }

  // Display input parameters.
  printf("\n*-----------------------------------------------------\n"\
  	  	  "*\tInput parameters read from file %s:\n"\
		  "*-----------------------------------------------------\n\n"\
		  "beta_a = %lf\n"\
		  "beta_b = %lf\n"\
		  "d_min = %lf\n"\
		  "d_max = %lf\n"\
		  "C_rep = %lf\n"\
		  "N_p_min = %d\n"\
		  "N_pc_min = %d\n"\
		  "N_c = %d\n"\
		  "class_uni = %d\n",\
		  InputFileName, beta_a, beta_b, d_min, d_max, C_rep_target, N_p_min,\
		  	  N_pc_min, N_c, class_uni);

  
  // ------------------------------------------------------------------------------------
  // Generate grading curve, populate classes and generate PSD

  // Allocate memory for analytical grading curve and PSD class populations.
  //h = malloc( (N_c+1) * sizeof(double) );
  //N = malloc( N_c * sizeof(int) );
  dia = (double*)calloc(N_c+1, sizeof(double));
  h = (double*)calloc(N_c+1, sizeof(double));
  N = (int*)calloc(N_c, sizeof(int));

  // Calculate size class boundaries
  dia[0] = d_min;
  dia[N_c] = d_max;
  if(class_uni){
	  // Uniform class width
	  double dia_class_width = (d_max - d_min) / N_c;
	  for(int i=1; i<N_c; i++){
		  dia[i] = dia[i-1] + dia_class_width;
	  }
  }
  else{
	  // Geometric class width
	  double geo_ratio = exp( (log(d_max) - log(d_min)) / N_c );
	  double dia_class_width = ((1 - geo_ratio) / (1 - pow(geo_ratio, N_c))) * (d_max - d_min);
	  for(int i=1; i<N_c; i++){
		  dia[i] = dia[i-1] + dia_class_width * pow(geo_ratio, i-1);
	  }
  }
  printf("\n*-----------------------------------------------------\n"\
		  "Size class boundaries:\n");
  for(int i=0; i<=N_c; i++){
    printf("dia[%d] = %lf\n", i, dia[i]);
  }

  // Generate analytical grading curve with BetaDist.
  BetaDist(beta_a, beta_b, N_c+1, dia, h);
  printf("\nAnalytical grading curve:\n");
  for(int i=0; i<=N_c; i++){
    printf("h[%d] = %lf\n", i, h[i]);
  }

  // Estimate populations from grading curve with constraints.
  GenPopulations(&N_tot, N_c, N_pc_min, N_p_min, N, h, dia);
  printf("\nGrading class populations:\n");
  for(int i=0; i<N_c; i++){
    printf("N[%d] = %d\n", i, N[i]);
  }
  printf("\nTotal number of particles in PSD %d.\n", N_tot);
  
  y = (double*)calloc(N_tot, sizeof(double));

  // Generate particle size distribution.
  GenPSD(N_tot, N_c, N, y, dia, time(NULL));

  // Check statistical degree of representativity between generated and target PSD
  C_rep = CheckPSD(y, N_tot, h, N_c, dia);
  if( C_rep <= C_rep_target ){
	printf("Representativity criteria met with C_r = %1.3e.\n", C_rep);
  }
  else{
	printf("Representativity criteria not met with C_r = %1.3e. Program will exit.\n", C_rep);
	return 1;
  }


  // ------------------------------------------------------------------------------------
  // Write PSD to data file.

  OutputFileName = "psd.dat";
  OutputFile = fopen(OutputFileName,"w");
  if(OutputFile == NULL){
      printf("Output file %s could not be opened. Program will exit.\n", OutputFileName);
      return 1;
  }
  else{
      printf("Writing PSD to file %s.\n", OutputFileName);

      // Get current time for output
	  time_t rawtime;
	  struct tm *timeinfo;
	  char buffer[80];
	  time(&rawtime);
	  timeinfo = localtime(&rawtime);
	  strftime(buffer, 80, "%d-%m-%Y %I:%M:%S", timeinfo);
	  std::string time_str(buffer);

	  fprintf(OutputFile, "# Particle size distribution file %s generated with "\
			  "GenPSD on %s\n", OutputFileName, time_str.c_str());

	  double SphereVol, SumVol=0.0;
	  for(int i=0; i<N_tot; i++){
		SphereVol = M_PI * y[i] * y[i] * y[i] / 6.0;
		SumVol += SphereVol;
		fprintf(OutputFile, "%1.15e %1.15e %1.15e\n", y[i], SphereVol, SumVol);
	  }
  }
  fclose(OutputFile);
  // ------------------------------------------------------------------------------------
  // Write grading to data file.

  OutputFileName = "grading.dat";
  OutputFile = fopen(OutputFileName,"w");
  if(OutputFile == NULL){
      printf("Output file %s could not be opened. Program will exit.\n", OutputFileName);
      return 1;
  }
  else{
      printf("Writing grading to file %s.\n", OutputFileName);

      // Get current time for output
	  time_t rawtime;
	  struct tm *timeinfo;
	  char buffer[80];
	  time(&rawtime);
	  timeinfo = localtime(&rawtime);
	  strftime(buffer, 80, "%d-%m-%Y %I:%M:%S", timeinfo);
	  std::string time_str(buffer);

	  fprintf(OutputFile, "# Particle grading file %s generated with "\
			  "GenPSD on %s\n", OutputFileName, time_str.c_str());

	  double SphereVol, SumVol=0.0;
	  for(int i=0; i<N_c; i++){
		fprintf(OutputFile, "%1.15e %1.15e %1.15e\n", dia[i], dia[i+1], h[i+1]-h[i]);
	  }
  }
  fclose(OutputFile);

  return 0;
}
