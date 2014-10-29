#ifndef GEN_PSD_H
#define GEN_PSD_H

/*
  Function BetaDist

  Description: BetaDist populates the array containing the classical grading
               curve date according to the expression for the beta distribution
               described in [VOI07]. Eqn. 5.15 p.139 in [RAD11].
*/
void BetaDist(double, double, int, double*);

/*
  Function: ReducedDiam

  Description: ReducedDiam takes the minimum and maximum diameters present in
               the PSD as first and second arguments respectively. The third
               argument is the value of diameter currently considered. The
               function calculates and returns the reduced diameter according
               to eqn. 5.18 p.140 [RAD11].
*/
double ReducedDiam(double*, double*, double*);

#endif // GEN_PSD_H
