#ifndef RAND_GEN_H
#define RAND_GEN_H

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/generator_iterator.hpp>

template <class RealT>
class RandGenRealUni{

 public:
  boost::mt19937 Generator; // Mersenne twister.
  boost::uniform_real<RealT> Distribution; // Uniform distribution.
  boost::variate_generator< boost::mt19937, boost::uniform_real<RealT> > Dice;

  RandGenRealUni(double Seed, double RangeMin, double RangeMax):
    Generator(Seed), Distribution(RangeMin, RangeMax),
    Dice(Generator, Distribution){}; 

  RealT operator() (void){
    return Dice();
  }

};

#endif // RAND_GEN_H
