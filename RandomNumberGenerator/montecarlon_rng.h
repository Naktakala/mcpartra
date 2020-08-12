#ifndef _montecarlon_rng_h
#define _montecarlon_rng_h

#include "../chi_montecarlon.h"
//#define R123_UNIFORM_FLOAT_STORE 1
//#include<Random123/threefry.h>
//#include<Random123/uniform.hpp>

#include <random>

//#########################################################
/**Random number generator based on threefry.*/
class chi_montecarlon::RandomNumberGenerator
{
private:
//  r123::Threefry2x64 generator;
//  r123::Threefry2x64::ctr_type  ctr;
//  r123::Threefry2x64::key_type key;
//
//  bool    flipFlop;
//  double  storedNumber;

  std::mt19937_64    mt1993764_generator;
  std::uniform_real_distribution<double> distribution;

public:
  RandomNumberGenerator() :
    distribution(0.0,1.0)
  {
    mt1993764_generator.seed(0);
  }
  RandomNumberGenerator(int seed) :
    distribution(0.0,1.0)
  {
    mt1993764_generator.seed(seed);
  }
  double Rand()
  {
    return distribution(mt1993764_generator);
  }
//  //=================================== Default constructor
//  /**Default constructor. Seed=0, Key=0*/
//  RandomNumberGenerator()
//  {
//    ctr = {{0,0}};
//    key = {{0,0}};
//
//    flipFlop = false;
//    storedNumber=0.5;
//  }
//
//  /**Constructor for a given seed and key.*/
//  RandomNumberGenerator(uint64_t seed, uint64_t stream)
//  {
//    ctr = {{seed,seed}};
//    key = {{stream,0}};
//
//    flipFlop = false;
//    storedNumber=0.5;
//  }
//
//  /**Returns a random number from the counter*/
//  double Rand()
//  {
//    if (flipFlop==true)
//    {
//      flipFlop=false;
//      return storedNumber;
//    }
//    else
//    {
//      flipFlop=true;
//      ctr[0]++;
//      ctr[1]++;
//      r123::Threefry2x64::ctr_type rand = generator(ctr, key);
//      storedNumber = r123::u01<double>(rand.v[0]);
//      return         r123::u01<double>(rand.v[1]);
//    }
//
//  }
};


#endif
