#ifndef RANDOMREPAST_H
#define RANDOMREPAST_H
#include "randomizer.h"
#include "repast_hpc/Random.h"
#include "repast_hpc/initialize_random.h"
#define U_MAX 4294967294 // Maximum for unsigned integers = 2^32
//For reasons that are unclear using this wrapper instead of the direct calls to the repast::Random functions slowws things down by about a factor fo two!
class RandomRepast:public randomizer {
public:
    RandomRepast();
    ~RandomRepast();
    void SetSeed( const unsigned );

    void SetSeedFromSystemTime( );

    double GetUniform( );

    unsigned GetUint( );

    double GetNormal( const double& mean = 0.0, const double& standardDeviation = 1.0  );

    double GetLogNormal( const double&, const double& );

};

#endif

