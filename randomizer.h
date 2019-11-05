
#ifndef RANDOMIZER_H
#define	RANDOMIZER_H

#include <math.h>
#include <assert.h>
#include <chrono>

#define U_MAX 4294967294 // Maximum for unsigned integers = 2^32
/// <summary>
/// Interface for the random number generator - allows for switching between different versions.
/// </summary>

class randomizer {
public:
    /// <summary>
    /// Constructor the random number generator: sets default values
    /// </summary>
    randomizer( ){};
    virtual ~randomizer(){;}
    
    // The random generator seed can be set three ways:
    // 1) specifying one non-zero unsigned integer
    // 2) setting the seed from the system time


    /// <summary>
    /// Set the seed of the random number generator using the specified integer
    /// </summary>
    /// <param name="u">An integer</param>
    virtual void SetSeed(const unsigned  )=0;

    /// <summary>
    /// Set the seed of the random number generator based on the system time
    /// </summary>

    virtual void SetSeedFromSystemTime( )=0;

    /// <summary>
    /// A random draw from a uniform distribution between 0 and 1
    /// </summary>
    /// <returns>A random draw from a uniform distribution between 0 and 1</returns>
    /// <remarks>This will not return either 0 or 1</remarks>
    virtual double GetUniform( )=0;

    /// <summary>
    /// Get a random unsigned integer
    virtual unsigned GetUint( )=0;

    /// <summary>
    /// A random draw from a normal distribution with mean 0 and standard deviation 1
    /// </summary>
    /// <returns>A random draw from a normal distribution with default mean 0 and standard deviation 1</returns>
    virtual double GetNormal(const double& mean = 0.0, const double& standardDeviation = 1.0  )=0;


    /// <summary>
    /// A random draw from a lognormal distribution
    /// </summary>
    /// <param name="mu">Mean of the lognormal distribution</param>
    /// <param name="sigma">Standard deviation of the lognormal distribution</param>
    /// <returns>A random draw from a lognormal distribution</returns>
    virtual double GetLogNormal( const double&, const double& )=0;


};
#endif


   


