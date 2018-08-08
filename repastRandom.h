
#ifndef REPASTRANDOM_H
#define	REPASTRANDOM_H

#include <math.h>
#include <assert.h>
#include <chrono>

/// <summary>
/// SimpleRNG is a simple random number generator based on 
/// George Marsaglia's MWC (multiply with carry) generator.
/// Although it is very simple, it passes Marsaglia's DIEHARD
/// series of random number generator tests.
/// 
/// Written by John D. Cook 
/// http://www.johndcook.com
/// </summary>

class repastRandom {
public:
    /// <summary>
    /// Constructor the random number generator: sets default values for the integers for George Marsaglia's MWC algorithm
    /// </summary>
    repastRandom( );
    ~repastRandom(){;}
    void Reset( );
    // The random generator seed can be set three ways:
    // 1) specifying two non-zero unsigned integers
    // 2) specifying one non-zero unsigned integer and taking a default value for the second
    // 3) setting the seed from the system time

    /// <summary>
    /// Set the seed of the random number generator using the specified integers
    /// </summary>
    /// <param name="u">An integer</param>
    /// <param name="v">An integer</param>
    void SetSeed( unsigned u, unsigned v );

    /// <summary>
    /// Set the seed of the random number generator using the specified integer
    /// </summary>
    /// <param name="u">An integer</param>
    void SetSeed( unsigned u );

    /// <summary>
    /// Set the seed of the random number generator based on the system time
    /// </summary>

    void SetSeedFromSystemTime( );

    /// <summary>
    /// A random draw from a uniform distribution between 0 and 1
    /// </summary>
    /// <returns>A random draw from a uniform distribution between 0 and 1</returns>
    /// <remarks>This will not return either 0 or 1</remarks>
    double GetUniform( );

    /// <summary>
    /// Get a random unsigned integer using uses George Marsaglia's MWC algorithm
    /// </summary>
    /// <returns>A random unsigned integer using uses George Marsaglia's MWC algorithm</returns>
    /// <remarks>See http://www.bobwheeler.com/statistics/Password/MarsagliaPost.txt </remarks>
    unsigned GetUint( );

    /// <summary>
    /// A random draw from a normal distribution with mean 0 and standard deviation 1
    /// </summary>
    /// <returns>A random draw from a normal distribution with mean 0 and standard deviation 1</returns>
    double GetNormal( );

    /// <summary>
    /// A random draw from a normal distribution
    /// </summary>
    /// <param name="mean">The mean of the normal distribution</param>
    /// <param name="standardDeviation">The standard deviation of the normal distribution</param>
    /// <returns>A random draw from a normal distribution</returns>
    double GetNormal( double mean, double standardDeviation );


    /// <summary>
    /// A random draw from a lognormal distribution
    /// </summary>
    /// <param name="mu">Mean of the lognormal distribution</param>
    /// <param name="sigma">Standard deviation of the lognormal distribution</param>
    /// <returns>A random draw from a lognormal distribution</returns>
    double GetLogNormal( double mu, double sigma );


private:
    /// <summary>
    /// Integer for George Marsaglia's MWC algorithm 
    /// </summary>
    unsigned m_w;
    /// <summary>
    /// Integer for George Marsaglia's MWC algorithm
    /// </summary>
    unsigned m_z;
    double PI;
    double Temp2PI;
};
#endif
