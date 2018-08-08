#include "repastRandom.h"

repastRandom::repastRandom( ) {
    // These values are not magical, just the default values Marsaglia used.
    // Any pair of unsigned integers should be fine.
    m_w = 521288629;
    m_z = 362436069;

    PI = acos( -1. );
    Temp2PI = 2.0 * PI;
}

void repastRandom::Reset( ) {
    // These values are not magical, just the default values Marsaglia used.
    // Any pair of unsigned integers should be fine.
    m_w = 521288629;
    m_z = 362436069;
}

void repastRandom::SetSeed( unsigned u, unsigned v ) {
    if( u != 0 ) m_w = u;
    if( v != 0 ) m_z = v;
}

void repastRandom::SetSeed( unsigned u ) {
    m_w = u;
}

void repastRandom::SetSeedFromSystemTime( ) {
    //System.DateTime dt = System.DateTime.Now;
    //long x = dt.ToFileTime();
    //SetSeed((uint)(x >> 16), (uint)(x % 4294967296));
    unsigned seed = std::chrono::system_clock::now( ).time_since_epoch( ).count( );

    SetSeed( seed >> 16, seed % 4294967296 );
}

double repastRandom::GetUniform( ) {
    // 0 <= u < 2^32
    unsigned u = GetUint( );
    // The magic number below is 1/(2^32 + 2).
    // The result is strictly between 0 and 1.
    return ( u + 1.0 ) * 2.328306435454494e-10;
}

unsigned repastRandom::GetUint( ) {
    m_z = 36969 * ( m_z & 65535 ) + ( m_z >> 16 );
    m_w = 18000 * ( m_w & 65535 ) + ( m_w >> 16 );
    return (m_z << 16 ) +m_w;
}

double repastRandom::GetNormal( ) {
    // Use Box-Muller algorithm
    double u1 = GetUniform( );
    double u2 = GetUniform( );
    double r = sqrt( -2.0 * log( u1 ) );
    //double theta = 2.0 * Math.PI * u2;
    //return r * Math.Sin(theta);
    return r * sin( Temp2PI * u2 );
}

double repastRandom::GetNormal( double mean, double standardDeviation ) {
    //if (standardDeviation <= 0.0)
    //{
    //string msg = string.Format("Shape must be positive. Received {0}.", standardDeviation);
    //throw new ArgumentOutOfRangeException(msg);
    assert( standardDeviation > 0. );
    //}
    return mean + standardDeviation * GetNormal( );
}


double repastRandom::GetLogNormal( double mu, double sigma ) {
    return exp( GetNormal( mu, sigma ) );
}


