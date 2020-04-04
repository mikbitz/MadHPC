 /* 
 * UtilityFunctions.h
 * Created on: May 7, 2018
 *      Author: Mike Bithell
 * 
 *  Derived from Origianl C# code by
 * Emergent Global Patterns of Ecosystem Structure and Function from a Mechanistic General Ecosystem Model , Michael B. J. Harfoot , Tim Newbold, Derek P. Tittensor,  Stephen Emmott, Jon Hutton, Vassily Lyutsarev, Matthew J. Smith, Jörn P. W. Scharlemann, Drew W. Purves PLOS Biology 12(4) 2014 e1001841 https://doi.org/10.1371/journal.pbio.1001841
 */
#ifndef UTILITYFUNCTIONS
#define UTILITYFUNCTIONS

#include "Parameters.h"
#include "Convertor.h"

#include <math.h>
#include <algorithm>
#include <vector>
#include <map>
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock
#include <assert.h>
#include <iostream>

using namespace std;

class UtilityFunctions {
public:

    /** \brief If longitudinal cell coordinates run from 0 to 360, the convert to -180 to 180 values
    @param lons The longitudinal coorindates of the cells in the model grid */
    void ConvertToM180To180( std::vector<double>& lons );


    /** \briefConverts values per square km to per square degree, given cell latitude
    @param valueToConvert The value per square km
    @param latitude The latitude of the grid cell
    @return The specified value converted to per square degree */
    double ConvertSqMToSqDegrees( double valueToConvert, double latitude );

    /** \brief Calculates the probability of a particular value under a log-normal distribution with specified mean and standard deviation
    @param xValue The value to return the probability of under the log-normal distribtuion, in identity space
    @param meanIdentity The mean of the log-normal distribution, in identity space
    @param standardDeviation The standard deviation of the log-normal distribution, in log space
    @return The probability of the specified value under the specified log-normal distribution
     */
    double LogNormalPDF( double xValue, double meanIdentity, double standardDeviation );

    /** \brief Calculates the probability of a particular value under a normal distribution with specified mean and standard deviation
    @param xValue The value to return the probability of under the normal distribtuion
    @param meanValue The mean of the normal distribution
    @param standardDeviation The standard deviation of the normal distribution
    @return The probability of the specified value under the specified normal distribution */
    double NormalPDF( double xValue, double meanValue, double standardDeviation );

    /** \brief Calculate the area of a grid cell in square km, given its dimensions and geographical position
    @param latitude The latitude of the bottom-left corner of the grid cell
    @param lonCellSize The longitudinal dimension of the grid cell
    @param latCellSize The latitudinal dimension of the grid cell
    @return The area in square km of the grid cell */
    double CalculateGridCellArea( double latitude, double cellSize );

    /** \brief Calculate the length of a degree of latitude at a particular latitude /
    @param latitude The latitude of the bottom-left corner of the grid cell MB ?? This is not consistent with cellarea in the original code
    @return The length of a degree of latitude in kilometres */
    double CalculateLengthOfDegreeLatitude( float latitude );

    /** \brief Calculate the length of a degree of longitude at a particular latitude
    @param latitude The latitude of the bottom-left corner of the grid cell
    @return The length of a degree of longitude in kilometres MB ?? This is not consistent with cellarea in the original code */
    double CalculateLengthOfDegreeLongitude( float latitude );
    
    /** \brief Calculate the approximate distance between two points on a sphere but report result in degrees
    @param lon1 The longitude of the first point
    @param lat1 The latitude of the first point
    @param lon2 The longitude of the second point
    @param lat2 The latitude of the second point
    @return The distance between two points in kilometres */

    double HaversineDistanceInDegrees(double lon1,double lat1,double lon2,double lat2);
    /** \brief Calculate the approximate distance between two points on a sphere
    @param lon1 The longitude of the first point
    @param lat1 The latitude of the first point
    @param lon2 The longitude of the second point
    @param lat2 The latitude of the second point
    @return The distance between two points in kilometres */
    double HaversineDistance(double lon1,double lat1,double lon2,double lat2);
    
    /** \brief Convert from degrees to radians
    @param degrees The value in degrees to convert
    @return The value converted to radians</returns> */
    double DegreesToRadians( double degrees );
};
#endif
