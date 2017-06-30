#include "UtilityFunctions.h"

void UtilityFunctions::ConvertToM180To180( std::vector<double>& lons ) {
    // Loop over longitudinal coordinates of the model grid cells
    for( unsigned jj = 0; jj < lons.size( ); jj++ ) {
        // If longitudinal coorindates exceed 180, then subtrarct 360 to correct the coorindates
        if( lons[jj] >= 180.0 ) {
            lons[jj] -= 360.0;
        }
    }
    // Re-sort the longitudinal coordinates
    sort( lons.begin( ), lons.end( ) );
}


double UtilityFunctions::ConvertSqMToSqDegrees( double valueToConvert, double latitude ) {
    // Convert the value to per square degree using the cosine of latitude and assuming cell dimensions of 110km by 110km at the Equator
    return valueToConvert * 110000.0 * 110000.0 * cos( DegreesToRadians( latitude ) );
}

double UtilityFunctions::CalculateGridCellArea( double latitude, double cellSize ) {
    const double PI = acos( -1. );
    // Convert from degrees to radians
    double latitudeRad = DegreesToRadians( latitude );

    // Equatorial radius in metres
    double EquatorialRadius = 6378137;

    // Polar radius in metres
    double PolarRadius = 6356752.3142;

    // Angular eccentricity
    double AngularEccentricity = acos( DegreesToRadians( PolarRadius / EquatorialRadius ) );

    // First eccentricity squared
    //double ESquared = pow(sin(DegreesToRadians(AngularEccentricity)), 2);

    // Flattening
    //double Flattening = 1 - cos(DegreesToRadians(AngularEccentricity));

    // Temporary value to save computations
    double TempVal = pow( ( EquatorialRadius * cos( latitudeRad ) ), 2 ) + pow( ( PolarRadius * sin( latitudeRad ) ), 2 );

    // Meridional radius of curvature
    double MPhi = pow( EquatorialRadius * PolarRadius, 2 ) / pow( TempVal, 1.5 );

    // Normal radius of curvature
    double NPhi = pow( EquatorialRadius, 2 ) / sqrt( TempVal );

    // Length of latitude (km)
    double LatitudeLength = PI / 180 * MPhi / 1000;

    // Length of longitude (km)
    double LongitudeLength = PI / 180 * cos( latitudeRad ) * NPhi / 1000;

    // Return the cell area in km^2
    return LatitudeLength * cellSize * LongitudeLength * cellSize;
}

double UtilityFunctions::CalculateLengthOfDegreeLatitude( float latitude ) {
    const double PI = acos( -1. );
    // Convert from degrees to radians
    double latitudeRad = DegreesToRadians( latitude );

    // Equatorial radius in metres
    double EquatorialRadius = 6378137;

    // Polar radius in metres
    double PolarRadius = 6356752.3142;

    // Angular eccentricity
    //double AngularEccentricity = acos(DegreesToRadians(PolarRadius / EquatorialRadius));

    // First eccentricity squared
    //double ESquared = pow(sin(DegreesToRadians(AngularEccentricity)), 2);

    // Flattening
    //double Flattening = 1 - cos(DegreesToRadians(AngularEccentricity));

    // Temporary value to save computations
    double TempVal = pow( ( EquatorialRadius * cos( latitudeRad ) ), 2 ) + pow( ( PolarRadius * sin( latitudeRad ) ), 2 );

    // Meridional radius of curvature
    double MPhi = pow( EquatorialRadius * PolarRadius, 2 ) / pow( TempVal, 1.5 );

    // Length of latitude (km)
    return PI / 180 * MPhi / 1000;
}

double UtilityFunctions::CalculateLengthOfDegreeLongitude( float latitude ) {
    const double PI = acos( -1. );
    // Convert from degrees to radians
    double latitudeRad = DegreesToRadians( latitude );

    // Equatorial radius in metres
    double EquatorialRadius = 6378137;

    // Polar radius in metres
    double PolarRadius = 6356752.3142;

    // Temporary value to save computations
    double TempVal = pow( ( EquatorialRadius * cos( latitudeRad ) ), 2 ) + pow( ( PolarRadius * sin( latitudeRad ) ), 2 );

    // Normal radius of curvature
    double NPhi = pow( EquatorialRadius, 2 ) / sqrt( TempVal );

    // Length of longitude (km)
    return PI / 180 * cos( latitudeRad ) * NPhi / 1000;
}

double UtilityFunctions::DegreesToRadians( double degrees ) {
    const double PI = acos( -1. );
    return (degrees * PI / 180.0 );
}
