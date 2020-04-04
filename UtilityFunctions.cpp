 /* 
 * UtilityFunctions.cpp
 * Created on: May 7, 2018
 *      Author: Mike Bithell
 * 
 *  Derived from Origianl C# code by
 * Emergent Global Patterns of Ecosystem Structure and Function from a Mechanistic General Ecosystem Model , Michael B. J. Harfoot , Tim Newbold, Derek P. Tittensor,  Stephen Emmott, Jon Hutton, Vassily Lyutsarev, Matthew J. Smith, Jörn P. W. Scharlemann, Drew W. Purves PLOS Biology 12(4) 2014 e1001841 https://doi.org/10.1371/journal.pbio.1001841
 */
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
double UtilityFunctions::HaversineDistanceInDegrees(double lon1,double lat1,double lon2,double lat2){
    //inputs in degrees.
    //approximate distance between two points
    //more accurate version would use something like Vincenty's method
    double phi1=DegreesToRadians( lat1 );
    double phi2=DegreesToRadians( lat2 );
    double lam1=DegreesToRadians( lon1 );
    double lam2=DegreesToRadians( lon2 );
    double r=180./acos(-1.); //180 over pi to convert radians to degrees
    //use min function to ensure sqrt argument doesn't try to exceed 1
    return 2 * r * asin( sqrt( min( 1., ( pow( sin( (phi2-phi1)/2 ), 2 ) + cos(phi1) * cos(phi2) * pow( sin( (lam2-lam1)/2 ),2 ) ) ) ) );
}
double UtilityFunctions::HaversineDistance(double lon1,double lat1,double lon2,double lat2){
    //inputs in degrees.
    //approximate distance between two points
    //more accurate version would use something like Vincenty's method
    double phi1=DegreesToRadians( lat1 );
    double phi2=DegreesToRadians( lat2 );
    double lam1=DegreesToRadians( lon1 );
    double lam2=DegreesToRadians( lon2 );
    double r=(6378137+6356752.3142)/2./1000.; //value in km, assuming a sphere
    //use min function to ensure sqrt argument doesn't try to exceed 1
    return 2 * r * asin( sqrt( min( 1., ( pow( sin( (phi2-phi1)/2 ), 2 ) + cos(phi1) * cos(phi2) * pow( sin( (lam2-lam1)/2 ),2 ) ) ) ) );
}
double UtilityFunctions::DegreesToRadians( double degrees ) {
    const double PI = acos( -1. );
    return (degrees * PI / 180.0 );
}
