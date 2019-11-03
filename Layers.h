#ifndef LAYER_H
#define	LAYER_H

#include "TimeStep.h"
#include "Maths.h"

/**
\file Layers.h
\brief This is the Layers class header file
*/

/**
\brief 2D layers with and without time. 

These classes are used to store 2D driving Layers read in from Netcdf.
Layers are stored with names in a DataLayerSet, from which they are accessed through the Environment class
Layers are added to the DataLayerSet by the Filereader class.
Layers are assumed to be in the same co-ordinate system as the main model grid,
and also gridded, but possibly with a different resolution. 
Grid co-ordinates are stored in longitude and latitude vectors. 
Time is also stored in a vector when present: units currently unspecified!
However an assumption is made that this is monthly data, with one time entry per month.
At present this co-ordinate system is assumed to be Longitude-Latitude.
Data is stored as a 1D c-style floating point array, with the index calculated as
longitudeIndex + latitudeIndex*LongitudeSize + time*LongitudeSize*LatitudeSize.
This is for convenience in reading the data from Netcdf, which supplies the data fields in this form.


@author Mike Bithell <mb425@cam.ac.uk>
*/

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
class Layer{
protected:
    float* _data;
public:
//------------------------------------------------------------------------------------------------------------
/**
 * Data is a pointer to a 2D or 3D array stored as a one-dimensional array,
 * so it needs to be indexed carefully to extract values correctly.
 */
    Layer(float* data):_data(data){}
//------------------------------------------------------------------------------------------------------------
/**
 * Storage for the data is allocated in Filereader.cpp
 */
    virtual ~Layer(){delete[] _data;}
//------------------------------------------------------------------------------------------------------------
    virtual float GetDataAtLonLat(const double,
                                  const double) = 0;      
/**
\brief Find the location of a point in a vector using binary chopping
This is used to return data values at a given point (Longitude,Latitude)
by indexing each co-ordinate separately

@param coord The vector to be searched
@param point Value of the data to be located
@return lo The point should end up between coord[lo] and coord[lo+1]
*/
    unsigned binaryChop(std::vector<float> coord, 
                        float point){
        unsigned lo=0,hi=coord.size()-1;
        float sgn=1;
        if (coord[hi]<coord[lo]){sgn=-1;}
        while (hi > lo+1){
         if (sgn*point<sgn*coord[(hi+lo)/2]) hi=(hi+lo)/2; else lo=(hi+lo)/2;
        }
        return lo;
    }
};
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
class Layer2D:public Layer{
protected:
    std::vector<float>_lon,_lat;
public:
//------------------------------------------------------------------------------------------------------------
    Layer2D(float* data, 
            std::vector<float>lon,
            std::vector<float>lat):Layer(data),_lon(lon),_lat(lat){}
//------------------------------------------------------------------------------------------------------------
    float GetDataAtLonLat(const double Longitude,
                          const double Latitude){
        
        unsigned lonix=binaryChop(_lon,Longitude);
        unsigned latix=binaryChop(_lat,Latitude);
        unsigned cellIndex=lonix+latix*_lon.size();

        return _data[cellIndex];}
};
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
class Layer2DWithTime:public Layer2D{
    std::vector<float>_time;
public:
//------------------------------------------------------------------------------------------------------------
    Layer2DWithTime(float* data,std::vector<float>lon,std::vector<float>lat,std::vector<float>time):Layer2D(data,lon,lat){}
//------------------------------------------------------------------------------------------------------------
    float GetDataAtLonLat(const double Longitude,
                          const double Latitude){
        unsigned lonix=binaryChop(_lon,Longitude);
        unsigned latix=binaryChop(_lat,Latitude);
        unsigned cellIndex=lonix+latix*_lon.size();
        // FIX: The way the month index is calculated here is due to the use of an annual climatology - also assumes time units of month in the data
        unsigned monthIndex = TimeStep::Get( )->Get( Constants::cMonthlyTimeUnitName ) % 12 ;
        unsigned dataIndex = cellIndex + ( monthIndex * _lon.size() * _lat.size() );
        return _data[dataIndex];
        
    }
};
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

#endif
