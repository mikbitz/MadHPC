#ifndef LAYER_H
#define	LAYER_H

#include "TimeStep.h"
#include "Maths.h"

/**
\file Layers.h
\brief This is the Layers class header file. Current classes appear only in this header
*/

/**
\brief 2D layers with and without time. 

These classes are used to store 2D driving Layers read in from Netcdf plus possibly some fields calcuated from them.
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
/**
 * \brief Storage for the dataset
 */
    float* _data;
public:
//------------------------------------------------------------------------------------------------------------
/**
 * \brief Base-class constructor
 * 
 * The base class has just a single variable to hold data.
 * Data is a pointer to a 2D or 3D array stored as a one-dimensional array,
 * so it needs to be indexed carefully to extract values correctly.
 * 
 * @param data A pointer to a float array
 */
    Layer(float* data):_data(data){}
//------------------------------------------------------------------------------------------------------------
/**
 * \brief Destructor - removes the storage in _data
 * 
 * Storage for the data is allocated in Filereader.cpp
 */
    virtual ~Layer(){delete[] _data;}
//------------------------------------------------------------------------------------------------------------
/**
* \brief Function to return data at a point
*
* This is used to return data values at a given point (Longitude,Latitude)
* Since the base class has no co-ordinates, this is pure virtual function at this point

@param Longitude x-co-ordinate
@param Latitude  y-co-ordinate
@return _data at the required point in space
*/
    virtual float GetDataAtLonLat(const double,
                                  const double) = 0;      
/**
* \brief Find the location of a point in an (ordered) vector using binary chopping
*
* This is used to return data values at a given point (Longitude,Latitude)
* by indexing each co-ordinate separately. NB co-ordinates must be monotonic!

@param coord The vector to be searched
@param point Value of the data to be located
@return lo The point should end up between coord[lo] and coord[lo+1]
*/
unsigned binaryChop(std::vector<float> coord, 
                        float point){
        unsigned lo=0,hi=coord.size()-1;
        //need to make sure works for vector in increasing or decreasing order
        //so reverse signs if the data is decreasing
        float sgn=1;
        if (coord[hi]<coord[lo]){sgn=-1;}
        while (hi > lo+1){
         if (sgn*point<sgn*coord[(hi+lo)/2]) hi=(hi+lo)/2; else lo=(hi+lo)/2;
        }
        if (hi==lo+1 && sgn*point>=sgn*coord[hi])lo=hi;
        return lo;
}
virtual double bilin(const double lon, const double lat,const unsigned lonix,const unsigned latix,const unsigned tix=0)=0;
virtual float& operator()(unsigned,unsigned,unsigned tix=0)=0;
};
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
class Layer2D:public Layer{
/**
 * \brief A two dimnesional spatial Layer
 * 
 * Assumed to be indexed by two variables, nominally Longitude and Latitude, although
 * in principle could be any pair of doubles (e.g. metres in UTM). The sizes of these two
 * vectors need to be correctly commensurate with the size and shape of _data
 * Inherits storage _data from Layer base class
 */
protected:
    std::vector<float>_lon;
    std::vector<float>_lat;
public:
//------------------------------------------------------------------------------------------------------------
Layer2D(float* data, 
        std::vector<float>lon,
        std::vector<float>lat):Layer(data),_lon(lon),_lat(lat){}
//------------------------------------------------------------------------------------------------------------
float& operator()(unsigned lonix,unsigned latix,unsigned tix=0){
         assert (lonix<_lon.size());
         assert (latix<_lat.size());
         return _data[lonix+latix*_lon.size()];
}
//------------------------------------------------------------------------------------------------------------
double bilin(const double lon, const double lat,const unsigned lonix,const unsigned latix,const unsigned tix=0){
    unsigned lal=latix,lau=latix+1;
    if (lau>=_lat.size())lau=0;
    unsigned lol=lonix,lou=lonix+1;
    if (lou>=_lon.size())lou=_lon.size()-1;
    
    double data_ll=(*this)(lol,lal);
    double data_lu=(*this)(lol,lau);
    double data_ul=(*this)(lou,lal);
    double data_uu=(*this)(lou,lau);
    double data_av=0,dn=0;
    if (data_ll!=Constants::cMissingValue){data_av+=data_ll;dn+=1;}  
    if (data_lu!=Constants::cMissingValue){data_av+=data_lu;dn+=1;}  
    if (data_ul!=Constants::cMissingValue){data_av+=data_ul;dn+=1;}  
    if (data_uu!=Constants::cMissingValue){data_av+=data_uu;dn+=1;}  
    if (dn!=0)data_av=data_av/dn;else data_av=Constants::cMissingValue;
    if (data_ll==Constants::cMissingValue){data_ll=data_av;}  
    if (data_lu==Constants::cMissingValue){data_lu=data_av;}  
    if (data_ul==Constants::cMissingValue){data_ul=data_av;}  
    if (data_uu==Constants::cMissingValue){data_uu=data_av;}
    
    assert(data_ll != Constants::cMissingValue && data_lu != Constants::cMissingValue && data_ul != Constants::cMissingValue && data_ll != Constants::cMissingValue );
    double result=
     data_ll*(_lon[lou]-lon)*(_lat[lau]-lat)
    -data_ul*(_lon[lol]-lon)*(_lat[lau]-lat)
    -data_lu*(_lon[lou]-lon)*(_lat[lal]-lat)
    +data_uu*(_lon[lol]-lon)*(_lat[lal]-lat);
    double dlon=_lon[lou]-_lon[lol],dlat=_lat[lau]-_lat[lal];
    if (dlon==0)dlon=1;
    if (dlat==0)dlat=1;
    result=result/dlon/dlat;
    return result;
}
//------------------------------------------------------------------------------------------------------------
float GetDataAtLonLat(const double Longitude,
                      const double Latitude){
        
        unsigned lonix=binaryChop(_lon,Longitude);
        unsigned latix=binaryChop(_lat,Latitude);

        return (*this)(lonix,latix);
        //return bilin(Longitude,Latitude,lonix,latix);
}
};
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
class Layer2DWithTime:public Layer2D{
    /**
     * \brief A two dimensional spatial Layer with a time-varying component
     * 
     * Spatially indexed as for layer2D, but now with a time variable - it is assumed the current timevariable
     * is extracted from the TimeStep class when needed, so that the interface is the same as for Layer2D.
     */
    /** \brief A vector to hold the times (in days) corresponding to the data points */
    std::vector<float>_time;
    /** \brief Number of days between the last day in the time array in a given year and the year end */
    unsigned _YearBoundary;
public:
//------------------------------------------------------------------------------------------------------------
Layer2DWithTime(float* data,std::vector<float>lon,std::vector<float>lat,std::vector<float>time):Layer2D(data,lon,lat),_time(time){
       _YearBoundary=365-(unsigned(_time[time.size()-1])%365) + unsigned(_time[time.size()-1]);//units of time assumed days
}
//------------------------------------------------------------------------------------------------------------
float& operator()(unsigned lonix,unsigned latix,unsigned tix=0){
        assert (lonix<_lon.size());
        assert (latix<_lat.size());
        assert (tix  <_time.size());
        return _data[lonix+latix*_lon.size()+ tix  * _lon.size() * _lat.size()];
}
//------------------------------------------------------------------------------------------------------------
double bilin(const double lon, const double lat,const unsigned lonix,const unsigned latix,const unsigned tix=0){
    unsigned lal=latix,lau=latix+1;
    if (lau>=_lat.size())lau=0;
    unsigned lol=lonix,lou=lonix+1;
    if (lou>=_lon.size())lou=_lon.size()-1;
    
    double data_ll=(*this)(lol,lal,tix);
    double data_lu=(*this)(lol,lau,tix);
    double data_ul=(*this)(lou,lal,tix);
    double data_uu=(*this)(lou,lau,tix);
    
    double data_av=0,dn=0;
    if (data_ll!=Constants::cMissingValue){data_av+=data_ll;dn+=1;}  
    if (data_lu!=Constants::cMissingValue){data_av+=data_lu;dn+=1;}  
    if (data_ul!=Constants::cMissingValue){data_av+=data_ul;dn+=1;}  
    if (data_uu!=Constants::cMissingValue){data_av+=data_uu;dn+=1;}  
    if (dn!=0)data_av=data_av/dn;else data_av=Constants::cMissingValue;
    if (data_ll==Constants::cMissingValue){data_ll=data_av;}  
    if (data_lu==Constants::cMissingValue){data_lu=data_av;}  
    if (data_ul==Constants::cMissingValue){data_ul=data_av;}  
    if (data_uu==Constants::cMissingValue){data_uu=data_av;}
    assert(data_ll != Constants::cMissingValue && data_lu != Constants::cMissingValue && data_ul != Constants::cMissingValue && data_ll != Constants::cMissingValue );

    double result=
     data_ll*(_lon[lou]-lon)*(_lat[lau]-lat)
    -data_ul*(_lon[lol]-lon)*(_lat[lau]-lat)
    -data_lu*(_lon[lou]-lon)*(_lat[lal]-lat)
    +data_uu*(_lon[lol]-lon)*(_lat[lal]-lat);
    double dlon=_lon[lou]-_lon[lol],dlat=_lat[lau]-_lat[lal];
    if (dlon==0)dlon=1;
    if (dlat==0)dlat=1;
    result=result/dlon/dlat;
    return result;
}
//------------------------------------------------------------------------------------------------------------
float GetDataAtLonLat(const double Longitude,
                          const double Latitude){
        
        unsigned lonix=binaryChop(_lon,Longitude);
        unsigned latix=binaryChop(_lat,Latitude);
        
        double day=TimeStep::instance()->CurrentDataAccessDay();
        
        unsigned timeIndexu,timeIndexl;
        double tu,tl;
        if (day >= _time[_time.size()-1])day=unsigned(day) % _YearBoundary + (day-unsigned(day));
        if (day<_time[0]) day=day+_YearBoundary;
        if (day >= _time[_time.size()-1]){
            tu=_YearBoundary+_time[0];
            tl=_time[_time.size()-1];
            timeIndexl=_time.size()-1;
            timeIndexu=0;
        }else{
            timeIndexl=binaryChop(_time,day);
            timeIndexu=(timeIndexl + 1);
            tu=_time[timeIndexu];
            tl=_time[timeIndexl];
        }
        
        return ( (*this)(lonix,latix,timeIndexl) * (tu-day) + 
                 (*this)(lonix,latix,timeIndexu) * (day-tl) )/ ( tu - tl );
        //return ( bilin(Longitude,Latitude,lonix,latix,timeIndexl) * (tu-day) + 
        //         bilin(Longitude,Latitude,lonix,latix,timeIndexu) * (day-tl) )/ ( tu - tl );

        
}
//------------------------------------------------------------------------------------------------------------
//get a year of data starting with 0 as the year origin: y is in years. Return empty if beyond last time
std::vector<double> GetYearAtLonLat(const unsigned y,
                                    const double Longitude,
                                    const double Latitude){
            
            unsigned lonix = binaryChop(_lon,Longitude);
            unsigned latix = binaryChop(_lat,Latitude);  
            unsigned tix   = binaryChop(_time,y * 365);
            std::vector<double> result;

            if (tix>=_time.size()-1)return result;

            while (tix<_time.size() && _time[tix]<(y+1)*365){
                result.push_back((*this)(lonix,latix,tix));
                tix++;
            }
            return result;
}
};

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

#endif
