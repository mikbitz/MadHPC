#ifndef LAYER_H
#define	LAYER_H

#include "TimeStep.h"
#include "Maths.h"
#include "Parameters.h"

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
    bool _spatialInterpolation;
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
    Layer(float* data,bool spatialInterpolation=true):_data(data),_spatialInterpolation(spatialInterpolation){}
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
                                  const double,
                                        double day=-1) = 0;      
/**
* \brief Find the location of a point in an (ordered) vector using binary chopping
*
* This is used to return data values at a given point (Longitude,Latitude)
* by indexing each co-ordinate separately. NB co-ordinates must be monotonic!
https://en.wikipedia.org/wiki/Bert_Bithell
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
        //return max val for lo if data is beyond coord[coord.size()-1]
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
        std::vector<float>lat,
        bool spatialInterpolation=true):Layer(data,spatialInterpolation),_lon(lon),_lat(lat){}
//------------------------------------------------------------------------------------------------------------
float& operator()(unsigned lonix,unsigned latix,unsigned tix=0){
         assert (lonix<_lon.size());
         assert (latix<_lat.size());
         return _data[lonix+latix*_lon.size()];
}
//------------------------------------------------------------------------------------------------------------
double bilin(const double lon, const double lat,const unsigned lonix,const unsigned latix,const unsigned tix=0){
    unsigned lal=latix,lau=latix+1;
    if (lau>=_lat.size())lau=latix;
    unsigned lol=lonix,lou=lonix+1;
    unsigned ladj=0;
    //make sure the right data is used if wrapped - assumed that when beyond longitude bounds, only beyond upper bound 
    if (lou>=_lon.size()){
     if (Parameters::instance()->GetNoLongitudeWrap( )){ladj=-1;}else ladj=-lonix-1;
    }
    
    double data_ll=(*this)(lol,lal);
    double data_lu=(*this)(lol,lau);
    double data_ul=(*this)(lou+ladj,lal);
    double data_uu=(*this)(lou+ladj,lau);
    return bilnearWithMissing(data_ll,data_lu,data_ul,data_uu,lol,lal,lou,lau,lon,lat);
}
//------------------------------------------------------------------------------------------------------------
double bilnearWithMissing(double data_ll,double data_lu,double data_ul,double data_uu,unsigned lol,unsigned lal,unsigned lou,unsigned lau,double lon,double lat){
    double data_av=0,dn=0;

    double lplus=0;
    //wrapping may fall in the gap between max. lat and min lat...
    if (lou>=_lon.size()){
     if (Parameters::instance()->GetNoLongitudeWrap( )){lou=lou-1;}else {lou=0;lplus=360;}
    }
    if (lat<_lat[lal])lat=_lat[lal];
    if (lat>_lat[lau])lat=_lat[lau];
    //simplistic way to deal with miising data - but using cut quarter cells might be better...
    if (data_ll>Constants::cMissingValue){data_av+=data_ll;dn+=1;}else {if ((_lon[lol]-lon)==0 && (_lat[lal]-lat)==0)return Constants::cMissingValue;}  
    if (data_lu>Constants::cMissingValue){data_av+=data_lu;dn+=1;}else {if ((_lon[lol]-lon)==0 && (_lat[lau]-lat)==0)return Constants::cMissingValue;}  
    if (data_ul>Constants::cMissingValue){data_av+=data_ul;dn+=1;}else {if ((_lon[lou]+lplus-lon)==0 && (_lat[lal]-lat)==0)return Constants::cMissingValue;}  
    if (data_uu>Constants::cMissingValue){data_av+=data_uu;dn+=1;}else {if ((_lon[lou]+lplus-lon)==0 && (_lat[lau]-lat)==0)return Constants::cMissingValue;} 
    if (dn!=0)data_av=data_av/dn;else data_av=Constants::cMissingValue;
    if (data_ll<=Constants::cMissingValue){data_ll=data_av;}  
    if (data_lu<=Constants::cMissingValue){data_lu=data_av;}  
    if (data_ul<=Constants::cMissingValue){data_ul=data_av;}  
    if (data_uu<=Constants::cMissingValue){data_uu=data_av;}

     
   
    double dlon=abs(_lon[lou]+lplus-_lon[lol]),dlat=abs(_lat[lau]-_lat[lal]);
    //ensure nor division by zero!
    if (dlon==0 && dlat==0) return  data_ll;
    if (dlat==0)            return (data_ll*abs(_lon[lou]+lplus-lon)+data_ul*abs(_lon[lol]-lon))/dlon;
    if (dlon==0)            return (data_ll*abs(_lat[lau]-lat)      +data_lu*abs(_lat[lal]-lat))/dlat;

    double result=0;    
    result= data_ll*abs((_lon[lou]+lplus-lon)*(_lat[lau]-lat))
           +data_ul*abs((_lon[lol]-lon)      *(_lat[lau]-lat))
           +data_lu*abs((_lon[lou]+lplus-lon)*(_lat[lal]-lat))
           +data_uu*abs((_lon[lol]-lon)      *(_lat[lal]-lat));

    result=result/dlon/dlat;
    return result;
}
//------------------------------------------------------------------------------------------------------------
float GetDataAtLonLat(const double Longitude,
                      const double Latitude,
                            double day=-1){
    //binary chop will return 0 if beyond lower bound: no approaptiate for longitude wrapping so add 360 
    double lplus=0;
    if (Longitude< _lon[0] && !Parameters::instance()->GetNoLongitudeWrap( ))lplus=360;
    unsigned lonix=binaryChop(_lon,Longitude+lplus);
    unsigned latix=binaryChop(_lat,Latitude);
    if (_spatialInterpolation){
        return bilin(Longitude+lplus,Latitude,lonix,latix);
    }else{
        return (*this)(lonix,latix);
    }
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
Layer2DWithTime(float* data,
                std::vector<float>lon,
                std::vector<float>lat,
                std::vector<float>time,
                bool spatialInterpolation=true):Layer2D(data,lon,lat,spatialInterpolation),
                _time(time)
{
       //how many days are there between the end of the dataset and the end of the year following the end of the dataset?
       //add this to the final time in order to find the year-end just beyond the end of the dataset
       //needed in case we have to interpolate including a wrapping in time, e.g. for climatology
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
    if (lau>=_lat.size())lau=latix;
    unsigned lol=lonix,lou=lonix+1;
    unsigned ladj=0;
    if (lou>=_lon.size()){
     if (Parameters::instance()->GetNoLongitudeWrap( )){ladj=-1;}else ladj=-lonix-1;
    }
    double data_ll=(*this)(lol,lal,tix);
    double data_lu=(*this)(lol,lau,tix);
    double data_ul=(*this)(lou+ladj,lal,tix);
    double data_uu=(*this)(lou+ladj,lau,tix);
    return bilnearWithMissing(data_ll,data_lu,data_ul,data_uu,lol,lal,lou,lau,lon,lat);
}
//------------------------------------------------------------------------------------------------------------
float GetDataAtLonLat(const double Longitude,
                      const double Latitude,
                            double day=-1  ){
        double lplus=0;
        if (Longitude< _lon[0] && !Parameters::instance()->GetNoLongitudeWrap( ))lplus=360;
        unsigned lonix=binaryChop(_lon,Longitude+lplus);
        unsigned latix=binaryChop(_lat,Latitude);
        
        if (day<0) day=TimeStep::instance()->CurrentDataAccessDay();
        
        unsigned timeIndexu,timeIndexl;
        double tu,tl;
        //set time indices allowing for day to be any distance beyond the end of time
        //(e.g. day could be 100*365 in a climatology run of 100 years with only 1 year of (repeated) data)
        //also allow for day to possibly be fractional, and for t[0] to be at any arbitrary day
        if (day<_time[0])day=day+_YearBoundary;
        if (day >= _time[_time.size()-1])day=unsigned(day) % _YearBoundary + (day-unsigned(day));
        if (day<unsigned(_time[0])%365+(_time[0]-unsigned(_time[0]))) day=day+_YearBoundary;
        if (day >= _time[_time.size()-1]){
            //wrap across the end of the dataset in time
            tu=_YearBoundary+unsigned(_time[0])%365+(_time[0]-unsigned(_time[0]));
            tl=_time[_time.size()-1];
            timeIndexl=_time.size()-1;
            timeIndexu=0;
        }else{
            timeIndexl=binaryChop(_time,day);
            timeIndexu=(timeIndexl + 1);
            tu=_time[timeIndexu];
            tl=_time[timeIndexl];
        }
        if (_spatialInterpolation){
            return ( bilin(Longitude+lplus,Latitude,lonix,latix,timeIndexl) * (tu-day) + 
                     bilin(Longitude+lplus,Latitude,lonix,latix,timeIndexu) * (day-tl) )/ ( tu - tl );
        }else{
            return ( (*this)(lonix,latix,timeIndexl) * (tu-day) + 
                     (*this)(lonix,latix,timeIndexu) * (day-tl) )/ ( tu - tl );
        }

        
}
//------------------------------------------------------------------------------------------------------------
//get a year of data starting with 0 as the year origin: y is in years. Return empty if beyond last time
std::vector<double> GetYearAtLonLat(const unsigned y,
                                    const double Longitude,
                                    const double Latitude){
            double lplus=0;
            if (Longitude< _lon[0] && !Parameters::instance()->GetNoLongitudeWrap( ))lplus=360;
            unsigned lonix=binaryChop(_lon,Longitude+lplus);            
            unsigned latix = binaryChop(_lat,Latitude);  
            unsigned tix   = binaryChop(_time,y * 365);
            std::vector<double> result;

            if (tix>=_time.size()-1)return result;

            while (tix<_time.size() && _time[tix]<(y+1)*365){
                if (_spatialInterpolation){
                  result.push_back(bilin(Longitude+lplus,Latitude,lonix,latix,tix));
                }else{
                  result.push_back((*this)(lonix,latix,tix));
                }
                tix++;
            }
            return result;
}
};

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
class layerTester{

public:
static void tests(){
    float* data;
    std::vector<float>lon {-120,-60,0,60,120,180};
    std::vector<float>lat {-5,-4,-3,-2,-1,0,1,2,3,4,5};
    data=new float(lon.size()*lat.size());
    Layer* L=new Layer2D(data,lon,lat);
    std::cout<<"Layers test 1:binary chop ";
    std::vector<float>Ar1 {1,2,3,4,5};
    std::vector<float>Ar2 {-5,-4,-3,-2,-1,0,1,2,3,4,5};
    std::vector<float>Ar3 {6,3.7,1.5,0.8,-4,-100};
    assert (L->binaryChop(Ar1,-11)==0);
    assert (L->binaryChop(Ar1,0)==0);
    assert (L->binaryChop(Ar1,1)==0);
    assert (L->binaryChop(Ar1,1.1)==0);
    assert (L->binaryChop(Ar1,2)==1);
    assert (L->binaryChop(Ar1,2.7)==1);
    assert (L->binaryChop(Ar1,3)==2);
    assert (L->binaryChop(Ar1,3.01)==2);
    assert (L->binaryChop(Ar1,3.9999)==2);
    assert (L->binaryChop(Ar1,4)==3);
    assert (L->binaryChop(Ar1,4.8312)==3);
    assert (L->binaryChop(Ar1,5)==4);
    assert (L->binaryChop(Ar1,9)==4);
    assert (L->binaryChop(Ar2,-6)==0);
    assert (L->binaryChop(Ar2,-4.9)==0);
    assert (L->binaryChop(Ar2,0)==5);
    assert (L->binaryChop(Ar2,-0.001)==4);
    assert (L->binaryChop(Ar2,0.001)==5);
    assert (L->binaryChop(Ar2,5)==10);
    assert (L->binaryChop(Ar3,6)==0);
    assert (L->binaryChop(Ar3,-100)==5);
    assert (L->binaryChop(Ar3,-1000)==5);
    assert (L->binaryChop(Ar3,3.8)==0);
    assert (L->binaryChop(Ar3,3.65)==1);
    assert (L->binaryChop(Ar3,1.4)==2);
    assert (L->binaryChop(Ar3,0.99)==2);
    assert (L->binaryChop(Ar3,0.75)==3);
    assert (L->binaryChop(Ar3,-4)==4);
    assert (L->binaryChop(Ar3,-5)==4);
    assert (L->binaryChop(Ar3,-99)==4);
    std::cout<<"success"<<std::endl;
    std::cout<<"Layers test 2:bilinear interpolation ";
    unsigned lonix=L->binaryChop(lon,0);
    unsigned latix=L->binaryChop(lat,0);

    for (unsigned i=0;i<lon.size();i++)
        for (unsigned j=0;j<lat.size();j++ )(*L)(i,j)=1;
    assert(L->bilin(0,0,lonix,latix)==1);

    for (unsigned i=0;i<lon.size();i++)
        for (unsigned j=0;j<lat.size();j++ )(*L)(i,j)=-0.5;
    assert(L->bilin(0,0,lonix,latix)==-0.5);

    for (unsigned i=0;i<lon.size();i++)
        for (unsigned j=0;j<lat.size();j++ ){(*L)(i,j)=i;}

    lonix=L->binaryChop(lon,180);
    assert(L->bilin(+180,0,lonix,latix)==5);
    
    lonix=L->binaryChop(lon,0);
    assert(L->bilin(0,0,lonix,latix)==2);
    
    lonix=L->binaryChop(lon,150);
    assert(L->bilin(150,0,lonix,latix)==4.5);
        
    lonix=L->binaryChop(lon,-120);
    assert(L->bilin(-120,0,lonix,latix)==0);
    
    lonix=L->binaryChop(lon,210);
    assert(L->bilin(210,0,lonix,latix)==2.5);
    
    lonix=L->binaryChop(lon,-90);
    assert(L->bilin(-90,0,lonix,latix)==0.5);
    
    lonix=L->binaryChop(lon,90);
    latix=L->binaryChop(lat,4.15);
    assert(L->bilin(90,4.15,lonix,latix)==3.5);
    
    lonix=L->binaryChop(lon,135);
    latix=L->binaryChop(lat,-9);

    assert(L->bilin(135,-9,lonix,latix)==4.25);

    for (unsigned i=0;i<lon.size();i++)
        for (unsigned j=0;j<lat.size();j++ )(*L)(i,j)=2*(j-5)+10;

    lonix=L->binaryChop(lon,135);
    latix=L->binaryChop(lat,0);
    assert(L->bilin(135,0,lonix,latix)==10);

    lonix=L->binaryChop(lon,175);
    latix=L->binaryChop(lat,-5);
    assert(L->bilin(175,-5,lonix,latix)==0);


    lonix=L->binaryChop(lon,175);
    latix=L->binaryChop(lat,-6);
    assert(L->bilin(175,-6,lonix,latix)==0);
    
    lonix=L->binaryChop(lon,-113);
    latix=L->binaryChop(lat,4);
    assert(L->bilin(-113,4,lonix,latix)==18);
   
    lonix=L->binaryChop(lon,77);
    latix=L->binaryChop(lat,5);
    assert(L->bilin(77,5,lonix,latix)==20);
    
    assert(L->GetDataAtLonLat(77,5)==20);
    assert(L->GetDataAtLonLat(-120,3)==16);
    
    std::cout<<"success"<<std::endl;
    std::cout<<"Layers test 3:time varying fields ";

    std::vector<float>time = {130,230,330,430,530,630};
    delete L;
    data=new float(lon.size()*lat.size()*time.size());
    L=new Layer2DWithTime(data,lon,lat,time);
    for (unsigned i=0;i<lon.size();i++)
        for (unsigned j=0;j<lat.size();j++ )
            for (unsigned t= 0;t<time.size();t++)(*L)(i,j,t)=t;
    assert(L->GetDataAtLonLat(-120,3,130)==0);
    assert(L->GetDataAtLonLat(-110,0,230)==1);
    assert(L->GetDataAtLonLat( 120,5,530)==4);
    assert(L->GetDataAtLonLat(   0,0,180)==0.5);
    assert(L->GetDataAtLonLat(   0,0,305)==1.75);
    assert(L->GetDataAtLonLat(   0,0,630)==5);
    //630 is 100 days from year boundary - 115 days after 630 should be half way to wrapped year where the first data value is at 130
    assert(L->GetDataAtLonLat(   0,0,630+115)==2.5);
    assert(L->GetDataAtLonLat(   0,0,100*365+130)==0);    
    std::cout<<"success"<<std::endl;
    delete L;
}
    
};

#endif
