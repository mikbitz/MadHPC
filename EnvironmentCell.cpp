/*
 *
 * EnvironmentCell.cpp
 *
 *  Created on: May 7, 2018
 *      Author: Mike Bithell
 * 
 *  Derived from Original C# code by
 * Emergent Global Patterns of Ecosystem Structure and Function from a Mechanistic General Ecosystem Model , Michael B. J. Harfoot , Tim Newbold, Derek P. Tittensor,  Stephen Emmott, Jon Hutton, Vassily Lyutsarev, Matthew J. Smith, Jörn P. W. Scharlemann, Drew W. Purves PLOS Biology 12(4) 2014 e1001841 https://doi.org/10.1371/journal.pbio.1001841
 */
#include "EnvironmentCell.h"
#include "CohortMerger.h"

#include "ClimateVariablesCalculator.h"
#include "DataLayerSet.h"
#include "TimeStep.h"
#include "Parameters.h"
#include "UtilityFunctions.h"
#include "Stock.h"


EnvironmentCell::EnvironmentCell() {}

//------------------------------------------------------------------------------
EnvironmentCell::EnvironmentCell(int x,int y):_x(x),_y(y){
    UtilityFunctions Utility;
    _Area=   Utility.CalculateGridCellArea(Latitude(),cellSize());//in sqkm
    _Width=  Utility.CalculateLengthOfDegreeLatitude( Latitude()) *cellSize();
    _Height= Utility.CalculateLengthOfDegreeLongitude( Latitude())*cellSize();

    _OrganicPool=0;
    _RespiratoryCO2Pool=0;

    SetRealm( );
    SetTotalPrecip();
    SetAVGSDTemp();
    SetNPPSeasonality();
    SetBreeding();
    SetFrostandFire();

}
//------------------------------------------------------------------------------

void EnvironmentCell::SetRealm( ) {
    _Realm="none";
    if( DataLayerSet::Data( )->GetDataAtLonLatFor( "Realm", Longitude(),  Latitude() ) <= 1.5 ) {
          _Realm="marine";
    } else if( DataLayerSet::Data( )->GetDataAtLonLatFor( "Realm", Longitude(),  Latitude() ) > 1.5) {
          _Realm="terrestrial";
    }

}
//------------------------------------------------------------------------------
void EnvironmentCell::SetTotalPrecip(){
 double d = Constants::cMissingValue;
 if( _Realm=="terrestrial" ) {

 auto precip=DataLayerSet::Data( )->GetLayer("TerrestrialPre");
 // allow for multi-year data
 unsigned year=0;
 auto p=((Layer2DWithTime*)(precip))->GetYearAtLonLat( year, Longitude(),  Latitude() );
 while (p.size()>0){
    d=0;
    for (auto& pre:p)if(pre > Constants::cMissingValue)d+=pre;
     year++;
     p=((Layer2DWithTime*)(precip))->GetYearAtLonLat( year, Longitude(),  Latitude() );
 }
 _TotalPrecip.push_back(d);
 }

}
//------------------------------------------------------------------------------
void EnvironmentCell::SetAVGSDTemp(){
    Layer* Temp=nullptr;
    if( _Realm=="marine" ) {
        Temp=DataLayerSet::Data( )->GetLayer("MarineTemp");
    } else if( _Realm=="terrestrial" ) {
        Temp=DataLayerSet::Data( )->GetLayer("TerrestrialTemp");
    }
    assert (Temp!=nullptr);
    unsigned year=0;
    auto T=((Layer2DWithTime*)(Temp))->GetYearAtLonLat( year, Longitude(),  Latitude() );
    while (T.size()>0){
        double avg=0;
        for (auto& t:T){
            if( t != Constants::cMissingValue ) avg+=t;
        }
        avg=avg/T.size();
        double sota = 0, sumExp = 0;
        _exptdev.push_back(vector<double>(T.size(),0));
        for (unsigned i=0;i<T.size();i++){
            if( T[i] != Constants::cMissingValue ) {
                sota += ( T[i] - avg )*( T[i] - avg );
                _exptdev[year][i] = exp( -( T[i] - avg ) / 3 );
                sumExp += _exptdev[year][i];
            }
        }
        for( auto& e:_exptdev[year]) {
            e = e / sumExp;
        }
        _AnnualTemperature.push_back(avg);
        _SDTemperature.push_back(sqrt( sota / 12 ));
        year++;
        T=((Layer2DWithTime*)(Temp))->GetYearAtLonLat( year, Longitude(),  Latitude() );
    }

}
//----------------------------------------------------------------------------------------------

/** \brief Calculate monthly seasonality values of Net Primary Production - ignores missing values. If there is no NPP data (ie all zero or missing values)
then assign 1/12 for each month.
 */
void EnvironmentCell::SetNPPSeasonality( ) {
    Layer* N=nullptr;
    if( _Realm=="marine" ) {
        N= DataLayerSet::Data( )->GetLayer( "MarineNPP");
    } else if( _Realm=="terrestrial" ) {
        N= DataLayerSet::Data( )->GetLayer("TerrestrialNPP");
    }
    assert(N!=nullptr);
    unsigned year=0;
    auto NPP=((Layer2DWithTime*)(N))->GetYearAtLonLat( year, Longitude(),  Latitude() );
    while (NPP.size()!=0){
        _Seasonality.push_back(vector<double>(NPP.size(),0));
        // Loop over months and calculate total annual NPP
        
        double totalNPP = 0.0;
        for( auto& n:NPP ) {
            
            if( n != Constants::cMissingValue && n > 0 ) totalNPP += n;
        }
        if( totalNPP == 0 ) {
            // Loop over months and calculate seasonality
            // If there is no NPP value then assign a uniform flat seasonality
            for( auto& s:_Seasonality[year]) {
                s = 1.0 / 12.0;
            }
        } else {
            // Some NPP data exists for this grid cell so use that to infer the NPP seasonality
            // Loop over months and calculate seasonality
            for( int i = 0; i < NPP.size(); i++ ) {
                double n = NPP[i];
                if( n != Constants::cMissingValue && n > 0 ) {
                    _Seasonality[year][ i ] = n / totalNPP;
                } else {
                    _Seasonality[year][ i ] = 0.0;
                }
            }
        }
        year++;
        NPP=((Layer2DWithTime*)(N))->GetYearAtLonLat( year, Longitude(),  Latitude() );
    }
}
//----------------------------------------------------------------------------------------------
void EnvironmentCell::SetBreeding( ) {
    // Designate a breeding season for this grid cell, where a month is considered to be part of the breeding season if its NPP is at
    // least 80% of the maximum NPP throughout the whole year
    double maxSeason = -1;
    for (unsigned year=0;year<_Seasonality.size();year++){
        _Breeding_Season.push_back(vector<double>(_Seasonality[year].size(),0));
        for( int timeIndex = 0; timeIndex < _Seasonality[year].size(); timeIndex++ ) {
            maxSeason = std::max( maxSeason, _Seasonality[year][timeIndex]);
        }
        for( int timeIndex = 0; timeIndex < _Seasonality[year].size(); timeIndex++ ) {
            if(  _Seasonality[year][timeIndex] / maxSeason > 0.5 ) {
                _Breeding_Season[year][timeIndex] = 1.0;
            } else {
                _Breeding_Season[year][timeIndex] = 0.0;
            }
        }
    }
}
//----------------------------------------------------------------------------------------------
void EnvironmentCell::SetFrostandFire( ) {
    // Calculate other climate variables from temperature and precipitation
    // Declare an instance of the climate variables calculator
    if( _Realm=="terrestrial" ){
        ClimateVariablesCalculator CVC;
        auto Frost = DataLayerSet::Data( )->GetLayer("TerrestrialFrost");
        auto pre   = DataLayerSet::Data( )->GetLayer("TerrestrialPre");
        auto T     = DataLayerSet::Data( )->GetLayer("TerrestrialTemp");
        
        double AWC = DataLayerSet::Data()->GetDataAtLonLatFor( "TerrestrialAWC", Longitude(),  Latitude() );
        unsigned year=0;
        auto FrostDays=((Layer2DWithTime*)(Frost))->GetYearAtLonLat( year, Longitude(),  Latitude() );
        auto Temp     =((Layer2DWithTime*)(T    ))->GetYearAtLonLat( year, Longitude(),  Latitude() );
        auto Precip   =((Layer2DWithTime*)(pre  ))->GetYearAtLonLat( year, Longitude(),  Latitude() );

        while(Precip.size()!=0){
            _FractionYearFrost.push_back(CVC.GetNDF( FrostDays, Temp, Constants::cMissingValue ));
            
            std::tuple< std::vector< double >, double, double > TempTuple = CVC.MonthlyActualEvapotranspirationSoilMoisture( AWC, Precip, Temp );
            _TotalAET.push_back(0);
            _AET.push_back(vector<double>(Precip.size(),0));
            for( int timeIndex = 0; timeIndex < Precip.size(); timeIndex++ ) {
                _AET[year][timeIndex ] = std::get< 0 >( TempTuple )[ timeIndex ];
                _TotalAET[year] += std::get< 0 >( TempTuple )[ timeIndex ];
            }
            _FractionYearFire.push_back( ( std::get< 2 > ( TempTuple ) / 360 ));
            year++;
            FrostDays=((Layer2DWithTime*)(Frost))->GetYearAtLonLat( year, Longitude(),  Latitude() );
            Temp     =((Layer2DWithTime*)(T    ))->GetYearAtLonLat( year, Longitude(),  Latitude() );
            Precip   =((Layer2DWithTime*)(pre  ))->GetYearAtLonLat( year, Longitude(),  Latitude() );
        }
    }
}
//----------------------------------------------------------------------------------------------

double EnvironmentCell::Width(){return _Width;}
double EnvironmentCell::Height(){return _Height;}
double EnvironmentCell::Area(){return  _Area;}
//----------------------------------------------------------------------------------------------
double EnvironmentCell::TotalAET(){return _TotalAET[TimeStep::instance()->CurrentYear()];}
double EnvironmentCell::AET(){return _AET[TimeStep::instance()->CurrentYear()][TimeStep::instance( )->CurrentMonth()];}
//----------------------------------------------------------------------------------------------
double EnvironmentCell::FractionYearFrost(){return _FractionYearFrost[TimeStep::instance()->CurrentYear()];}
double EnvironmentCell::FractionYearFire(){return _FractionYearFire[TimeStep::instance()->CurrentYear()];}
double EnvironmentCell::Breeding_Season(){assert(TimeStep::instance( )->CurrentYear()<1);return _Breeding_Season[TimeStep::instance()->CurrentYear()][TimeStep::instance( )->CurrentMonth()];}
//------------------------------------------------------------------------------
double EnvironmentCell::Seasonality(){ return _Seasonality[TimeStep::instance()->CurrentYear()][TimeStep::instance( )->CurrentMonth()];}
//------------------------------------------------------------------------------
double EnvironmentCell::ExpTDevWeight(){ return _exptdev[TimeStep::instance()->CurrentYear()][TimeStep::instance( )->CurrentMonth()];}
//------------------------------------------------------------------------------
double EnvironmentCell::TerrestrialHANPP(){
 double d=0;
 if( _Realm=="terrestrial" ) {
   d = GetVariableFromDatasetNamed( "TerrestrialHANPP");
 }
 return d;
}
//------------------------------------------------------------------------------
double EnvironmentCell::Temperature(){

 double d = 0;
 
 if( _Realm=="marine" ) {
    d = GetVariableFromDatasetNamed("MarineTemp");
 } else if( _Realm=="terrestrial" ) {
    d = GetVariableFromDatasetNamed( "TerrestrialTemp");
 }

 return d;
}
//------------------------------------------------------------------------------
double EnvironmentCell::uVel(){

  double d = 0;

  if( _Realm=="marine"  ) {
    d = GetVariableFromDatasetNamed("MarineEastVel");
  }

  return d;
}
//------------------------------------------------------------------------------
double EnvironmentCell::vVel(){
    
  double d = 0;

  if( _Realm=="marine"  ) {
    d = GetVariableFromDatasetNamed("MarineNorthVel");
  }

  return d;
}
//------------------------------------------------------------------------------
double EnvironmentCell::DiurnalTemperatureRange(){

 double d = 0;
 if( _Realm=="terrestrial" ) {
   d = GetVariableFromDatasetNamed("TerrestrialDTR");
  }
  return d;
    
}
//------------------------------------------------------------------------------
double EnvironmentCell::NPP( ) {
 double d = 0;
 if( _Realm=="marine" ) {
   d= GetVariableFromDatasetNamed( "MarineNPP");
 } else if( _Realm=="terrestrial" ) {
   d= GetVariableFromDatasetNamed("TerrestrialNPP");
 }
 return d;
}
//------------------------------------------------------------------------------
double EnvironmentCell::Precipitation(){

 double d = Constants::cMissingValue;//currently no marine precip
 if( _Realm=="terrestrial" ) {
   d = GetVariableFromDatasetNamed("TerrestrialPre");
  }
  return d;
    
}
//-----------------------------------------------------------------------------
double EnvironmentCell::TotalPrecip(){return _TotalPrecip[TimeStep::instance()->CurrentYear()];}
//------------------------------------------------------------------------------
double EnvironmentCell::AnnualTemperature(){return _AnnualTemperature[TimeStep::instance()->CurrentYear()];}
//------------------------------------------------------------------------------
double EnvironmentCell::SDTemperature(){return _SDTemperature[TimeStep::instance()->CurrentYear()];}
//------------------------------------------------------------------------------
double EnvironmentCell::GetVariableFromDatasetNamed(std:: string s){

    double d = Constants::cMissingValue;
    d = DataLayerSet::Data( )->GetDataAtLonLatFor( s, Longitude(),  Latitude() );
    if( d == Constants::cMissingValue ) {
      std::cout << "Warning EnvironmentCell::GetVariableFromDatasetNamed- missing values in "<<s<<" field!! "<<Longitude()<<" E "<<Latitude()<<" N "<< std::endl;
    }
    return d;
}
//------------------------------------------------------------------------------
double EnvironmentCell::Latitude(){
    return Parameters::instance()->GetLatitudeAtIndex(_y);
}
//------------------------------------------------------------------------------
double EnvironmentCell::Longitude(){
    return Parameters::instance()->GetLongitudeAtIndex(_x);
}
//------------------------------------------------------------------------------
double EnvironmentCell::cellSize(){
    return Parameters::instance()->GetGridCellSize();
}



