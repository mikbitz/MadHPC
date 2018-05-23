#include "Environment.h"
#include "CohortMerger.h"

#include "ClimateVariablesCalculator.h"

#include "DataLayerSet.h"
#include "DataIndices.h"
#include "TimeStep.h"
#include "Parameters.h"
#include "UtilityFunctions.h"
#include "Stock.h"

bool fn(Cohort* a, Cohort* b){return a->_FunctionalGroupIndex * 10000000 + a->getId().id() <b->_FunctionalGroupIndex * 10000000 + b->getId().id();}

void Environment::merge(int& bap){
  bap+=CohortMerger::MergeToReachThresholdFast( this );
}

void Environment::step(){
    repast::relogo::AgentSet<Cohort> cohorts= turtlesHere<Cohort>();

    repast::relogo::AgentSet<Stock> stocks= turtlesHere<Stock>();
    //DEBUG
   // if (getId().id()==0 && getId().startingRank()==0){
         // cout<<"Eeee "<<xCor()<< " " <<yCor()<<endl;   

       // cout<<"Eeee "<<getId().id()<<" "<<getId().agentType()<<" ook ";
       // for (auto c:cohorts){cout<<c->getId().id()<<" "<<c->getId().agentType()<<" "<<c->xCor()<<" "<<c->yCor()<<endl;cout.flush();}
   // }
    //random_shuffle(cohorts.begin(),cohorts.end());
    //DEBUG
    //sort(cohorts.begin(),cohorts.end(),fn);

    for (auto c:cohorts){if(c->getId().currentRank()==repast::RepastProcess::instance()->rank())c->step(this,cohorts,stocks);}
}
void Environment::flap(){
    repast::relogo::AgentSet<Cohort> cohorts= turtlesHere<Cohort>();

    for (auto c:cohorts){if(c->getId().currentRank()==repast::RepastProcess::instance()->rank())c->moveIt();}
}
Environment::Environment(repast::AgentId id, repast::relogo::Observer* obs) : repast::relogo::Patch(id, obs) {}
//------------------------------------------------------------------------------
void Environment::setup(){
    _cellIndex=pxCor()+pyCor()*Parameters::Get()->GetLengthUserLongitudeArray( );//based on 1-D representation of 2D arrays
    name="slartifartbarst";
    UtilityFunctions Utility;
    _Area=   Utility.CalculateGridCellArea(Latitude(),cellSize());
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
    //make sure all time dependent fields set to the start

    TimeStep::Get( )->SetMonthly( 0);
}
//------------------------------------------------------------------------------

void Environment::SetRealm( ) {
    _Realm="none";
    if( DataLayerSet::Get( )->GetDataAtCellIndexFor( "Realm", _cellIndex ) == 1 ) {
          _Realm="marine";
    } else if( DataLayerSet::Get( )->GetDataAtCellIndexFor( "Realm", _cellIndex ) == 2 ) {
          _Realm="terrestrial";
    }
}
//------------------------------------------------------------------------------
void Environment::SetTotalPrecip(){
 double d = Constants::cMissingValue;
 if( _Realm=="terrestrial" ) {
  d=0;
  for( int timeIndex = 0; timeIndex < 12; timeIndex++ ) {
    TimeStep::Get( )->SetMonthly( timeIndex );
    d+=Precipitation();
  }
 }
 _TotalPrecip=d;
}
//------------------------------------------------------------------------------
void Environment::SetAVGSDTemp(){

 double avg = 0;
 for( int timeIndex = 0; timeIndex < 12; timeIndex++ ) {
     double d = Constants::cMissingValue;
    TimeStep::Get( )->SetMonthly( timeIndex );
    d=Temperature();
    if( d == Constants::cMissingValue ) d = 0;
    avg += d;
 }
 
 avg = avg / 12;
 double sota = 0, sumExp = 0;
 _exptdev.resize(12);
 for( int timeIndex = 0; timeIndex < 12; timeIndex++ ) {
    double d = Constants::cMissingValue;
    TimeStep::Get( )->SetMonthly( timeIndex );
    d=Temperature();

    if( d == Constants::cMissingValue )d = 0;
    sota += ( d - avg )*( d - avg );
    _exptdev[ timeIndex ] = exp( -( d - avg ) / 3 );
    sumExp += _exptdev[ timeIndex ];
 }
 for( int tm = 0; tm < 12; tm++ ) {
    _exptdev[tm] = _exptdev[tm] / sumExp;
 }
 _AnnualTemperature = avg;
 _SDTemperature=sqrt( sota / 12 );

}
//----------------------------------------------------------------------------------------------

/** \brief Calculate monthly seasonality values of Net Primary Production - ignores missing values. If there is no NPP data (ie all zero or missing values)
then assign 1/12 for each month.
 */
void Environment::SetNPPSeasonality( ) {
  _Seasonality.resize(12);
 // Loop over months and calculate total annual NPP
 double totalNPP = 0.0;
 for( int timeIndex = 0; timeIndex < 12; timeIndex++ ) {
     TimeStep::Get( )->SetMonthly( timeIndex );
     double N = NPP();
     if( N != Constants::cMissingValue && N > 0 ) totalNPP += N;
 }
 if( totalNPP == 0 ) {
     // Loop over months and calculate seasonality
     // If there is no NPP value then assign a uniform flat seasonality
     for( int timeIndex = 0; timeIndex < 12; timeIndex++ ) {
        _Seasonality[ timeIndex ] = 1.0 / 12.0;
     }
  } else {
     // Some NPP data exists for this grid cell so use that to infer the NPP seasonality
     // Loop over months and calculate seasonality
     for( int timeIndex = 0; timeIndex < 12; timeIndex++ ) {
         TimeStep::Get( )->SetMonthly( timeIndex );
         double N = NPP();
         if( N != Constants::cMissingValue && N > 0 ) {
             _Seasonality[ timeIndex ] = N / totalNPP;
         } else {
             _Seasonality[ timeIndex  ] = 0.0;
         }
     }
 }

}
//----------------------------------------------------------------------------------------------
void Environment::SetBreeding( ) {
// Designate a breeding season for this grid cell, where a month is considered to be part of the breeding season if its NPP is at
// least 80% of the maximum NPP throughout the whole year
    double maxSeason = -1;
    _Breeding_Season.resize(12);
    for( int timeIndex = 0; timeIndex < 12; timeIndex++ ) {
        maxSeason = std::max( maxSeason, _Seasonality[timeIndex]);
    }
    for( int timeIndex = 0; timeIndex < 12; timeIndex++ ) {
      if(  _Seasonality[timeIndex] / maxSeason > 0.5 ) {
          _Breeding_Season[timeIndex] = 1.0;
      } else {
          _Breeding_Season[timeIndex] = 0.0;
      }
    }

}
//----------------------------------------------------------------------------------------------
void Environment::SetFrostandFire( ) {
    // Calculate other climate variables from temperature and precipitation
    // Declare an instance of the climate variables calculator
    ClimateVariablesCalculator CVC;
    _AET.resize(12);

        // Calculate the fraction of the year that experiences frost
        std::vector< double > FrostDays( 12 ), Temp( 12 ), Precip( 12 );
        for( int timeIndex = 0; timeIndex < 12; timeIndex++ ) {

            TimeStep::Get( )->SetMonthly( timeIndex );
            if(_Realm=="terrestrial" ) {
                FrostDays[timeIndex] = DataLayerSet::Get( )->GetDataAtCellIndexFor( "TerrestrialFrost", _cellIndex );
                Precip[timeIndex] = Precipitation();
                Temp[timeIndex] = Temperature();
            }
        }
        _FractionYearFrost = CVC.GetNDF( FrostDays, Temp, Constants::cMissingValue );

        //std::vector< double > MonthDays = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
        //not used??
        //for( int timeIndex = 0; timeIndex < 12; timeIndex++ ) {
        //    _FractionMonthFrost[timeIndex ] = std::min( FrostDays[timeIndex] / MonthDays[timeIndex], ( double )1.0 );
        //}
        double AWC = DataLayerSet::Get( )->GetDataAtCellIndexFor( "TerrestrialAWC", _cellIndex );
        
        std::tuple< std::vector< double >, double, double > TempTuple = CVC.MonthlyActualEvapotranspirationSoilMoisture( AWC, Precip, Temp );
        _TotalAET = 0;
        for( int timeIndex = 0; timeIndex < 12; timeIndex++ ) {
            _AET[timeIndex ] = std::get< 0 >( TempTuple )[ timeIndex ];
            _TotalAET += std::get< 0 >( TempTuple )[ timeIndex ];
        }
        _FractionYearFire = ( std::get< 2 > ( TempTuple ) / 360 );
    
}
//----------------------------------------------------------------------------------------------

double Environment::Width(){return _Width;}
double Environment::Height(){return _Height;}
double Environment::Area(){return  _Area;}
//----------------------------------------------------------------------------------------------
double Environment::TotalAET(){return _TotalAET;}
double Environment::AET(){return _AET[TimeStep::Get( )->Get(Constants::cMonthlyTimeUnitName)];}
//----------------------------------------------------------------------------------------------
double Environment::FractionYearFrost(){return _FractionYearFrost;}
double Environment::FractionYearFire(){return _FractionYearFire;}
double Environment::Breeding_Season(){ return _Breeding_Season[TimeStep::Get( )->Get(Constants::cMonthlyTimeUnitName)];}
//------------------------------------------------------------------------------
double Environment::Seasonality(){ return _Seasonality[TimeStep::Get( )->Get(Constants::cMonthlyTimeUnitName)];}
//------------------------------------------------------------------------------
double Environment::ExpTDevWeight(){ return _exptdev[TimeStep::Get( )->Get(Constants::cMonthlyTimeUnitName)];}
//------------------------------------------------------------------------------
double Environment::TerrestrialHANPP(){
 double d=0;
 if( _Realm=="terrestrial" ) {
   d = GetVariableFromDatasetNamed( "TerrestrialHANPP");
 }
 return d;
}
//------------------------------------------------------------------------------
double Environment::Temperature(){

 double d = 0;
 
 if( _Realm=="marine" ) {
    d = GetVariableFromDatasetNamed("MarineTemp");
 } else if( _Realm=="terrestrial" ) {
    d = GetVariableFromDatasetNamed( "TerrestrialTemp");
 }

 return d;
}
//------------------------------------------------------------------------------
double Environment::uVel(){

  double d = 0;

  if( _Realm=="marine"  ) {
    d = GetVariableFromDatasetNamed("MarineEastVel");
  }

  return d;
}
//------------------------------------------------------------------------------
double Environment::vVel(){
    
  double d = 0;

  if( _Realm=="marine"  ) {
    d = GetVariableFromDatasetNamed("MarineNorthVel");
  }

  return d;
}
//------------------------------------------------------------------------------
double Environment::DiurnalTemperatureRange(){

 double d = 0;
 if( _Realm=="terrestrial" ) {
   d = GetVariableFromDatasetNamed("TerrestrialDTR");
  }
  return d;
    
}
//------------------------------------------------------------------------------
double Environment::NPP( ) {
 double d = 0;
 if( _Realm=="marine" ) {
   d= GetVariableFromDatasetNamed( "MarineNPP");
 } else if( _Realm=="terrestrial" ) {
   d= GetVariableFromDatasetNamed("TerrestrialNPP");
 }
 return d;
}
//------------------------------------------------------------------------------
double Environment::Precipitation(){

 double d = Constants::cMissingValue;//currently no marine precip
 if( _Realm=="terrestrial" ) {
   d = GetVariableFromDatasetNamed("TerrestrialPre");
  }
  return d;
    
}
//-----------------------------------------------------------------------------
double Environment::TotalPrecip(){return _TotalPrecip;}
//------------------------------------------------------------------------------
double Environment::AnnualTemperature(){return _AnnualTemperature;}
//------------------------------------------------------------------------------
double Environment::SDTemperature(){return _SDTemperature;}
//------------------------------------------------------------------------------
double Environment::GetVariableFromDatasetNamed(std:: string s){

    double d = Constants::cMissingValue;
    d = DataLayerSet::Get( )->GetDataAtCellIndexFor( s, _cellIndex );
    if( d == Constants::cMissingValue ) {
      std::cout << "Warning Environment::GetVariableFromDatasetNamed- missing values in "<<s<<" field!!"<< std::endl;
    }
    return d;
}
//------------------------------------------------------------------------------
double Environment::Latitude(){
    return Parameters::Get()->GetUserLatitudeAtIndex(LatitudeIndex());
}
//------------------------------------------------------------------------------
double Environment::Longitude(){
    return Parameters::Get()->GetUserLongitudeAtIndex(LongitudeIndex());
}
//------------------------------------------------------------------------------
unsigned Environment::LatitudeIndex(){
    Types::DataIndicesPointer indices = Parameters::Get( )->GetDataIndicesFromCellIndex( _cellIndex );
    return indices->GetY( );
}
//------------------------------------------------------------------------------
unsigned Environment::LongitudeIndex(){
    Types::DataIndicesPointer indices = Parameters::Get( )->GetDataIndicesFromCellIndex( _cellIndex );
    return indices->GetX( );
}
//------------------------------------------------------------------------------
double Environment::cellSize(){
    return Parameters::Get()->GetGridCellSize();
}



