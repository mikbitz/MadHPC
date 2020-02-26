#include "TimeStep.h"
#include "Constants.h"
#include "Parameters.h"
Types::TimeStepPointer TimeStep::_Instance = NULL;

//------------------------------------------------------------------------------------------------------------
TimeStep::TimeStep() {

}
//------------------------------------------------------------------------------------------------------------
TimeStep::~TimeStep( ) {
}
//------------------------------------------------------------------------------------------------------------
Types::TimeStepPointer TimeStep::instance( ) {
    if( _Instance == NULL ) {
        _Instance = new TimeStep( );
    }
    return _Instance;
}
//------------------------------------------------------------------------------------------------------------
void TimeStep::Initialise( std::string TimeStepUnits, double TimeStepLength, double SimulationLength, unsigned startStep) {

    _SimulationLength= SimulationLength;
    _TimeStepLength= TimeStepLength;
    _TimeStepUnits = TimeStepUnits;
    
    assert(_TimeStepUnits=="second" || 
           _TimeStepUnits=="minute" || 
           _TimeStepUnits=="hour"   ||
           _TimeStepUnits=="day"    ||
           _TimeStepUnits=="month"  ||
           _TimeStepUnits=="year");

    _NumberOfDaysPer["second"]=1./24/3600;
    _NumberOfDaysPer["minute"]=1./24/60;
    _NumberOfDaysPer["hour"]=1./24;
    _NumberOfDaysPer["day"]=1;
    _NumberOfDaysPer["month"]=365/12.;//roughly speaking!
    _NumberOfDaysPer["year"]= 365;
    _MonthsPerTimeStep =_NumberOfDaysPer[_TimeStepUnits]/_NumberOfDaysPer["month"]*_TimeStepLength;
    _DaysPerTimeStep   =_NumberOfDaysPer[_TimeStepUnits]*_TimeStepLength;
    _YearsPerTimeStep  =_NumberOfDaysPer[_TimeStepUnits] / _NumberOfDaysPer["year"];
    
    _DayZero=Parameters::instance()->GetInitialDay();
    _BaseStep=Parameters::instance()->GetBurninSteps();
    _Day=_DayZero + startStep*_DaysPerTimeStep;
    _DataAccessDay=_DayZero + startStep*_DaysPerTimeStep;
    _DataRepeatDay=Parameters::instance()->GetRepeatDay();
    
    _CurrentYear=unsigned(_DataAccessDay)/_NumberOfDaysPer["year"];
    _CurrentMonth=unsigned(_DataAccessDay/_NumberOfDaysPer["month"]) % 12;

    double tm= _Day;
    for (unsigned i=0;i<unsigned(SimulationLength/TimeStepLength);i++){_TimeStepArray.push_back(tm);tm+=_DaysPerTimeStep;}
    

}
//------------------------------------------------------------------------------------------------------------
void TimeStep::SetTime( const unsigned& CurrentTimeStep ) {
    _TimeStep=CurrentTimeStep;

    if (CurrentTimeStep>0){
        //_BaseStep determine{s the number of burn-in steps 
        if (_BaseStep>0){
            _BaseStep--;
        }else{
            _Day += _DaysPerTimeStep;
            _DataAccessDay += _DaysPerTimeStep;
            if (_DataRepeatDay> 0 && _DataAccessDay  >= _DataRepeatDay)while(_DataAccessDay  >= _DataRepeatDay){_DataAccessDay-=_DataRepeatDay;}
            _CurrentYear=unsigned(_DataAccessDay)/_NumberOfDaysPer["year"];
            _CurrentMonth=unsigned(_DataAccessDay/_NumberOfDaysPer["month"]) % 12;//months run from 0 to 11
        }
    }
}
//------------------------------------------------------------------------------------------------------------
void TimeStep::SetDay( const double& Day ) {
     _Day = Day;
     _DataAccessDay = Day;
     if (_DataRepeatDay> 0 && _DataAccessDay  > _DataRepeatDay)while(_DataAccessDay  >= _DataRepeatDay){_DataAccessDay-=_DataRepeatDay;}
}
//------------------------------------------------------------------------------------------------------------
double TimeStep::CurrentDay(){return _Day;}
//------------------------------------------------------------------------------------------------------------
double TimeStep::CurrentDataAccessDay(){return _DataAccessDay;}
//------------------------------------------------------------------------------------------------------------
unsigned TimeStep::CurrentMonth() {
        return _CurrentMonth;
}
//------------------------------------------------------------------------------------------------------------
unsigned TimeStep::CurrentYear(){return _CurrentYear;}
