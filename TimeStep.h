#ifndef TIMESTEP
#define	TIMESTEP

#include "Types.h"
/**
\file TimeStep.h
\brief This is the TimeStep class header file.
*/

/**
\brief Timestep class largely for access to data files, since the data times might not match the model timestep points 

For the present time units are assumed to be days and for the convenience of the plant model one data value per month


@author Mike Bithell <mb425@cam.ac.uk>
*/
class TimeStep {
public:
    ~TimeStep( );

    static Types::TimeStepPointer instance( );  
    void SetTime( const unsigned& );
    void SetDay( const double& ); 
    void SetMonthly( const unsigned& );
    void Initialise( std::string, double,double,unsigned);
    double MonthsPerTimeStep( ) const {
       return _MonthsPerTimeStep;
    }
    double DaysPerTimeStep() const{
       return _DaysPerTimeStep;
    }
    double YearsPerTimeStep() const{
       return _YearsPerTimeStep;
    }
    std::vector<float> TimeStepArray(){return _TimeStepArray;}
    std::string TimeStepUnits(){return _TimeStepUnits;}
    double TimeStepLength(){return _TimeStepLength;}
    double SimulationLength(){return _SimulationLength;}
    unsigned TotalSimulationSteps(){return _TimeStepArray.size();}
    double CurrentDay();
    double CurrentDataAccessDay();
    unsigned CurrentMonth();
    unsigned CurrentYear();

private:
    TimeStep( );
    
    static Types::TimeStepPointer _Instance;
    
    unsigned _TimeStep;
    unsigned _CurrentMonth;
    unsigned _CurrentYear;
    std::map<std::string,double>_NumberOfDaysPer;
    double _Day;
    double _DataAccessDay;
    double _DataRepeatDay;
    double _BaseStep;
    double _DaysPerTimeStep; 
    double _DayZero;
    double _MonthsPerTimeStep;
    double _YearsPerTimeStep;
    std::string _TimeStepUnits;
    double _TimeStepLength;
    std::vector<float> _TimeStepArray;
    double _SimulationLength;

    
};

#endif

