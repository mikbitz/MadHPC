#include "Parameters.h"

#include "Constants.h"
#include "Convertor.h"
#include "repast_hpc/Utilities.h"
#include "TimeStep.h"
#include "Types.h"
#include <iostream>

Types::ParametersPointer Parameters::_Instance = NULL;

Types::ParametersPointer Parameters::instance( ) {
    if( _Instance == NULL ) {
        _Instance = new Parameters( );
    }
    return _Instance;
}

Parameters::~Parameters( ) {

    delete[ ] _LongitudeArray;
    delete[ ] _LatitudeArray;

    if( _Instance != NULL ) {
        delete _Instance;
    }
}

Parameters::Parameters( ) {
}
bool Parameters::Initialise( repast::Properties& props ) {
    bool success=false;
    try{
                
        SetMinimumLongitude      (repast::strToDouble(props.getProperty("simulation.minimumLongitude")));
        SetMaximumLongitude      (repast::strToDouble(props.getProperty("simulation.maximumLongitude")));
        SetMinimumLatitude       (repast::strToDouble(props.getProperty("simulation.minimumLatitude")));
        SetMaximumLatitude       (repast::strToDouble(props.getProperty("simulation.maximumLatitude")));
        SetGridCellSize              (repast::strToDouble(props.getProperty("simulation.GridCellSize")));
        SetExtinctionThreshold       (repast::strToDouble(props.getProperty("simulation.ExtinctionThreshold")));
        SetMaximumNumberOfCohorts    (repast::strToDouble(props.getProperty("simulation.MaximumNumberOfCohorts")));
        SetPlanktonSizeThreshold     (repast::strToDouble(props.getProperty("simulation.PlanktonSizeThreshold")));
        SetHumanNPPExtractionScale   (repast::strToDouble(props.getProperty("simulation.HumanNPPExtractionScale")));
        SetHumanNPPScenarioDuration  (repast::strToDouble(props.getProperty("simulation.HumanNPPScenarioDuration")));
        SetBurninSteps               (repast::strToInt(props.getProperty("simulation.BurninSteps")));
        SetImpactSteps               (repast::strToDouble(props.getProperty("simulation.ImpactSteps")));
        SetRecoverySteps             (repast::strToDouble(props.getProperty("simulation.RecoverySteps")));
        SetDrawRandomly              (                    props.getProperty("simulation.DrawRandomly"));
        SetHumanNPPScenarioType      (                    props.getProperty("simulation.HumanNPPScenarioType"));

        SetRootDataDirectory         (                    props.getProperty("input.RootNcDataDirectory"));

       _DataDescriptorFileName=props.getProperty("input.DataDirectory")+"/"+props.getProperty("input.EnvironmentVariablesFile");
       
        CalculateLonLat( );
        //throw "Parameters.cpp: Something bad happened when trying to read or calculate input parameters: check the model.props file? Can't continue. Exiting...";

        //make sure the dimensions are consistent with the environmental data files
        props.putProperty("min.x",0);
        props.putProperty("min.y",0);
        props.putProperty("max.x",GetLengthLongitudeArray( )-1);
        props.putProperty("max.y",GetLengthLatitudeArray( )-1);
        //only don't wrap in longitude if set to false here 
        std::string wrap=props.getProperty("simulation.LongitudeWrap");
        if (wrap=="")wrap="true";
        assert(wrap=="false" || wrap=="true");
        props.putProperty("noLongitudeWrap",(wrap=="false"));
        _noLongitudeWrap=(wrap=="false");

        std::string interpol=props.getProperty("simulation.interpolateInputData");
        if (interpol=="")interpol="true";
        assert(interpol=="false" || interpol=="true");
        props.putProperty("interpolateInputData",(interpol=="true"));
        _interpolateInputData=(interpol=="true");

        SetUpTimeStep(props);

        success = true;
    }catch(char* e){
        std::cout<<e<<std::endl;
        exit(1);
    }
    return success;
}

void Parameters::SetUpTimeStep(repast::Properties& props){

    unsigned startStep=0;
    std::string rstrt=props.getProperty("simulation.RestartStep");
    if (rstrt!="")startStep=repast::strToInt(rstrt);
    std::string instep=props.getProperty("simulation.InitialStep");
    if (startStep>0 && instep!="")startStep=repast::strToInt(instep);
    
    double TimeStepLength;
    std::string TimeStepLengthStr = props.getProperty("simulation.TimeStepLength");
    if  (TimeStepLengthStr.length()==0) TimeStepLength=1; 
    else TimeStepLength=repast::strToDouble(TimeStepLengthStr);
        
    double SimulationLength=repast::strToDouble(props.getProperty("simulation.LengthOfSimulation"));
                
    //in Repast you can use stop.at onthe command line model.props to determine run-length of simulation
    //scheduler will begin at startStep+1, so stop.at needs to be increased accordingly
    if (props.getProperty("stop.at").length()==0)
          props.putProperty("stop.at",unsigned(SimulationLength/TimeStepLength)+startStep);
    else
          assert( (repast::strToInt(props.getProperty("stop.at")) <= unsigned(SimulationLength/TimeStepLength)+startStep) );
    
    std::string InitialDayStr = props.getProperty("simulation.InitialDay");
    if   (InitialDayStr.length()==0) _InitialDay=0; 
    else _InitialDay=repast::strToDouble(InitialDayStr);
    
    std::string DayStr = props.getProperty("simulation.DataRepeatDay");    
    if   (DayStr.length()==0) _DataRepeatDay=0; 
    else _DataRepeatDay=repast::strToDouble(DayStr);
    
    TimeStep::instance()->Initialise(props.getProperty("simulation.TimeStepUnits"),TimeStepLength,SimulationLength,startStep);

}

void Parameters::CalculateLonLat( ) {

    _LengthLongitudeArray = ( _MaximumLongitude - _MinimumLongitude )/_GridCellSize + 1;

    _LongitudeArray = new float[ _LengthLongitudeArray ];
    for( unsigned userLongitudeIndex = 0; userLongitudeIndex < _LengthLongitudeArray; ++userLongitudeIndex ) {
        _LongitudeArray[ userLongitudeIndex ] = (_MinimumLongitude  ) + ( userLongitudeIndex * ( float )_GridCellSize );

    }

    _LengthLatitudeArray = ( _MaximumLatitude - _MinimumLatitude )/_GridCellSize + 1;

    _LatitudeArray = new float[ _LengthLatitudeArray ];
    for( unsigned userLatitudeIndex = 0; userLatitudeIndex < _LengthLatitudeArray; ++userLatitudeIndex ) {
        _LatitudeArray[ userLatitudeIndex ] = (_MinimumLatitude  ) + ( userLatitudeIndex * ( float )_GridCellSize );
    }

    _NumberOfGridCells = _LengthLongitudeArray * _LengthLatitudeArray;


}
double Parameters::GetInitialDay() const{
    return _InitialDay;
}
double Parameters::GetRepeatDay() const{
    return _DataRepeatDay;
}
std::string Parameters::GetDataDecriptorFileName(){
    return _DataDescriptorFileName;
}

std::string Parameters::GetRootDataDirectory( ) const {
    return _RootDataDirectory;
}

int Parameters::GetMinimumLongitude( ) const {
    return _MinimumLongitude;
}

int Parameters::GetMaximumLongitude( ) const {
    return _MaximumLongitude;
}

int Parameters::GetMinimumLatitude( ) const {
    return _MinimumLatitude;
}

int Parameters::GetMaximumLatitude( ) const {
    return _MaximumLatitude;
}

double Parameters::GetGridCellSize( ) const {
    return _GridCellSize;
}

float Parameters::GetExtinctionThreshold( ) const {
    return _ExtinctionThreshold;
}

unsigned Parameters::GetMaximumNumberOfCohorts( ) const {
    return _MaximumNumberOfCohorts;
}

float Parameters::GetPlanktonSizeThreshold( ) const {
    return _PlanktonSizeThreshold;
}

bool Parameters::GetDrawRandomly( ) const {
    return _DrawRandomly;
}

std::string Parameters::GetHumanNPPScenarioType( ) const {
    return _HumanNPPScenarioType;
}
double Parameters::GetHumanNPPExtractionScale( ) const{
    return _HumanNPPExtractionScale;
}
double Parameters::GetHumanNPPScenarioDuration( ) const{
    return _HumanNPPScenarioDuration;
}
unsigned Parameters::GetBurninSteps( ) const{

    return _BurninSteps;
}
unsigned Parameters::GetImpactSteps( ) const{
    return _ImpactSteps;
}
unsigned Parameters::GetRecoverySteps( ) const{
    return _RecoverySteps;
}

void Parameters::SetRootDataDirectory( const std::string& rootDataDirectory ) {
    _RootDataDirectory = rootDataDirectory;
}

void Parameters::SetMinimumLongitude( const int& userMinimumLongitude ) {
    _MinimumLongitude = userMinimumLongitude;
}

void Parameters::SetMaximumLongitude( const int& userMaximumLongitude ) {
    _MaximumLongitude = userMaximumLongitude;
}

void Parameters::SetMinimumLatitude( const int& userMinimumLatitude ) {
    _MinimumLatitude = userMinimumLatitude;
}

void Parameters::SetMaximumLatitude( const int& userMaximumLatitude ) {
    _MaximumLatitude = userMaximumLatitude;
}

void Parameters::SetGridCellSize( const double& gridCellSize ) {
    _GridCellSize = gridCellSize;
}

void Parameters::SetExtinctionThreshold( const float& extinctionThreshold ) {
    _ExtinctionThreshold = extinctionThreshold;
}

void Parameters::SetMaximumNumberOfCohorts( const unsigned& maximumNumberOfCohorts ) {
    _MaximumNumberOfCohorts = maximumNumberOfCohorts;
}

void Parameters::SetPlanktonSizeThreshold( const float& planktonSizeThreshold ) {
    _PlanktonSizeThreshold = planktonSizeThreshold;
}

void Parameters::SetDrawRandomly( const std::string& drawRandomlyString ) {
    if( drawRandomlyString == "yes" )
        _DrawRandomly = true;
    else
        _DrawRandomly = false;
}

void Parameters::SetHumanNPPScenarioType(const std::string& humanNPPScenarioType){
    _HumanNPPScenarioType=humanNPPScenarioType;
}
void Parameters::SetHumanNPPExtractionScale(const double& humanNPPExtractionScale ){
    _HumanNPPExtractionScale=humanNPPExtractionScale;
}
void Parameters::SetHumanNPPScenarioDuration(const double & humanNPPScenarioDuration){
    _HumanNPPScenarioDuration=humanNPPScenarioDuration;
}
void Parameters::SetBurninSteps(const unsigned& burninSteps){

    _BurninSteps=burninSteps;

}
void Parameters::SetImpactSteps(const unsigned& impactSteps){
    _ImpactSteps=impactSteps;
}
void Parameters::SetRecoverySteps(const unsigned& recoverySteps){
    _RecoverySteps=recoverySteps;
}

unsigned Parameters::GetNumberOfGridCells( ) const {
    return _NumberOfGridCells;
}

unsigned Parameters::GetLengthLongitudeArray( ) const {
    return _LengthLongitudeArray;
}

unsigned Parameters::GetLengthLatitudeArray( ) const {
    return _LengthLatitudeArray;
}

float Parameters::GetLongitudeAtIndex( const unsigned& index ) const {
    return _LongitudeArray[ index ];
}

float Parameters::GetLatitudeAtIndex( const unsigned& index ) const {
    return _LatitudeArray[ index ];
}

float* Parameters::GetLongitudeArray( ) const {
    return _LongitudeArray;
}

float* Parameters::GetLatitudeArray( ) const {
    return _LatitudeArray;
}
bool Parameters::GetNoLongitudeWrap( ) const {
    return _noLongitudeWrap;
}
bool Parameters::GetSpatialInterpolation( ) const {
    return _interpolateInputData;
}
