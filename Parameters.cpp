#include "Parameters.h"

#include "Constants.h"
#include "Convertor.h"
#include "Maths.h"
#include "Processor.h"
#include "DataCoords.h"
#include "DataIndices.h"
#include "repast_hpc/Utilities.h"
#include <iostream>

Types::ParametersPointer Parameters::mThis = NULL;

Types::ParametersPointer Parameters::Get( ) {
    if( mThis == NULL ) {
        mThis = new Parameters( );
    }
    return mThis;
}

Parameters::~Parameters( ) {

    delete[ ] mMonthlyTimeStepArray;
    delete[ ] mAnnualTimeStepArray;
    delete[ ] mDataLongitudeArray;
    delete[ ] mDataLatitudeArray;
    delete[ ] mUserLongitudeArray;
    delete[ ] mUserLatitudeArray;

    if( mThis != NULL ) {
        delete mThis;
    }
}

Parameters::Parameters( ) {
}
bool Parameters::Initialise( repast::Properties& props ) {
    bool success=false;
    try {
        SetTimeStepUnits             (                    props.getProperty("simulation.TimeStepUnits"));
        SetLengthOfSimulation        (repast::strToDouble(props.getProperty("simulation.LengthOfSimulation")));
        SetTimeStepLength            (                   (props.getProperty("simulation.TimeStepLength")));
        SetUserMinimumLongitude      (repast::strToDouble(props.getProperty("simulation.minimumLongitude")));
        SetUserMaximumLongitude      (repast::strToDouble(props.getProperty("simulation.maximumLongitude")));
        SetUserMinimumLatitude       (repast::strToDouble(props.getProperty("simulation.minimumLatitude")));
        SetUserMaximumLatitude       (repast::strToDouble(props.getProperty("simulation.maximumLatitude")));
        SetGridCellSize              (repast::strToDouble(props.getProperty("simulation.GridCellSize")));
        SetExtinctionThreshold       (repast::strToDouble(props.getProperty("simulation.ExtinctionThreshold")));
        SetMaximumNumberOfCohorts    (repast::strToDouble(props.getProperty("simulation.MaximumNumberOfCohorts")));
        SetPlanktonSizeThreshold     (repast::strToDouble(props.getProperty("simulation.PlanktonSizeThreshold")));
        SetHumanNPPExtractionScale   (repast::strToDouble(props.getProperty("simulation.HumanNPPExtractionScale")));
        SetHumanNPPScenarioDuration  (repast::strToDouble(props.getProperty("simulation.HumanNPPScenarioDuration")));
        SetBurninSteps               (repast::strToDouble(props.getProperty("simulation.BurninSteps")));
        SetImpactSteps               (repast::strToDouble(props.getProperty("simulation.ImpactSteps")));
        SetRecoverySteps             (repast::strToDouble(props.getProperty("simulation.RecoverySteps")));
        SetDrawRandomly              (                    props.getProperty("simulation.DrawRandomly"));
        SetHumanNPPScenarioType      (                    props.getProperty("simulation.HumanNPPScenarioType"));
        SetRootDataDirectory         (                    props.getProperty("simulation.RootDataDirectory"));

        std::map<std::string,double>t;
        t["second"]=1./24/3600;
        t["minute"]=1./24/60;
        t["hour"]=1./24;
        t["day"]=1;
        t["month"]=30;
        t["year"]=12*30;

        //in Repast you can use stop.at onthe command line model.props to determine run-length of simulation 
        if (props.getProperty("stop.at").length()==0)
            props.putProperty("stop.at",unsigned(mLengthOfSimulation/mTimeStepLength));
        else
            SetLengthOfSimulation(repast::strToInt(props.getProperty("stop.at")));
        
        //legacy values needed because data reading all assumes timesteps must be months
        SetLengthOfSimulationInYears (mLengthOfSimulation/t["year"] *t[mTimeStepUnits]);
        SetLengthOfSimulationInMonths(mLengthOfSimulation/t["month"]*t[mTimeStepUnits]);

        SetMonthsPerTimeStep(t[mTimeStepUnits]/t["month"]*mTimeStepLength);
        
        CalculateParameters( );
        //throw "Parameters.cpp: Something bad happened when trying to read or calculate input parameters: check the model.props file? Can't continue. Exiting...";
        //make sure the dimensions are consistent with the environmental data files
        props.putProperty("min.x",0);
        props.putProperty("min.y",0);
        props.putProperty("max.x",GetLengthUserLongitudeArray( )-1);
        props.putProperty("max.y",GetLengthUserLatitudeArray( )-1);
        //only wrap in longitude if the domain size is exactly 360 degrees! Not super robust...
        props.putProperty("noLongitudeWrap",((mUserMaximumLongitude-mUserMinimumLongitude)%360 !=0));

  
        success = true;
    }catch(char* e){
        std::cout<<e<<std::endl;
        exit(1);
    }
    return success;
}


void Parameters::CalculateParameters( ) {


    mMonthlyTimeStepArray = new unsigned[ mLengthOfSimulationInMonths ];
    for( unsigned monthIndex = 0; monthIndex < mLengthOfSimulationInMonths; ++monthIndex ) {
        mMonthlyTimeStepArray[ monthIndex ] = monthIndex;
    }

    mAnnualTimeStepArray = new unsigned[ mLengthOfSimulationInYears ];
    for( unsigned yearIndex = 0; yearIndex < mLengthOfSimulationInYears; ++yearIndex ) {
        mAnnualTimeStepArray[ yearIndex ] = yearIndex;
    }

    // Calculate spatial parameters
    mLengthDataLongitudeArray = 360 / mGridCellSize;

    //there are some assumptions here about lat. and long ranges - this code may not work for grid boxes spanning the dateline.
    mDataLongitudeArray = new float[ mLengthDataLongitudeArray ];
    for( unsigned longitudeIndex = 0; longitudeIndex < mLengthDataLongitudeArray; ++longitudeIndex ) {
        mDataLongitudeArray[ longitudeIndex ] = ( -180 + ( ( float )mGridCellSize / 2 ) ) + ( longitudeIndex * ( float )mGridCellSize );
    }
    mLengthDataLatitudeArray = 180 / mGridCellSize;
    mDataLatitudeArray = new float[ mLengthDataLatitudeArray ];
    for( unsigned latitudeIndex = 0; latitudeIndex < mLengthDataLatitudeArray; ++latitudeIndex ) {
        mDataLatitudeArray[ latitudeIndex ] = ( -90 + ( ( float )mGridCellSize / 2 ) ) + ( latitudeIndex * ( float )mGridCellSize );
    }

    mDataIndexOfUserMinimumLongitude = Processor::Get( )->CalculateArrayIndexOfValue( mDataLongitudeArray, mLengthDataLongitudeArray, mUserMinimumLongitude );
    mDataIndexOfUserMaximumLongitude = Processor::Get( )->CalculateArrayIndexOfValue( mDataLongitudeArray, mLengthDataLongitudeArray, mUserMaximumLongitude );
    mLengthUserLongitudeArray = ( mDataIndexOfUserMaximumLongitude - mDataIndexOfUserMinimumLongitude ) + 1;

    mUserLongitudeArray = new float[ mLengthUserLongitudeArray ];
    for( unsigned userLongitudeIndex = 0; userLongitudeIndex < mLengthUserLongitudeArray; ++userLongitudeIndex ) {
        mUserLongitudeArray[ userLongitudeIndex ] = mDataLongitudeArray[ userLongitudeIndex + mDataIndexOfUserMinimumLongitude ];
    }

    mDataIndexOfUserMinimumLatitude = Processor::Get( )->CalculateArrayIndexOfValue( mDataLatitudeArray, mLengthDataLatitudeArray, mUserMinimumLatitude );
    mDataIndexOfUserMaximumLatitude = Processor::Get( )->CalculateArrayIndexOfValue( mDataLatitudeArray, mLengthDataLatitudeArray, mUserMaximumLatitude );
    mLengthUserLatitudeArray = ( mDataIndexOfUserMaximumLatitude - mDataIndexOfUserMinimumLatitude ) + 1;

    mUserLatitudeArray = new float[ mLengthUserLatitudeArray ];
    for( unsigned userLatitudeIndex = 0; userLatitudeIndex < mLengthUserLatitudeArray; ++userLatitudeIndex ) {
        mUserLatitudeArray[ userLatitudeIndex ] = mDataLatitudeArray[ userLatitudeIndex + mDataIndexOfUserMinimumLatitude ];
    }

    mNumberOfGridCells = mLengthUserLongitudeArray * mLengthUserLatitudeArray;
    mSizeOfMonthlyGridDatum = mNumberOfGridCells * mLengthOfSimulationInMonths;
    mSizeOfAnnualGridDatum = mNumberOfGridCells * mLengthOfSimulationInYears;

    unsigned cellIndex = 0;
    mCoordsIndicesLookup.resize( mNumberOfGridCells );
    for( unsigned latitudeIndex = 0; latitudeIndex < mLengthUserLatitudeArray; ++latitudeIndex ) {
        for( unsigned longitudeIndex = 0; longitudeIndex < mLengthUserLongitudeArray; ++longitudeIndex ) {

            float longitude = mUserLongitudeArray[ longitudeIndex ];
            float latitude = mUserLatitudeArray[ latitudeIndex ];

            Types::DataCoordsPointer coords = new DataCoords( longitude, latitude );
            Types::DataIndicesPointer indices = new DataIndices( longitudeIndex, latitudeIndex );

            mCoordsIndicesLookup[ cellIndex ] = std::make_pair( coords, indices );

            cellIndex += 1;
        }
    }
}

void Parameters::SetTimeStepLength( const std::string& lengthString){
    if (lengthString.length()==0)
        mTimeStepLength=1;
    else
       mTimeStepLength=repast::strToDouble(lengthString); 
}
void Parameters::SetMonthsPerTimeStep(float monthsPerTimeStep){
    mMonthsPerTimeStep=monthsPerTimeStep;
}
float Parameters::MonthsPerTimeStep( ) const {
    return mMonthsPerTimeStep;
}
float Parameters::DaysPerTimeStep() const{
    return 30.*mMonthsPerTimeStep;
}
std::string Parameters::GetRootDataDirectory( ) const {
    return mRootDataDirectory;
}

std::string Parameters::GetTimeStepUnits( ) const {
    return mTimeStepUnits;
}

unsigned Parameters::GetLengthOfSimulationInYears( ) const {
    return mLengthOfSimulationInYears;
}

int Parameters::GetUserMinimumLongitude( ) const {
    return mUserMinimumLongitude;
}

int Parameters::GetUserMaximumLongitude( ) const {
    return mUserMaximumLongitude;
}

int Parameters::GetUserMinimumLatitude( ) const {
    return mUserMinimumLatitude;
}

int Parameters::GetUserMaximumLatitude( ) const {
    return mUserMaximumLatitude;
}

double Parameters::GetGridCellSize( ) const {
    return mGridCellSize;
}

float Parameters::GetExtinctionThreshold( ) const {
    return mExtinctionThreshold;
}

unsigned Parameters::GetMaximumNumberOfCohorts( ) const {
    return mMaximumNumberOfCohorts;
}

float Parameters::GetPlanktonSizeThreshold( ) const {
    return mPlanktonSizeThreshold;
}

bool Parameters::GetDrawRandomly( ) const {
    return mDrawRandomly;
}


std::string Parameters::GetHumanNPPScenarioType( ) const {
    return mHumanNPPScenarioType;
}
double Parameters::GetHumanNPPExtractionScale( ) const{
    return mHumanNPPExtractionScale;
}
double Parameters::GetHumanNPPScenarioDuration( ) const{
    return mHumanNPPScenarioDuration;
}
unsigned Parameters::GetBurninSteps( ) const{
    return mBurninSteps;
}
unsigned Parameters::GetImpactSteps( ) const{
    return mImpactSteps;
}
unsigned Parameters::GetRecoverySteps( ) const{
    return mRecoverySteps;
}

void Parameters::SetRootDataDirectory( const std::string& rootDataDirectory ) {
    mRootDataDirectory = rootDataDirectory;
}

void Parameters::SetTimeStepUnits( const std::string& timeStepUnits ) {
    mTimeStepUnits = timeStepUnits;
}
void Parameters::SetLengthOfSimulation( const unsigned& lengthOfSimulation ) {
    mLengthOfSimulation = lengthOfSimulation;
}
void Parameters::SetLengthOfSimulationInYears( const unsigned& lengthOfSimulationInYears ) {
    mLengthOfSimulationInYears = lengthOfSimulationInYears;
}
void Parameters::SetLengthOfSimulationInMonths( const unsigned& lengthOfSimulationInMonths ) {
    mLengthOfSimulationInMonths = lengthOfSimulationInMonths;
}

void Parameters::SetUserMinimumLongitude( const int& userMinimumLongitude ) {
    mUserMinimumLongitude = userMinimumLongitude;
}

void Parameters::SetUserMaximumLongitude( const int& userMaximumLongitude ) {
    mUserMaximumLongitude = userMaximumLongitude;
}

void Parameters::SetUserMinimumLatitude( const int& userMinimumLatitude ) {
    mUserMinimumLatitude = userMinimumLatitude;
}

void Parameters::SetUserMaximumLatitude( const int& userMaximumLatitude ) {
    mUserMaximumLatitude = userMaximumLatitude;
}

void Parameters::SetGridCellSize( const double& gridCellSize ) {
    mGridCellSize = gridCellSize;
}

void Parameters::SetExtinctionThreshold( const float& extinctionThreshold ) {
    mExtinctionThreshold = extinctionThreshold;
}

void Parameters::SetMaximumNumberOfCohorts( const unsigned& maximumNumberOfCohorts ) {
    mMaximumNumberOfCohorts = maximumNumberOfCohorts;
}

void Parameters::SetPlanktonSizeThreshold( const float& planktonSizeThreshold ) {
    mPlanktonSizeThreshold = planktonSizeThreshold;
}

void Parameters::SetDrawRandomly( const std::string& drawRandomlyString ) {
    if( drawRandomlyString == "yes" )
        mDrawRandomly = true;
    else
        mDrawRandomly = false;
}

void Parameters::SetHumanNPPScenarioType(const std::string& humanNPPScenarioType){
    mHumanNPPScenarioType=humanNPPScenarioType;
}
void Parameters::SetHumanNPPExtractionScale(const double& humanNPPExtractionScale ){
    mHumanNPPExtractionScale=humanNPPExtractionScale;
}
void Parameters::SetHumanNPPScenarioDuration(const double & humanNPPScenarioDuration){
    mHumanNPPScenarioDuration=humanNPPScenarioDuration;
}
void Parameters::SetBurninSteps(const unsigned& burninSteps){
    mBurninSteps=burninSteps;
}
void Parameters::SetImpactSteps(const unsigned& impactSteps){
    mImpactSteps=impactSteps;
}
void Parameters::SetRecoverySteps(const unsigned& recoverySteps){
    mRecoverySteps=recoverySteps;
}

unsigned Parameters::GetNumberOfGridCells( ) const {
    return mNumberOfGridCells;
}

unsigned Parameters::GetLengthOfSimulationInMonths( ) const {
    return mLengthOfSimulationInMonths;
}

unsigned Parameters::GetLengthDataLongitudeArray( ) const {
    return mLengthDataLongitudeArray;
}

unsigned Parameters::GetLengthDataLatitudeArray( ) const {
    return mLengthDataLatitudeArray;
}

unsigned Parameters::GetDataIndexOfUserMinimumLongitude( ) const {
    return mDataIndexOfUserMinimumLongitude;
}

unsigned Parameters::GetDataIndexOfUserMaximumLongitude( ) const {
    return mDataIndexOfUserMaximumLongitude;
}

unsigned Parameters::GetDataIndexOfUserMinimumLatitude( ) const {
    return mDataIndexOfUserMinimumLatitude;
}

unsigned Parameters::GetDataIndexOfUserMaximumLatitude( ) const {
    return mDataIndexOfUserMaximumLatitude;
}

unsigned Parameters::GetLengthUserLongitudeArray( ) const {
    return mLengthUserLongitudeArray;
}

unsigned Parameters::GetLengthUserLatitudeArray( ) const {
    return mLengthUserLatitudeArray;
}

unsigned Parameters::GetSizeOfMonthlyGridDatum( ) const {
    return mSizeOfMonthlyGridDatum;
}

unsigned Parameters::GetSizeOfAnnualGridDatum( ) const {
    return mSizeOfAnnualGridDatum;
}

float Parameters::GetDataLongitudeAtIndex( const unsigned& index ) const {
    return mDataLongitudeArray[ index ];
}

float Parameters::GetDataLatitudeAtIndex( const unsigned& index ) const {
    return mDataLatitudeArray[ index ];
}

float Parameters::GetUserLongitudeAtIndex( const unsigned& index ) const {
    return mUserLongitudeArray[ index ];
}

float Parameters::GetUserLatitudeAtIndex( const unsigned& index ) const {
    return mUserLatitudeArray[ index ];
}

//float* Parameters::GetDataLongitudeArray( ) const {
//    return mDataLongitudeArray;
//}
//
//float* Parameters::GetDataLatitudeArray( ) const {
//    return mDataLatitudeArray;
//}

unsigned* Parameters::GetMonthlyTimeStepArray( ) const {
    return mMonthlyTimeStepArray;
}

unsigned* Parameters::GetAnnualTimeStepArray( ) const {
    return mAnnualTimeStepArray;
}

float* Parameters::GetUserLongitudeArray( ) const {
    return mUserLongitudeArray;
}

float* Parameters::GetUserLatitudeArray( ) const {
    return mUserLatitudeArray;
}

int Parameters::GetCellIndexFromDataIndices( const unsigned& longitudeIndex, const unsigned& latitudeIndex ) const {

    int cellIndex = Constants::cMissingValue;
    for( unsigned index = 0; index < mNumberOfGridCells; ++index ) {
        Types::DataIndicesPointer indices = mCoordsIndicesLookup[ index ].second;

        if( indices->GetX( ) == longitudeIndex && indices->GetY( ) == latitudeIndex ) {
            cellIndex = index;
            break;
        }
    }
    return cellIndex;
}

Types::DataCoordsPointer Parameters::GetDataCoordsFromCellIndex( const unsigned& cellIndex ) const {
    return mCoordsIndicesLookup[ cellIndex ].first;
}

Types::DataIndicesPointer Parameters::GetDataIndicesFromCellIndex( const unsigned& cellIndex ) const {
    return mCoordsIndicesLookup[ cellIndex ].second;
}
