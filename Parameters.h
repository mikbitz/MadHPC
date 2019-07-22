#ifndef PARAMETERS
#define PARAMETERS

#include "Types.h"

#include "repast_hpc/Properties.h"

class Parameters {
public:
    ~Parameters( );
    static Types::ParametersPointer Get( );

    bool Initialise( repast::Properties& );
    // User defined parameters
    std::string GetRootDataDirectory( ) const;
    std::string GetTimeStepUnits( ) const;
    unsigned GetLengthOfSimulationInYears( ) const;
    int GetUserMinimumLatitude( ) const;
    int GetUserMaximumLatitude( ) const;
    int GetUserMinimumLongitude( ) const;
    int GetUserMaximumLongitude( ) const;
    double GetGridCellSize( ) const;
    float GetExtinctionThreshold( ) const;
    unsigned GetMaximumNumberOfCohorts( ) const;
    float GetPlanktonSizeThreshold( ) const;
    bool GetDrawRandomly( ) const;

    std::string GetHumanNPPScenarioType( ) const;
    double GetHumanNPPExtractionScale( ) const;
    double GetHumanNPPScenarioDuration( ) const;
    unsigned GetBurninSteps( ) const;
    unsigned GetImpactSteps( ) const;
    unsigned GetRecoverySteps( ) const;
    float MonthsPerTimeStep() const;
    float DaysPerTimeStep() const;   
    void SetRootDataDirectory( const std::string& );
    void SetTimeStepUnits( const std::string& );
    void SetTimeStepLength( const std::string& );
    void SetMonthsPerTimeStep(float);
    void SetLengthOfSimulation( const unsigned& );
    void SetLengthOfSimulationInYears( const unsigned& );
    void SetLengthOfSimulationInMonths( const unsigned& );
    void SetUserMinimumLongitude( const int& );
    void SetUserMaximumLongitude( const int& );
    void SetUserMinimumLatitude( const int& );
    void SetUserMaximumLatitude( const int& );
    void SetGridCellSize( const double& );
    void SetExtinctionThreshold( const float& );
    void SetMaximumNumberOfCohorts( const unsigned& );
    void SetPlanktonSizeThreshold( const float& );
    void SetDrawRandomly( const std::string& );

    void SetHumanNPPScenarioType(const std::string& humanNPPScenarioType);
    void SetHumanNPPExtractionScale(const double& humanNPPExtractionScale );
    void SetHumanNPPScenarioDuration(const double & humanNPPScenarioDuration);
    void SetBurninSteps(const unsigned& burninSteps);
    void SetImpactSteps(const unsigned& impactSteps);
    void SetRecoverySteps(const unsigned& recoverySteps);
    
    // Calculated parameters
    unsigned GetNumberOfGridCells( ) const;
    unsigned GetLengthOfSimulationInMonths( ) const;
    unsigned GetLengthDataLongitudeArray( ) const;
    unsigned GetLengthDataLatitudeArray( ) const;
    unsigned GetDataIndexOfUserMinimumLongitude( ) const;
    unsigned GetDataIndexOfUserMaximumLongitude( ) const;
    unsigned GetDataIndexOfUserMinimumLatitude( ) const;
    unsigned GetDataIndexOfUserMaximumLatitude( ) const;
    unsigned GetLengthUserLongitudeArray( ) const;
    unsigned GetLengthUserLatitudeArray( ) const;

    unsigned GetSizeOfAnnualGridDatum( ) const;
    unsigned GetSizeOfMonthlyGridDatum( ) const;

    float GetDataLongitudeAtIndex( const unsigned& ) const;
    float GetDataLatitudeAtIndex( const unsigned& ) const;
    float GetUserLongitudeAtIndex( const unsigned& ) const;
    float GetUserLatitudeAtIndex( const unsigned& ) const;

    unsigned* GetTimeStepArray( ) const;
    unsigned* GetMonthlyTimeStepArray( ) const;
    unsigned* GetAnnualTimeStepArray( ) const;
    //float* GetDataLongitudeArray( ) const;
    //float* GetDataLatitudeArray( ) const;
    float* GetUserLongitudeArray( ) const;
    float* GetUserLatitudeArray( ) const;
    
    int GetCellIndexFromDataIndices( const unsigned&, const unsigned& ) const;
    Types::DataCoordsPointer GetDataCoordsFromCellIndex( const unsigned& ) const;
    Types::DataIndicesPointer GetDataIndicesFromCellIndex( const unsigned& ) const;

private:
    Parameters( );
    void CalculateParameters( );

    static Types::ParametersPointer mThis;

    // User defined parameters
    std::string mRootDataDirectory;
    std::string mTimeStepUnits;
    float mTimeStepLength;
    int mUserMinimumLongitude;
    int mUserMaximumLongitude;
    int mUserMinimumLatitude;
    int mUserMaximumLatitude;
    double mGridCellSize;
    float mExtinctionThreshold;
    unsigned mMaximumNumberOfCohorts;
    float mPlanktonSizeThreshold;
    bool mDrawRandomly;

    std::string mHumanNPPScenarioType;
    double mHumanNPPExtractionScale;
    double mHumanNPPScenarioDuration;
    float mMonthsPerTimeStep;
    unsigned mBurninSteps;
    unsigned mImpactSteps;
    unsigned mRecoverySteps;
    // Calculated parameters
    unsigned mLengthOfSimulationInMonths;
    unsigned mLengthOfSimulationInYears;
    unsigned mLengthOfSimulation;
    unsigned mNumberOfGridCells;
    unsigned mLengthDataLongitudeArray;
    unsigned mLengthDataLatitudeArray;
    unsigned mLengthUserLongitudeArray;
    unsigned mLengthUserLatitudeArray;
    unsigned mDataIndexOfUserMinimumLongitude;
    unsigned mDataIndexOfUserMaximumLongitude;
    unsigned mDataIndexOfUserMinimumLatitude;
    unsigned mDataIndexOfUserMaximumLatitude;
    unsigned mSizeOfMonthlyGridDatum;
    unsigned mSizeOfAnnualGridDatum;
    unsigned* mTimeStepArray;
    unsigned* mMonthlyTimeStepArray;
    unsigned* mAnnualTimeStepArray;
    float* mDataLongitudeArray;
    float* mDataLatitudeArray;
    float* mUserLongitudeArray;
    float* mUserLatitudeArray;

    Types::CoordsIndicesVector mCoordsIndicesLookup;
};

#endif

