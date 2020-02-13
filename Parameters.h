#ifndef PARAMETERS
#define PARAMETERS

#include "Types.h"

#include "repast_hpc/Properties.h"

class Parameters {
public:
    ~Parameters( );
    static Types::ParametersPointer instance( );

    bool Initialise( repast::Properties& );
    // User defined parameters
    std::string GetRootDataDirectory( ) const;
    std::string GetTimeStepUnits( ) const;
    unsigned GetLengthOfSimulationInYears( ) const;
    int GetMinimumLatitude( ) const;
    int GetMaximumLatitude( ) const;
    int GetMinimumLongitude( ) const;
    int GetMaximumLongitude( ) const;
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

    void SetMinimumLongitude( const int& );
    void SetMaximumLongitude( const int& );
    void SetMinimumLatitude( const int& );
    void SetMaximumLatitude( const int& );
    void SetGridCellSize( const double& );

    std::string GetDataDecriptorFileName();
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

    unsigned GetDataIndexOfMinimumLongitude( ) const;
    unsigned GetDataIndexOfMaximumLongitude( ) const;
    unsigned GetDataIndexOfMinimumLatitude( ) const;
    unsigned GetDataIndexOfMaximumLatitude( ) const;
    unsigned GetLengthLongitudeArray( ) const;
    unsigned GetLengthLatitudeArray( ) const;

    float GetLongitudeAtIndex( const unsigned& ) const;
    float GetLatitudeAtIndex( const unsigned& ) const;

    unsigned* GetTimeStepArray( ) const;
    unsigned* GetMonthlyTimeStepArray( ) const;
    unsigned* GetAnnualTimeStepArray( ) const;
    double GetInitialDay() const;
    double GetRepeatDay() const;

    float* GetLongitudeArray( ) const;
    float* GetLatitudeArray( ) const;
    bool GetNoLongitudeWrap( ) const;
    bool GetSpatialInterpolation( ) const;


private:
    Parameters( );
    void SetUpTimeStep(repast::Properties& props);
    void CalculateLonLat( );

    static Types::ParametersPointer _Instance;

    // User defined parameters
    std::string _RootDataDirectory;
    std::string _TimeStepUnits;
    float _TimeStepLength;
    int _MinimumLongitude;
    int _MaximumLongitude;
    int _MinimumLatitude;
    int _MaximumLatitude;
    double _GridCellSize;

    float _ExtinctionThreshold;
    unsigned _MaximumNumberOfCohorts;
    float _PlanktonSizeThreshold;
    bool _DrawRandomly;

    std::string _HumanNPPScenarioType;
    double _HumanNPPExtractionScale;
    double _HumanNPPScenarioDuration;

    double _InitialDay;
    double _DataRepeatDay;
    unsigned _BurninSteps;
    unsigned _ImpactSteps;
    unsigned _RecoverySteps;
    bool _noLongitudeWrap;
    bool _interpolateInputData;
    // Calculated parameters

    unsigned _NumberOfGridCells;

    unsigned _LengthLongitudeArray;
    unsigned _LengthLatitudeArray;

    float* _LongitudeArray;
    float* _LatitudeArray;
    
    std::string _DataDescriptorFileName;

};

#endif

