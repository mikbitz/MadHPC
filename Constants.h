#ifndef CONSTANTS
#define CONSTANTS

#include <string>

namespace Constants {
    
    enum eOutputControlParametersMetadata {
        eDatumName,
        eDatumType,
        eTimeUnit,
        eDataUnit
    };

    enum eVariableTypes {
        eLongitude,
        eLatitude,
        eTime,
        eDepth,
        eOther
    };

    const std::string cLongitudeVariableNames[ ] = { "lon", "long", "longitude", "x" };
    const std::string cLatitudeVariableNames[ ] = { "lat", "lats", "latitude", "y" };
    const std::string cDepthVariableNames[ ] = { "dep", "depth", "z" };
    const std::string cTimeVariableNames[ ] = { "time", "month", "t" };

    const std::string cBasicDatumTypeName = "basic";
    const std::string cGridDatumTypeName = "grid";

    const std::string cOutputBaseDirectory = "./output/";
    const std::string cDataSetNameFormat = "%Y-%m-%d_%H-%M-%S";
    const std::string cCompleteDateFormat = "%c";
    const std::string cUnitsString = "units";
    const std::string cTimeVariableUnits = "months";

    ///////////////////////////////////////////////////////////////////////////
    const std::string cAnnualBasicOutputsFileName = "AnnualBasicOutputs.nc";
    const std::string cMonthlyBasicOutputsFileName = "MonthlyBasicOutputs.nc";
    const std::string cAnnualGridOutputsFileName = "AnnualGridOutputs.nc";
    const std::string cMonthlyGridOutputsFileName = "MonthlyGridOutputs.nc";
    const std::string cMonthlyTimeUnitName = "month";
    const std::string cAnnualTimeUnitName = "year";
    const std::string cLongitudeVariableUnit = "degrees east";
    const std::string cLatitudeVariableUnit = "degrees north";
    ///////////////////////////////////////////////////////////////////////////

    const int cMissingValue = -9999;
    const int cOutputFolderPermissions = 0775;
    const int cDateTimeBufferSize = 25;

    const char cDataDelimiterValue = ',';
    const char cTextFileCommentCharacter = '#';
    const char cFolderDelimiter = '/';
    const char cWhiteSpaceCharacter = ' ';
    const double cMonth = 1;//assumes value for timestep is 1 month!
    const double cDay = 30;// should be cMonth*365./12.; but the original model assumes 30 days in a month!
    const double cYear = cMonth*12.;
}

#endif
