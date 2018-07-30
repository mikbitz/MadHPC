#ifndef ENVIRONMENT
#define ENVIRONMENT
#include "Constants.h"
#include "Types.h"
#include <vector>
class Environment {
public:
	Environment();
    Environment(int,int);
    void SetRealm( );


    std::string _Realm;


	virtual ~Environment() {}
	
	void merge(int&);
    double Temperature();
    double DiurnalTemperatureRange();
    double Width();
    double Height();
    double Area();
    double Breeding_Season();
    double AnnualTemperature();
    double SDTemperature();
    double NPP();
    double Precipitation();
    double TotalPrecip();
    double TotalAET();
    double AET();
    double FractionYearFrost();
    double FractionYearFire();
    double Seasonality();
    double ExpTDevWeight();
    double TerrestrialHANPP();
    double uVel();
    double vVel();
    unsigned index(){return _cellIndex;}
    double cellSize();
    double Latitude();
    double Longitude();
    unsigned LatitudeIndex();
    unsigned LongitudeIndex();
    unsigned _cellIndex;
    void addToOrganicPool(double d){_OrganicPool+=d;}
    void addToRespiratoryCO2Pool(double d){_RespiratoryCO2Pool+=d;}
    double organicPool(){return _OrganicPool;};
    double respiratoryPool(){return _RespiratoryCO2Pool;};
    void zeroPools(){    _OrganicPool=0; _RespiratoryCO2Pool=0;}
private:
    double GetVariableFromDatasetNamed(std:: string s);
    void   SetTotalPrecip();
    void   SetAVGSDTemp();
    void   SetNPPSeasonality();
    void   SetBreeding();
    void   SetFrostandFire();
    double _OrganicPool;
    double _RespiratoryCO2Pool;
    std::vector< double > _exptdev;
    std::vector< double > _Seasonality;
    std::vector< double > _Breeding_Season;
    std::vector< double >  _AET;
    double _Width;
    double _Height;
    double _Area;
    double _AnnualTemperature;
    double _SDTemperature;
    double _TotalPrecip;
    double _TotalAET;

    double _FractionYearFrost;
    double _FractionYearFire;

    

};
#endif


