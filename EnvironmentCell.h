/*
 *
 * EnvironmentCell.h
 *
 *  Created on: May 7, 2018
 *      Author: Mike Bithell
 * 
 *  Derived from Origianl C# code by
 * Emergent Global Patterns of Ecosystem Structure and Function from a Mechanistic General Ecosystem Model , Michael B. J. Harfoot , Tim Newbold, Derek P. Tittensor,  Stephen Emmott, Jon Hutton, Vassily Lyutsarev, Matthew J. Smith, Jörn P. W. Scharlemann, Drew W. Purves PLOS Biology 12(4) 2014 e1001841 https://doi.org/10.1371/journal.pbio.1001841
 */
#ifndef ENVIRONMENTCELL
#define ENVIRONMENTCELL
#include "Constants.h"
#include "Types.h"
#include <vector>
class EnvironmentCell {
public:
	EnvironmentCell();
    EnvironmentCell(int,int);
    void SetRealm( );


    std::string _Realm;


	virtual ~EnvironmentCell() {}
	
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
    double cellSize();
    double Latitude();
    double Longitude();
    unsigned _x,_y;
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


