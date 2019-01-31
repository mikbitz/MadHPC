/*
 *
 * Cohort.h
 *
 *  Created on: May 7, 2018
 *      Author: Mike Bithell
 * 
 *  Derived from Original C# code by
 * Emergent Global Patterns of Ecosystem Structure and Function from a Mechanistic General Ecosystem Model , Michael B. J. Harfoot , Tim Newbold, Derek P. Tittensor,  Stephen Emmott, Jon Hutton, Vassily Lyutsarev, Matthew J. Smith, Jörn P. W. Scharlemann, Drew W. Purves PLOS Biology 12(4) 2014 e1001841 https://doi.org/10.1371/journal.pbio.1001841
 */

#ifndef COHORT_H_
#define COHORT_H_
#include "agent.h"
#include "repast_hpc/AgentId.h"
#include "repast_hpc/SharedContext.h"
#include "repast_hpc/SharedDiscreteSpace.h"
#include "AgentPackage.h"
#include "Environment.h"
#include "randomizer.h"

class MadModel;

class Cohort: public MadAgent  {

public:

    
    unsigned _FunctionalGroupIndex;   

    double _JuvenileMass;     
    double _AdultMass;               
    double _IndividualBodyMass;
    double _MaximumAchievedBodyMass;
    double _IndividualReproductivePotentialMass;
    double _MinimumMass;
    double _MaximumMass;
    
    double _CohortAbundance;
    unsigned _BirthTimeStep;            
    unsigned _MaturityTimeStep;            
    double _LogOptimalPreyBodySizeRatio; 
     
    bool _Merged;
    //bool _alive;
    //bool _moved;

  
	bool _Heterotroph;   
    bool _Autotroph; 
    bool _Endotherm; 
    bool _Ectotherm;
    
    std::string _Realm;      

    bool _Iteroparous;
    bool _Semelparous;
    bool _Herbivore;
    bool _Carnivore;
    bool _Omnivore;
    bool _IsPlanktonic;
    bool _IsFilterFeeder;

    
    double _ProportionSuitableTimeActive;
    
    bool _IsMature;
    
    double _AssimilationEfficiency_H;
    double _AssimilationEfficiency_C;

    double _ProportionTimeActive;

    double _PotentialAbundanceEaten;
    unsigned _CurrentTimeStep;
    
    //herbivores

    static const double _edibleFractionMarine           ;// =1.0;
    static const double _AttackRateExponentMarine       ;// =2.0;
    static const double _HandlingTimeExponentMarine     ;// =0.7;
    static const double _HandlingTimeScalarMarine       ;// =0.7;
    static const double _edibleFractionTerrestrial      ;// =0.1;
    static const double _AttackRateExponentTerrestrial  ;// =2.0;
    static const double _HandlingTimeExponentTerrestrial;// =0.7;
    static const double _HandlingTimeScalarTerrestrial  ;// =0.7;
    static const double _HerbivoryRateMassExponent      ;// =1.0;
    static const double _HerbivoryRateConstant          ;// =1.0e-11;
    static const double _ReferenceMass                  ;// =1.0;
 
    //Carnivores

    static const double _HandlingTimeScalar_C               ;//= 0.5;
    static const double _HandlingTimeExponent_C             ;//= 0.7;
    static const double _SearchRateConstant                 ;//= 1E-6;
    static const double _FeedingPreferenceStandardDeviation ;//= 0.7;
    static const double _NumberOfBins                       ;//= 12;

    //Dispersal constants
    static const double _DispersalSpeedBodyMassScalar       ;//= 0.0278;
    static const double _DispersalSpeedBodyMassExponent     ;//= 0.48;
    //Advective dispersal
    static const double _HorizontalDiffusivity                       ;//= 100;
    static const double _AdvectiveModelTimeStepLengthHours           ;//= 18;
    static const double _HorizontalDiffusivityKmSqPerADTimeStep      ;//= _HorizontalDiffusivity / ( 1000 * 1000 ) * 60 * 60 * AdvectiveModelTimeStepLengthHours;
    static const double _AdvectionTimeStepsPerModelTimeStep          ;//= Constants::cDay * 24 / _AdvectiveModelTimeStepLengthHours;
    static const double _VelocityUnitConversion                      ;//= 60 * 60 * 24 * Constants::cDay * Constants::cMonth  / 1000;
    //Responsive dispersal
    static const double _DensityThresholdScaling                     ;//= 50000;
    static const double _StarvationDispersalBodyMassThreshold        ;//=0.8;

    //ectotherms
    static const double _TerrestrialWarmingToleranceIntercept ;//= 6.61;
    static const double _TerrestrialWarmingToleranceSlope     ;//= 1.6;
    static const double _TerrestrialTSMIntercept              ;//= 1.51;
    static const double _TerrestrialTSMSlope                  ;//= 1.53;

    //metabolism
    //Ectotherms
    static const double _MetabolismMassExponentEcto            ;//= 0.88;
    static const double _NormalizationConstantEcto             ;//= 148984000000; //1.4898373851E+11;
    static const double _ActivationEnergyEcto                  ;//= 0.69; // includes endotherms in hibernation and torpor
    static const double _NormalizationConstantBMR              ;//= 41918272883; //exp( 20 )*60 * 60 * 24 / 1000;
    static const double _BasalMetabolismMassExponent           ;//= 0.69;
    static const double _EnergyScalarEcto                      ;//= 0.036697248; //1 / 27.25;
    //Endotherrms
    static const double _MetabolismMassExponentEndo           ;//= 0.7;
    static const double _NormalizationConstantEndo            ;//= 9.0809083973E+11;
    static const double _ActivationEnergyEndo                 ;//= 0.69; // includes endotherms in hibernation and torpor
    static const double _EnergyScalarEndo                     ;//= 1 / 27.25;

    static const double _BoltzmannConstant                    ;//= 8.617e-5;
    static const double _TemperatureUnitsConvert              ;// = 273.0;

    //reproduction
    static const double _MassRatioThreshold                   ;//= 1.5;
    static const double _MassEvolutionProbabilityThreshold    ;//= 0.95;
    static const double _MassEvolutionStandardDeviation       ;//= 0.05;
    static const double _SemelparityAdultMassAllocation       ;//= 0.5;
    //mortality
    static const double _MortalityRateBackground              ;//= 0.001;
    static const double _MortalityRateMature                  ;//= 0.003;
    static const double _LogisticInflectionPoint              ;//= 0.6;
    static const double _LogisticScalingParameter             ;//= 0.05;
    static const double _MaximumStarvationRate                ;//= 1;
    
    //Pi!
    static const double _Pi ;//=acos(-1.);
    //area units conversion
    static const double _CellAreaToHectares;//=100;
    
    //temporary store for within timestep changes
    static std::map < std::string, std::map<std::string,double> > _Accounting;



public:

    static unsigned _NextID;
    Cohort* _newH;
    Cohort(repast::AgentId id): MadAgent(id), _Merged(false){_NextID++;_newH=NULL; _location={0,0};}
	Cohort(repast::AgentId id, const AgentPackage& package): MadAgent(id){PullThingsOutofPackage(package);_newH=NULL;}
    void set(int currentRank, const AgentPackage& package){_id.currentRank(currentRank);PullThingsOutofPackage(package);}
	void setup(unsigned,unsigned,Environment*,randomizer*);

	virtual ~Cohort() {}

	void step(Environment* ,vector<Cohort*>&,vector<Stock*>&,const unsigned);

    void metabolize(Environment*);
    void assignTimeActive(Environment*);
    void reproduce(Environment*);
    void eat(Environment*,vector<Cohort*>&,vector<Stock*>&);
    void moveIt(Environment*,MadModel*);
    void mort();
    void markForDeath();
    void applyEcology(Environment*);
    void updatePools(Environment*);
    void setupOffspring( Cohort* , double , double , double , double , unsigned  );
    void TryToDisperse(double,Environment*,MadModel* );
    void TryToDisperse(double,double,Environment*,MadModel* );
    vector<int> _location,_destination;
    void setLocation(int x, int y){_location={x,y};}

void PushThingsIntoPackage( AgentPackage& );
void PullThingsOutofPackage( const AgentPackage& );
void ResetAccounts();
};
#endif /* COHORT_H_ */
