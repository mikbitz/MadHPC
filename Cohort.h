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
#include "repast_hpc/Properties.h"
#include "repast_hpc/Utilities.h"

#include "AgentPackage.h"
#include "Environment.h"
#include "randomizer.h"

class MadModel;

class Cohort: public MadAgent  {

public:

    int _sequencer;//used to co-ordinate ordering of cohort behaviour across threads
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
    
    static void setParameters(repast::Properties*);

    //herbivores

    static double _edibleFractionMarine           ;
    static double _AttackRateExponentMarine       ;
    static double _HandlingTimeExponentMarine     ;
    static double _HandlingTimeScalarMarine       ;
    static double _edibleFractionTerrestrial      ;
    static double _AttackRateExponentTerrestrial  ;
    static double _HandlingTimeExponentTerrestrial;
    static double _HandlingTimeScalarTerrestrial  ;
    static double _HerbivoryRateMassExponent      ;
    static double _HerbivoryRateConstant          ;
    static double _ReferenceMass                  ;
 
    //Carnivores

    static double _HandlingTimeScalar_C               ;
    static double _HandlingTimeExponent_C             ;
    static double _SearchRateConstant                 ;
    static double _FeedingPreferenceStandardDeviation ;
    static double _NumberOfBins                       ;

    //Dispersal constants
    static double _DispersalSpeedBodyMassScalar       ;
    static double _DispersalSpeedBodyMassExponent     ;
    //Advective dispersal
    static double _HorizontalDiffusivity                       ;
    static double _AdvectiveModelTimeStepLengthHours           ;
    static double _HorizontalDiffusivityKmSqPerADTimeStep      ;
    static double _AdvectionTimeStepsPerModelTimeStep          ;
    static double _VelocityUnitConversion                      ;
    //Responsive dispersal
    static double _DensityThresholdScaling                     ;
    static double _StarvationDispersalBodyMassThreshold        ;

    //ectotherms
    static double _TerrestrialWarmingToleranceIntercept ;
    static double _TerrestrialWarmingToleranceSlope     ;
    static double _TerrestrialTSMIntercept              ;
    static double _TerrestrialTSMSlope                  ;

    //metabolism
    //Ectotherms
    static double _MetabolismMassExponentEcto            ;
    static double _NormalizationConstantEcto             ;
    static double _ActivationEnergyEcto                  ; 
    static double _NormalizationConstantBMR              ;
    static double _BasalMetabolismMassExponent           ;
    static double _EnergyScalarEcto                      ;
    //Endotherrms
    static double _MetabolismMassExponentEndo           ;
    static double _NormalizationConstantEndo            ;
    static double _ActivationEnergyEndo                 ; 
    static double _EnergyScalarEndo                     ;

    static double _BoltzmannConstant                    ;
    static double _TemperatureUnitsConvert              ;

    //reproduction
    static double _MassRatioThreshold                   ;
    static double _MassEvolutionProbabilityThreshold    ;
    static double _MassEvolutionStandardDeviation       ;
    static double _SemelparityAdultMassAllocation       ;
    //mortality
    static double _MortalityRateBackground              ;
    static double _MortalityRateMature                  ;
    static double _LogisticInflectionPoint              ;
    static double _LogisticScalingParameter             ;
    static double _MaximumStarvationRate                ;
    
    //Pi!
    static double _Pi ;
    //area units conversion
    static double _CellAreaToHectares;
    
    //temporary store for within timestep changes
    static std::map < std::string, std::map<std::string,double> > _Accounting;



public:

    static unsigned _NextID;
    Cohort* _newH;
    Cohort(repast::AgentId id): MadAgent(id), _Merged(false){_NextID++;_newH=NULL;_sequencer=0;}
    //for copy across threads (needs increaseNextID=false) or restore from file (set increaseNextID to true)
	Cohort(repast::AgentId id, const AgentPackage& package,bool increaseNextID=false): MadAgent(id){PullThingsOutofPackage(package);_newH=NULL;if (increaseNextID)_NextID++;_sequencer=0;}
    void set(int currentRank, const AgentPackage& package){_id.currentRank(currentRank);PullThingsOutofPackage(package);}
	void setup(unsigned,unsigned,Environment*,randomizer*);
    void setPropertiesFromCohortDefinitions(unsigned);
	virtual ~Cohort() {}

	void step(Environment* ,vector<Cohort*>&,vector<Stock*>&,const unsigned,MadModel*);

    void metabolize(Environment*);
    void assignTimeActive(Environment*);
    void reproduce(Environment*);
    void eat(Environment*,vector<Cohort*>&,vector<Stock*>&,MadModel*);
    void moveIt(Environment*,MadModel*);
    void mort();
    void markForDeath();
    void applyEcology(Environment*);
    void updatePools(Environment*);
    void setupOffspring( Cohort* , double , double , double , double , unsigned  );
    void TryToDisperse(double,Environment*,MadModel* );
    void TryToDisperse(double,double,Environment*,MadModel*);
    vector<double> dProb(double,double,Environment*);
    vector<double> dDirect(double,double,Environment*);
    double distance(MadAgent*, MadAgent*,MadModel *);
    void PushThingsIntoPackage( AgentPackage& );
    void PullThingsOutofPackage( const AgentPackage& );
    void ResetAccounts();
};
#endif /* COHORT_H_ */
