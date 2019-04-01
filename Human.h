/*
 *
 *
 * Human.h
 *  Created on: Aug 8, 2018
 *  Author: Mike Bithell
 *  Partially derived from Original C# code by
 * Emergent Global Patterns of Ecosystem Structure and Function from a Mechanistic General Ecosystem Model , Michael B. J. Harfoot , Tim Newbold, Derek P. Tittensor,  Stephen Emmott, Jon Hutton, Vassily Lyutsarev, Matthew J. Smith, Jörn P. W. Scharlemann, Drew W. Purves PLOS Biology 12(4) 2014 e1001841 https://doi.org/10.1371/journal.pbio.1001841
 *
 */

#ifndef HUMAN_H_
#define HUMAN_H_

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

class Human: public MadAgent  {

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
   
    std::string _Realm;      

   
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
    //Responsive dispersal
    static double _DensityThresholdScaling                     ;
    static double _StarvationDispersalBodyMassThreshold        ;

    //metabolism
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
    Human* _newH;
    Human(repast::AgentId id): MadAgent(id), _Merged(false){_NextID++;_newH=NULL;}
    //for copy across threads (needs increaseNextID=false) or restore from file (set increaseNextID to true)
	Human(repast::AgentId id, const AgentPackage& package,bool increaseNextID=false): MadAgent(id){PullThingsOutofPackage(package);_newH=NULL;if (increaseNextID)_NextID++;}
    void set(int currentRank, const AgentPackage& package){_id.currentRank(currentRank);PullThingsOutofPackage(package);}
	void setup(unsigned,unsigned,Environment*,randomizer*);

	virtual ~Human() {}

	void step(Environment* ,vector<Human*>&,vector<Stock*>&,const unsigned);

    void metabolize(Environment*);
    void assignTimeActive(Environment*);
    void reproduce(Environment*);
    void eat(Environment*,vector<Human*>&,vector<Stock*>&);
    void moveIt(Environment*,MadModel*);
    void mort();
    void markForDeath();
    void applyEcology(Environment*);
    void updatePools(Environment*);
    void setupOffspring( Human* , double , double , double , double , unsigned  );
    void TryToDisperse(double,Environment*,MadModel* );
    void TryToDisperse(double,double,Environment*,MadModel* );

void PushThingsIntoPackage( AgentPackage& );
void PullThingsOutofPackage( const AgentPackage& );
void ResetAccounts();
};


#endif /* HUMAN_H_ */
