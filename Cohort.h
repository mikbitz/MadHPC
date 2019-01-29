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

    static std::map < std::string, std::map<std::string,double> > _Accounting;
    //std::map < std::string, std::map<std::string,double> > _Accounting;


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
