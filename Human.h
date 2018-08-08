/*
 *   Repast for High Performance Computing (Repast HPC)
 *
 *   Copyright (c) 2010 Argonne National Laboratory
 *   All rights reserved.
 *
 *   Redistribution and use in source and binary forms, with
 *   or without modification, are permitted provided that the following
 *   conditions are met:
 *
 *     Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *
 *     Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *
 *     Neither the name of the Argonne National Laboratory nor the names of its
 *     contributors may be used to endorse or promote products derived from
 *     this software without specific prior written permission.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *   PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE TRUSTEES OR
 *   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 *   EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 * Human.h
 *
 *  Created on: Sep 1, 2010
 *      Author: nick
 */

#ifndef HUMAN_H_
#define HUMAN_H_

/*
 *
 *
 * Human.h
 *
 *  Created on: Aug 8, 2018
 *      Author: Mike Bithell
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

    const double _edibleFractionMarine            =1.0;
    const double _AttackRateExponentMarine        =2.0;
    const double _HandlingTimeExponentMarine      =0.7;
    const double _HandlingTimeScalarMarine        =0.7;
    const double _edibleFractionTerrestrial       =0.1;
    const double _AttackRateExponentTerrestrial   =2.0;
    const double _HandlingTimeExponentTerrestrial =0.7;
    const double _HandlingTimeScalarTerrestrial   =0.7;
    const double _HerbivoryRateMassExponent       =1.0;
    const double _HerbivoryRateConstant           =1.0e-11;
    const double _ReferenceMass                   =1.0;
 
    //Carnivores

    const double _HandlingTimeScalar_C = 0.5;
    const double _HandlingTimeExponent_C = 0.7;
    const double _SearchRateConstant = 1E-6;
    const double _FeedingPreferenceStandardDeviation = 0.7;
    const double _NumberOfBins = 12;

    std::map < std::string, std::map<std::string,double> > _MassAccounting;


public:

    static unsigned _NextID;
    Human* _newH;
    Human(repast::AgentId id): MadAgent(id), _Merged(false){_NextID++;_newH=NULL;}
	Human(repast::AgentId id, const AgentPackage& package): MadAgent(id){PullThingsOutofPackage(package);_newH=NULL;}
    void set(int currentRank, const AgentPackage& package){_id.currentRank(currentRank);PullThingsOutofPackage(package);}
	void setup(unsigned,unsigned,Environment*,randomizer&);

	virtual ~Human() {}

	void step(Environment* ,vector<Cohort*>&,vector<Stock*>&,const unsigned);

    void metabolize(Environment*);
    void assignTimeActive(Environment*);
    void reproduce(Environment*);
    void eat(Environment*,vector<Cohort*>&,vector<Stock*>&);
    void moveIt(Environment*,MadModel*);
    void relocateBy(int,int, MadModel*);
    void mort();
    void markForDeath();
    void applyEcology(Environment*);
    void updatePools(Environment*);
    void setupOffspring( Human* , double , double , double , double , unsigned  );
    vector<int> TryToDisperse(double,Environment*,MadModel* );
    vector<int> TryToDisperse(double,double,Environment*,MadModel* );
    vector<int> _displacement;

void PushThingsIntoPackage( AgentPackage& );
void PullThingsOutofPackage( const AgentPackage& );
void ResetMassFluxes();
};



#endif /* HUMAN_H_ */
