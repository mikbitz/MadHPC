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
 * Cohort.h
 *
 *  Created on: Sep 1, 2010
 *      Author: nick
 */

#ifndef COHORT_H_
#define COHORT_H_

#include "relogo/Turtle.h"
#include "AgentPackage.h"

class Cohort : public repast::relogo::Turtle {

public:
    bool _alive;
    std::string _Realm;
    bool _Merged;
    unsigned _FunctionalGroupIndex;
	bool  _Heterotroph, _Autotroph,_Endotherm,_Ectotherm, _Iteroparous,_Semelparous;
    bool _Herbivore,_Carnivore,_Omnivore;
    bool _IsMature,_IsPlanktonic,_IsFilterFeeder;

    double _ProportionTimeActive,_ProportionSuitableTimeActive;
    double _MinimumMass,_MaximumMass;
    double _IndividualBodyMass,_CohortAbundance,_AdultMass,_JuvenileMass,_IndividualReproductivePotentialMass,_MaximumAchievedBodyMass;
    double _LogOptimalPreyBodySizeRatio;
    double _PotentialAbundanceEaten;
    unsigned _MaturityTimeStep,_BirthTimeStep,_CurrentTimeStep;
    //herbivores
    double _AssimilationEfficiency_H;
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
    double _AssimilationEfficiency_C;
    const double _HandlingTimeScalar_C = 0.5;
    const double _HandlingTimeExponent_C = 0.7;
    const double _SearchRateConstant = 1E-6;
    const double _FeedingPreferenceStandardDeviation = 0.7;
    const double _NumberOfBins = 12;

    std::map < std::string, std::map<std::string,double> > _MassAccounting;

public:
//	Cohort(repast::AgentId id, repast::relogo::Observer* obs): repast::relogo::Turtle(id, obs), _infected(false), _infectionTime(0) {}
//	Cohort(repast::AgentId id, repast::relogo::Observer* obs, const AgentPackage& package): repast::relogo::Turtle(id, obs), _infected(package.infected),
//			_infectionTime(package.infectionTime) {}
    Cohort(repast::AgentId id, repast::relogo::Observer* obs): repast::relogo::Turtle(id, obs),_alive(true),_Merged(false){}
	Cohort(repast::AgentId id, repast::relogo::Observer* obs, const AgentPackage& package): repast::relogo::Turtle(id, obs){}

	void setup(unsigned,unsigned);

	virtual ~Cohort() {}

	void step();

    void metabolize();
    void assignTimeActive();
    void reproduce();
    void eat();
    void moveIt();
    void mort();
    void expire();
    void applyEcology();
    void setupOffspring( Cohort* , double , double , double , double , unsigned  );
    void TryToDisperse(double);
    void TryToDisperse(double,double);


};

#endif /* COHORT_H_ */
