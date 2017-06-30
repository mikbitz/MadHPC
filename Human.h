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

#include "relogo/Turtle.h"
#include "AgentPackage.h"

class Human : public repast::relogo::Turtle {

private:
	bool _infected, _Heterotroph, _Endotherm,_Ectotherm,_Marine,_Terrestrial, _Iteroparous;
    bool _IsMature,_IsPlanktonic,_IsFilterFeeder;
    bool _alive;
    bool _Merged;
    bool _Herbivore,_Carnivore,_Omnivore;
    double _ProportionTimeActive,_ProportionSuitableTimeActive;
    double _MinimumMass,_MaximumMass;
    double _IndividualBodyMass,_CohortAbundance,_AdultMass,_JuvenileMass,_IndividualReproductivePotentialMass,_MaximumAchievedBodyMass;
    double _LogOptimalPreyBodySizeRatio;
    double _PotentialAbundanceEaten;
    unsigned _MaturityTimeStep,_BirthTimeStep,_CurrentTimeStep;
    unsigned _FunctionalGroupIndex;
//	int _infectionTime;
    double _edibleFractionMarine;
    double _AttackRateExponentMarine;
    double _HandlingTimeExponentMarine;
    double _HandlingTimeScalarMarine;
    double _edibleFractionTerrestrial;
    double _AttackRateExponentTerrestrial;
    double  _HandlingTimeExponentTerrestrial;
    double  _HandlingTimeScalarTerrestrial;
    double _HerbivoryRateMassExponent;
    double _HerbivoryRateConstant;
    double _ReferenceMass;
    double _AssimilationEfficiency;
    std::map < std::string, std::map<std::string,double> > _MassAccounting;

public:
//	Human(repast::AgentId id, repast::relogo::Observer* obs): repast::relogo::Turtle(id, obs), _infected(false), _infectionTime(0) {}
//	Human(repast::AgentId id, repast::relogo::Observer* obs, const AgentPackage& package): repast::relogo::Turtle(id, obs), _infected(package.infected),
//			_infectionTime(package.infectionTime) {}
	Human(repast::AgentId id, repast::relogo::Observer* obs): repast::relogo::Turtle(id, obs) {}
	Human(repast::AgentId id, repast::relogo::Observer* obs, const AgentPackage& package): repast::relogo::Turtle(id, obs){}
	virtual ~Human() {}

	void step();

//	void infect();
    void metabolize();
    void assignTimeActive();
    void reproduce();
    void eat();
    void moveIt();
    void mort();
    void expire();
    void applyEcology();
    void setOffspring( Human* , double , double , double , double , unsigned  );
    void TryToDisperse(double);
    void TryToDisperse(double,double);

//	bool infected() const {
//		return _infected;
//	}

//	int infectionTime() const {
//		return _infectionTime;
//	}

 // void set(bool infected, int infectionTime){
 //   _infected = infected;
 //   _infectionTime = infectionTime;
 // }
};

#endif /* HUMAN_H_ */
