/*
*Repast for High Performance Computing (Repast HPC)
*
*   Copyright (c) 2010 Argonne National Laboratory
*   All rights reserved.
*  
*   Redistribution and use in source and binary forms, with 
*   or without modification, are permitted provided that the following 
*   conditions are met:
*  
*  	 Redistributions of source code must retain the above copyright notice,
*  	 this list of conditions and the following disclaimer.
*  
*  	 Redistributions in binary form must reproduce the above copyright notice,
*  	 this list of conditions and the following disclaimer in the documentation
*  	 and/or other materials provided with the distribution.
*  
*  	 Neither the name of the Argonne National Laboratory nor the names of its
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
*/

/*
 * AgentContent.h
 *
 *  Created on: Sep 2, 2010
 *      Author: nick
 */

#ifndef AGENTPACKAGE_H
#define AGENTPACKAGE_H

#include "repast_hpc/AgentId.h"

struct AgentPackage {
    AgentPackage(){}
    AgentPackage(int id, int rank, int type, int currentRank):_id(id), _rank(rank), _type(type), _currentRank(currentRank){ }
	template<class Archive>
	void serialize(Archive& ar, const unsigned int version) {
		ar & _id;
		ar & _rank;
		ar & _type;
		ar & _currentRank;

		ar & _FunctionalGroupIndex;   

		ar & _JuvenileMass;     
		ar & _AdultMass;               
		ar & _IndividualBodyMass;
		ar & _MaximumAchievedBodyMass;
		ar & _IndividualReproductivePotentialMass;
		ar & _MinimumMass;
		ar & _MaximumMass;
    
		ar & _CohortAbundance;
		ar & _BirthTimeStep;            
		ar & _MaturityTimeStep;            
		ar & _LogOptimalPreyBodySizeRatio; 
     
		ar & _Merged;
		ar & _alive;                     

  
		ar & _Heterotroph;   
		ar & _Autotroph; 
		ar & _Endotherm; 
		ar & _Ectotherm;
    
		ar & _Realm;      

		ar & _Iteroparous;
		ar & _Semelparous;
		ar & _Herbivore;
		ar & _Carnivore;
		ar & _Omnivore;
		ar & _IsPlanktonic;
		ar & _IsFilterFeeder;

    
		ar & _ProportionSuitableTimeActive;
		ar & _IsMature;
    
		ar & _AssimilationEfficiency_H;
		ar & _AssimilationEfficiency_C;
	}

	int _id, _rank, _type, _currentRank;


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
    bool _alive;                     

  
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

	repast::AgentId getId() const {
		return repast::AgentId(_id, _rank, _type, _currentRank);
	}
};



#endif /* AGENTCONTENT_H_ */
