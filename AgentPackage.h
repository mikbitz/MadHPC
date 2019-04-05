/*
 *  AgentPackage.h
 *  Created on: May 7, 2018
 *      Author: Mike Bithell
 * 
 * based in part on Repast for High Performance Computing (Repast HPC)
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
 * AgentContent.h
 *
 *  Created on: Sep 2, 2010
 *      Author: nick
 *      Customised for Cohorts and Stocks MB 2019
 */

#ifndef AGENTPACKAGE_H
#define AGENTPACKAGE_H

#include <vector>
#include "repast_hpc/AgentId.h"
//content is a bit overspecified at the moment as it tries to cover cohorts, stocks and huamns
//however, not obvious to me how to change this! (tried with polymorphic pointers, but massive memeory leaks or seg. faults)
struct content {
    content(){}

    unsigned _FunctionalGroupIndex;   
    //needed for stocks
    double _TotalBiomass;
    bool _Marine;
    bool _Deciduous;
    //---
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
    
    bool _moved;
    std::vector<double> _location,_destination;
    
	template<class Archive>
	void serialize(Archive& ar, const unsigned int version) {
        ar & _FunctionalGroupIndex;   
        ar & _TotalBiomass;
        ar & _Marine;
        ar & _Deciduous;
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
        
        ar & _moved;
        ar & _location;
        ar & _destination;
    }
};

struct AgentPackage {
    AgentPackage(){}
    AgentPackage(repast::AgentId ID):_id(ID){}
	repast::AgentId getId() const {
		return _id;
	}
	void setId(repast::AgentId ID)  {
      _id=ID;
	}
    repast::AgentId _id;
    
    content _contents;
    
	template<class Archive>
	void serialize(Archive& ar, const unsigned int version) {
        ar & _id;
        ar & _contents;

	}

};


#endif /* AGENTCONTENT_H_ */
