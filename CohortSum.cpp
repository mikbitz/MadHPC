/*
 * CohortSum.cpp
 * Created on: May 7 2018
 * Author: Mike Bithell
 *
 *   Based on Infectionsum.cpp
 * Repast for High Performance Computing (Repast HPC)
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
 *  Created on: Sep 2, 2010
 *      Author: nick
 */

#include "CohortSum.h"
#include "model.h"

CohortSum::CohortSum(MadModel* Mobs) : obs(Mobs){}

CohortSum::~CohortSum() {}

int CohortSum::getData() {
	return obs->CohortCount();
}
//----------------------------------------------------------
StockSum::StockSum(MadModel* Mobs) : obs(Mobs){}

StockSum::~StockSum() {}

int StockSum::getData() {
	return obs->StockCount();
}
//----------------------------------------------------------
CohortAbundanceSum::CohortAbundanceSum(MadModel* Mobs) : obs(Mobs){}

CohortAbundanceSum::~CohortAbundanceSum() {}

double CohortAbundanceSum::getData() {
	return obs->CohortAbundance();
}
//----------------------------------------------------------
CohortBiomassSum::CohortBiomassSum(MadModel* Mobs) : obs(Mobs){}

CohortBiomassSum::~CohortBiomassSum() {}

double CohortBiomassSum::getData() {
	return obs->CohortBiomass();
}
//----------------------------------------------------------
StockBiomassSum::StockBiomassSum(MadModel* Mobs) : obs(Mobs){}

StockBiomassSum::~StockBiomassSum() {}

double StockBiomassSum::getData() {
	return obs->StockBiomass();
}
//----------------------------------------------------------
CohortOrganicPool::CohortOrganicPool (MadModel* Mobs) : obs(Mobs){}
CohortOrganicPool::~CohortOrganicPool(){}

double CohortOrganicPool::getData(){
    return obs->OrganicPool();
}
//----------------------------------------------------------
CohortResp::CohortResp (MadModel* Mobs) : obs(Mobs){}
CohortResp::~CohortResp(){}

double CohortResp::getData(){
   return obs->RespiratoryCO2Pool() ;
    
}
//----------------------------------------------------------
DispersalSum::DispersalSum (MadModel* Mobs) : obs(Mobs){}
DispersalSum::~DispersalSum(){}

double DispersalSum::getData(){
   return obs->Moved() ;
    
}
//----------------------------------------------------------
ExtinctionSum::ExtinctionSum (MadModel* Mobs) : obs(Mobs){}
ExtinctionSum::~ExtinctionSum(){}

double ExtinctionSum::getData(){
   return obs->Deaths() ;
    
}
//----------------------------------------------------------
ProductionSum::ProductionSum (MadModel* Mobs) : obs(Mobs){}
ProductionSum::~ProductionSum(){}

double ProductionSum::getData(){
   return obs->Reproductions() ;
    
}
//----------------------------------------------------------
CombinationSum::CombinationSum (MadModel* Mobs) : obs(Mobs){}
CombinationSum::~CombinationSum(){}

double CombinationSum::getData(){
   return obs->Merged() ;
    
}
