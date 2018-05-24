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
 * MadObserver.cpp
 *
 *  Created on: Sep 1, 2010
 *      Author: nick
 */
#include <sstream>

#include "relogo/RandomMove.h"
#include "relogo/grid_types.h"
#include "relogo/Patch.h"
#include "repast_hpc/AgentRequest.h"
#include "repast_hpc/SVDataSet.h"
#include "repast_hpc/SVDataSetBuilder.h"

#ifndef _WIN32
#include "repast_hpc/NCDataSetBuilder.h"
#endif

#include "MadObserver.h"
#include "CohortMerger.h"
#include "Environment.h"
#include "Groups.h"
#include "Stock.h"
#include "Cohort.h"
#include "FileReader.h"
#include "CohortSum.h"
#include "TimeStep.h"
#include "Constants.h"

using namespace std;
using namespace repast;
using namespace relogo;

//----------------------------------------------------------------------------------------------------------
void MadObserver::go() {
  if (_rank == 0) {
    //Log4CL::instance()->get_logger("root").log(INFO, "TICK BEGINS: " + boost::lexical_cast<string>(RepastProcess::instance()->getScheduleRunner().currentTick()));
  }
  //
  
   synchronize<AgentPackage>(*this, *this, *this, RepastProcess::USE_LAST_OR_USE_CURRENT);

  //sets the value used to get environmental fields (assumed monthly data) - NB data month currently starts at 0, but tick at 1 so subtract 1
  TimeStep::Get( )->SetMonthly( RepastProcess :: instance ()->getScheduleRunner ().currentTick () -1 );

  AgentSet<Environment> Env=patches<Environment>(); 

  //AgentSet<Stock> stocks;
  //get(stocks);
  
  _totalStockBiomass=0;
  //for (auto s:stocks)_totalStockBiomass+=s->_TotalBiomass/1000; //g to kg

  //stocks.ask(&Stock::step);

  AgentSet<Cohort> cohorts;

 // get(cohorts);
//cout<<"number "<<cohorts.size()<<endl;
//    std::map<int,int>counts;
        
//  _totalCohorts=cohorts.size();
//  _totalCohortAbundance=0;
//  _totalCohortBiomass=0;
//  for (auto c:cohorts){
//   _totalCohortAbundance += c->_CohortAbundance;
//   _totalCohortBiomass += ( c->_IndividualBodyMass + c->_IndividualReproductivePotentialMass ) * c->_CohortAbundance / 1000.;
//   counts[c->_FunctionalGroupIndex]++;
 // }
  //for (auto n:counts)cout<<n.first<<" "<<n.second<<endl;

 // cohorts.clear();
  for (auto e: Env) if (e->getId().currentRank()==repast::RepastProcess::instance()->rank()) e->step();
    for (auto e: Env) if (e->getId().currentRank()==repast::RepastProcess::instance()->rank()) e->flap();

  //NB - the folliwing would not be randomised: this is done in the above, which also localises to cells
  //cohorts.ask(&Cohort::step);
  //re-get the cohorts so as to include new births -otherwise merging acts on subsets in Environment that include new births, but not visible here

  get(cohorts);

  //cohorts.ask(&Cohort::markForDeath);
  //Merge cohorts, if necessary

  int totalMerged=0;
  //for (auto e: Env) e->merge(totalMerged);

 // cout<<"merged "<<totalMerged<<endl;
  //for (auto c:cohorts){if (c->getId().currentRank()==repast::RepastProcess::instance()->rank()) c->moveIt();}
  //cohorts.ask(&Cohort::moveIt);
  //cohorts.ask(&Cohort::expire);

  cohorts.clear();

  synchronize<AgentPackage>(*this, *this, *this, RepastProcess::USE_LAST_OR_USE_CURRENT);



  if (_rank == 0) {
    //Log4CL::instance()->get_logger("root").log(INFO, "TICK ENDS: " + boost::lexical_cast<string>(RepastProcess::instance()->getScheduleRunner().currentTick()));
  }

}

//----------------------------------------------------------------------------------------------------------
void MadObserver::setup(Properties& props) {
  //(SimRunner calls _setup, which calls this after all initialization is done)

  repast::Timer initTimer;
  initTimer.start();

  FileReader F;
  F.ReadFiles();
 Env.resize(x*y);
 forxy{Environment* E=new Environment(x,y);Env[x+xlen*y]=E;}
  for (auto E : Env){E->setup();}

  StockDefinitions::Initialise(Constants::cStockDefinitionsFileName);
  CohortDefinitions::Initialise(Constants::cCohortDefinitionsFileName);
  
  unsigned numCohortGroups=CohortDefinitions::Get()->size();

  unsigned cohortCount = strToInt(props.getProperty("cohort.count"));

  unsigned numStockGroups = StockDefinitions::Get()->size();;
  
  unsigned totalCohorts=0,totalStocks=0;
  for (auto E : Env){
    for (unsigned i=0;i<numCohortGroups;i++){

          if (E->_Realm==CohortDefinitions::Get()->Trait(i,"realm"))totalCohorts+=cohortCount;//cohortCount per functional group per cell
    }
    for (unsigned i=0;i<numStockGroups;i++){
          if (E->_Realm==StockDefinitions::Get()->Trait(i,"realm"))totalStocks+=1;//one stock per functional group
    }
  }
  //cout<<Env.size()<<" "<<totalCohorts<<endl;
  //DEBUG
  cohortType = create<Cohort> (totalCohorts);
  stockType  = create<Stock> (totalStocks);

  AgentSet<Cohort> cohorts;
  get(cohorts);
  AgentSet<Stock> stocks;
  get(stocks);

  unsigned cNum=0,sNum=0;
  for (auto E : Env){
      //DEBUG
  //if (E->getId().id()==0 && E->getId().startingRank()==0 && _rank==0) {
    unsigned totalCohortsThisCell=0;
    for (unsigned i=0;i<numCohortGroups;i++) if (E->_Realm==CohortDefinitions::Get()->Trait(i,"realm"))totalCohortsThisCell+=cohortCount;
    
    for (unsigned i=0;i<numCohortGroups;i++){
      if (E->_Realm==CohortDefinitions::Get()->Trait(i,"realm")){
        for (unsigned j=0;j<cohortCount;j++){
            cohorts[cNum]->moveTo(E);
            cohorts[cNum]->setup(i,totalCohortsThisCell);cNum++;
        }
      }
    }

    unsigned totalStocksThisCell=0;
    for (unsigned i=0;i<numStockGroups;i++){
     // if (E->_Realm==StockDefinitions::Get()->Trait(i,"realm")){stocks[sNum]->moveTo(E);stocks[sNum]->setup(i);totalStocksThisCell++;sNum++;}//one stock per functional group
    }
  //}
  }
  
//The things added to the datasetbuilder will be accumulated over cores each timestep and output to data.csv
//	SVDataSetBuilder svbuilder("./output/data.csv", ",", repast::RepastProcess::instance()->getScheduleRunner().schedule());
//	CohortSum* cSum = new CohortSum(this);
//	svbuilder.addDataSource(repast::createSVDataSource("Total Cohorts", cSum, std::plus<int>()));
//    StockBiomassSum* sSum = new StockBiomassSum(this);
//    CohortAbundanceSum* caSum = new CohortAbundanceSum(this);
//  	svbuilder.addDataSource(repast::createSVDataSource("Total Cohort Abundance", caSum, std::plus<double>()));
//  	svbuilder.addDataSource(repast::createSVDataSource("Total Stock Biomass", sSum, std::plus<double>()));
//    CohortBiomassSum* cbSum = new CohortBiomassSum(this);
//  	svbuilder.addDataSource(repast::createSVDataSource("Total Cohort Biomass", cbSum, std::plus<double>()));

//	addDataSet(svbuilder.createDataSet());
#ifndef _WIN32
	// no netcdf under windows?
//	NCDataSetBuilder builder("./output/data.ncf", RepastProcess::instance()->getScheduleRunner().schedule());
//	InfectionSum* infectionSum = new InfectionSum(this);
//	builder.addDataSource(repast::createNCDataSource("number_infected", infectionSum, std::plus<int>()));
//	addDataSet(builder.createDataSet());
#endif

	long double t = initTimer.stop();
	std::stringstream ss;
	ss << t;
	props.putProperty("init.time", ss.str());
    synchronize<AgentPackage>(*this, *this, *this, RepastProcess::USE_LAST_OR_USE_CURRENT);

}


//----------------------------------------------------------------------------------------------------------
//These methods are used in communication across processes
//----------------------------------------------------------------------------------------------------------
RelogoAgent* MadObserver::createAgent(const AgentPackage& content) {
	if (content.type == stockType) {
		return new Stock(content.getId(), this);
	} else if (content.type == cohortType) {
		return new Cohort(content.getId(), this, content);
	} else {
		// it's a patch.
      return new Environment(content.getId(), this);
	}

}
//----------------------------------------------------------------------------------------------------------
void MadObserver::provideContent(const repast::AgentRequest& request, std::vector<AgentPackage>& out) {
	const vector<AgentId>& ids = request.requestedAgents();
	for (int i = 0, n = ids.size(); i < n; i++) {
		AgentId id = ids[i];
        AgentPackage content = { id.id(), id.startingRank(), id.agentType(), id.currentRank()};
		if (id.agentType() == cohortType) {
			Cohort* cohort = who<Cohort> (id);
            cohort->SqodgeThingsIntoPackage(content);
		}
		out.push_back(content);
	}

}

//----------------------------------------------------------------------------------------------------------
void MadObserver::updateAgent(AgentPackage content){
  repast::AgentId id(content.id, content.proc, content.type);
  if (id.agentType() == cohortType) {
    Cohort * cohort = who<Cohort> (id);
    cohort->SuckThingsOutofPackage(content);
  }
}

//this never seems to get used? might be needed for other types of code
//----------------------------------------------------------------------------------------------------------
void MadObserver::provideContent(RelogoAgent* agent, std::vector<AgentPackage>& out) {
	AgentId id = agent->getId();
	AgentPackage content = { id.id(), id.startingRank(), id.agentType(), id.currentRank(), 0, false };
	if (id.agentType() == cohortType) {
		Cohort* cohort = static_cast<Cohort*> (agent);
        cohort->SqodgeThingsIntoPackage(content);
	}
	out.push_back(content);
}
//need once grid.buffer is greater than 0
//----------------------------------------------------------------------------------------------------------
void MadObserver::createAgents(std::vector<AgentPackage>& contents, std::vector<RelogoAgent*>& out) {
	for (size_t i = 0, n = contents.size(); i < n; ++i) {
		AgentPackage content = contents[i];
		if (content.type == stockType) {
			out.push_back(new Stock(content.getId(), this));
		} else if (content.type == cohortType) {
			out.push_back(new Cohort(content.getId(), this, content));
		} else {
			// it's a patch.
			out.push_back(new Environment(content.getId(), this));
		}
	}
}
