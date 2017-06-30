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
#include "InfectionSum.h"
#include "TimeStep.h"
#include "Constants.h"

using namespace std;
using namespace repast;
using namespace relogo;

//----------------------------------------------------------------------------------------------------------
void MadObserver::go() {
  if (_rank == 0) {
    Log4CL::instance()->get_logger("root").log(INFO, "TICK BEGINS: " + boost::lexical_cast<string>(RepastProcess::instance()->getScheduleRunner().currentTick()));
  }
  synchronize<AgentPackage>(*this, *this, *this, RepastProcess::USE_LAST_OR_USE_CURRENT);
  
  //sets the value used to get environemtal fields (assumed monthly data)
  TimeStep::Get( )->SetMonthly( RepastProcess :: instance ()->getScheduleRunner ().currentTick () );
  
  AgentSet<Stock> stocks;
  get(stocks);
  stocks.ask(&Stock::step);

  AgentSet<Cohort> cohorts;
  get(cohorts);

  cohorts.ask(&Cohort::step);
  //cohorts.ask(&Cohort::expire);
  AgentSet<Environment> Env=patches<Environment>();

  // Merge cohorts, if necessary
  //Env.ask(&Environment::merge);
  //cohorts.ask(&Cohort::expire);

  //cohorts.ask(&Cohort::moveIt);
  
  if (_rank == 0) {
    Log4CL::instance()->get_logger("root").log(INFO, "TICK ENDS: " + boost::lexical_cast<string>(RepastProcess::instance()->getScheduleRunner().currentTick()));
  }
}

//----------------------------------------------------------------------------------------------------------
void MadObserver::setup(Properties& props) {

  repast::Timer initTimer;
  initTimer.start();

  FileReader F;
  F.ReadFiles();

  AgentSet<Environment> Env=patches<Environment>();

  for (auto E : Env){E->setup();}

  StockDefinitions::Initialise(Constants::cStockDefinitionsFileName);
  CohortDefinitions::Initialise(Constants::cCohortDefinitionsFileName);
  
  unsigned numCohortGroups=CohortDefinitions::Get()->size();
  
  unsigned cohortCount = strToInt(props.getProperty("cohort.count"));

  unsigned numStockGroups = StockDefinitions::Get()->size();;
  
  unsigned totalCohorts=0,totalStocks=0;
  for (auto E : Env){
    for (unsigned i=0;i<numCohortGroups;i++){
          if (E->_Realm==CohortDefinitions::Get()->Trait(i,"Realm"))totalCohorts+=cohortCount;//one per functional group per cell
    }
    for (unsigned i=0;i<numStockGroups;i++){
          if (E->_Realm==StockDefinitions::Get()->Trait(i,"Realm"))totalStocks+=1;//one stock per functional group
    }
  }

  cohortType = create<Cohort> (totalCohorts);
  stockType  = create<Stock> (totalStocks);

  AgentSet<Cohort> cohorts;
  get(cohorts);
  AgentSet<Stock> stocks;
  get(stocks);
  unsigned cOffset=0,sOffset=0;
  
  for (auto E : Env){
  
    unsigned totalCohortsThisCell=0;
    for (unsigned i=0;i<numCohortGroups;i++) if (E->_Realm==CohortDefinitions::Get()->Trait(i,"Realm"))totalCohortsThisCell+=cohortCount;
    
    for (unsigned i=0;i<numCohortGroups;i++){
      if (E->_Realm==CohortDefinitions::Get()->Trait(i,"Realm")){
        for (unsigned j=0;j<cohortCount;j++){
            cohorts[i*cohortCount+j+cOffset]->moveTo(E);
            cohorts[i*cohortCount+j+cOffset]->setup(i,totalCohortsThisCell);
        }
      }
    }
    cOffset+=totalCohortsThisCell;

    unsigned totalStocksThisCell=0;
    for (unsigned i=0;i<numStockGroups;i++){
      if (E->_Realm==StockDefinitions::Get()->Trait(i,"Realm")){stocks[i+sOffset]->moveTo(E);stocks[i+sOffset]->setup(i);totalStocksThisCell++;}//one stock per functional group
    }
    sOffset+=totalStocksThisCell;
  }
//	SVDataSetBuilder svbuilder("./output/data.csv", ",", repast::RepastProcess::instance()->getScheduleRunner().schedule());
//	InfectionSum* iSum = new InfectionSum(this);
//	svbuilder.addDataSource(repast::createSVDataSource("number_infected", iSum, std::plus<int>()));
//	addDataSet(svbuilder.createDataSet());


#ifndef _WIN32
	// no netcdf under windows
//	NCDataSetBuilder builder("./output/data.ncf", RepastProcess::instance()->getScheduleRunner().schedule());
//	InfectionSum* infectionSum = new InfectionSum(this);
//	builder.addDataSource(repast::createNCDataSource("number_infected", infectionSum, std::plus<int>()));
//	addDataSet(builder.createDataSet());
#endif

	long double t = initTimer.stop();
	std::stringstream ss;
	ss << t;
	props.putProperty("init.time", ss.str());
}


//----------------------------------------------------------------------------------------------------------

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
		AgentPackage content = { id.id(), id.startingRank(), id.agentType(), id.currentRank(), 0, false };
		if (id.agentType() == cohortType) {
			Cohort* cohort = who<Cohort> (id);
//			content.infected = cohort->infected();
//			content.infectionTime = cohort->infectionTime();
		}
		out.push_back(content);
	}
}
//----------------------------------------------------------------------------------------------------------
void MadObserver::provideContent(RelogoAgent* agent, std::vector<AgentPackage>& out) {
	AgentId id = agent->getId();
	AgentPackage content = { id.id(), id.startingRank(), id.agentType(), id.currentRank(), 0, false };
	if (id.agentType() == cohortType) {
		Cohort* cohort = static_cast<Cohort*> (agent);
//		content.infected = cohort->infected();
//		content.infectionTime = cohort->infectionTime();
	}
	out.push_back(content);
}
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
//----------------------------------------------------------------------------------------------------------
void MadObserver::updateAgent(AgentPackage package){
  repast::AgentId id(package.id, package.proc, package.type);
  if (id.agentType() == cohortType) {
    Cohort * cohort = who<Cohort> (id);
//    cohort->set(package.infected, package.infectionTime);
  }
}
