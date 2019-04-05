/*
 *
 *
 * Stock.h
 *
 *  Created on: May 8, 2018
 *      Author: Mike Bithell
 */

#ifndef STOCK_H_
#define STOCK_H_

#include "repast_hpc/AgentId.h"
#include "agent.h"
#include "AgentPackage.h"

class Environment;

class Stock: public MadAgent  {
public:
    Stock(repast::AgentId id,Environment*) : MadAgent(id){ }
    Stock(repast::AgentId id,AgentPackage& package): MadAgent(id){PullThingsOutofPackage( package);}
	virtual ~Stock() {}
    void setup(unsigned,Environment* );
	void step(double&,Environment*,const unsigned);
    double _TotalBiomass;
    bool _Marine,_Deciduous;
    unsigned _FunctionalGroupIndex;
    double _IndividualBodyMass;

    void PushThingsIntoPackage( AgentPackage& );
    void PullThingsOutofPackage( const AgentPackage& );              
};

#endif /* STOCK_H_ */
