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

class EnvironmentCell;

class Stock: public MadAgent  {
public:
    Stock(repast::AgentId id,EnvironmentCell*) : MadAgent(id){ }
    Stock(repast::AgentId id,AgentPackage& package): MadAgent(id){PullThingsOutofPackage( package);}
	virtual ~Stock() {}
    void setup(unsigned,EnvironmentCell* );
	void step(double&,EnvironmentCell*,const unsigned);
    double _TotalBiomass;
    bool _Marine,_Deciduous;
    unsigned _FunctionalGroupIndex;
    double _IndividualBodyMass;
    void setPropertiesFromFunctionalGroupIndex( unsigned& );
    void PushThingsIntoPackage( AgentPackage& );
    void PullThingsOutofPackage( const AgentPackage& );              
};

#endif /* STOCK_H_ */
