/* Agent.h */
#ifndef AGENT
#define AGENT

#include "repast_hpc/AgentId.h"

class MadAgent{
	
protected:
    repast::AgentId   _id;

public:
    MadAgent(){_moved=false;_alive=true;}
    MadAgent(repast::AgentId id): _id(id){_moved=false;_alive=true;}
	bool _moved;
    bool _alive;
    void set(int currentRank){    _id.currentRank(currentRank);}//is this correct/needed?

    virtual ~MadAgent(){};
	
    /* Required Getters */
    virtual repast::AgentId& getId(){                   return _id;    }
    virtual const repast::AgentId& getId() const {      return _id;    }
	

};

#endif
