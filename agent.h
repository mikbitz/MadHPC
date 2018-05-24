/* Agent.h */
#ifndef AGENT
#define AGENT

#include "repast_hpc/AgentId.h"

class MadAgent{
	
protected:
    repast::AgentId   _id;

public:
    MadAgent(){}
    MadAgent(repast::AgentId id): _id(id){}
	
    void set(int currentRank){    _id.currentRank(currentRank);}

    virtual ~MadAgent(){};
	
    /* Required Getters */
    virtual repast::AgentId& getId(){                   return _id;    }
    virtual const repast::AgentId& getId() const {      return _id;    }
	

};

/* Serializable Agent Package */
/* struct MadAgentPackage {
	
public:
    int    id;
    int    rank;
    int    type;
    int    currentRank;

	
    MadAgentPackage(){}; // For serialization
    MadAgentPackage(int _id, int _rank, int _type, int _currentRank):id(_id), rank(_rank), type(_type), currentRank(_currentRank){ }


MadAgentPackage::
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version){
        ar & id;
        ar & rank;
        ar & type;
        ar & currentRank;

    }

};
*/
#endif
