/* Agent.h */
#ifndef AGENT
#define AGENT
#include <vector>
#include "repast_hpc/AgentId.h"

class MadAgent{
	
protected:
    repast::AgentId   _id;

public:
    MadAgent(){_moved=false;_alive=true;_location={0,0};_destination=_location;}
    MadAgent(repast::AgentId id): _id(id){_moved=false;_alive=true;_location={0,0};_destination=_location;}
	bool _moved;
    bool _alive;

    virtual ~MadAgent(){};
	
    /* Required Getters */
    virtual repast::AgentId& getId(){                   return _id;    }
    virtual const repast::AgentId& getId() const {      return _id;    }
	
    //locations in *fractions of a grid cell* at a given long/lat.
    std::vector<double> _location,_destination;
    std::vector<double> getLocation(){return _location;}
    void setLocation(double x, double y){_location={x,y};}
    void setDestination(double x, double y){_destination={x,y};}
    void setLocation(std::vector<double>d){_location=d;}
    void setDestination(std::vector<double>d){_destination=d;}
};

#endif
