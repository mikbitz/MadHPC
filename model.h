/* Model.h */

#ifndef MODEL
#define MODEL

#include <boost/mpi.hpp>
#include "repast_hpc/Schedule.h"
#include "repast_hpc/Properties.h"
#include "repast_hpc/SharedContext.h"
#include "repast_hpc/AgentRequest.h"
#include "repast_hpc/TDataSource.h"
#include "repast_hpc/SVDataSet.h"
#include "repast_hpc/SVDataSetBuilder.h"
#include "repast_hpc/SharedDiscreteSpace.h"
#include "repast_hpc/GridComponents.h"

class RepastHPCAgent;
class RepastHPCAgentPackage;
typedef repast::SharedDiscreteSpace<RepastHPCAgent, repast::WrapAroundBorders, repast::SimpleAdder<RepastHPCAgent> > modelSpaceType;
/* Data Collection */
class DataSource_AgentTotals : public repast::TDataSource<int>{
private:
	repast::SharedContext<RepastHPCAgent>* context;
    
public:
	DataSource_AgentTotals(repast::SharedContext<RepastHPCAgent>* c);
	int getData();
};


class DataSource_AgentCTotals : public repast::TDataSource<int>{
private:
	repast::SharedContext<RepastHPCAgent>* context;
	
public:
	DataSource_AgentCTotals(repast::SharedContext<RepastHPCAgent>* c);
	int getData();
};

/* Agent Package Provider */
class RepastHPCAgentPackageProvider {
	
private:
    repast::SharedContext<RepastHPCAgent>* agents;
	
public:
	
    RepastHPCAgentPackageProvider(repast::SharedContext<RepastHPCAgent>* agentPtr);
	
    void providePackage(RepastHPCAgent * agent, std::vector<RepastHPCAgentPackage>& out);
	
    void provideContent(repast::AgentRequest req, std::vector<RepastHPCAgentPackage>& out);
	
};

/* Agent Package Receiver */
class RepastHPCAgentPackageReceiver {
	
private:
    repast::SharedContext<RepastHPCAgent>* agents;
	
public:
	
    RepastHPCAgentPackageReceiver(repast::SharedContext<RepastHPCAgent>* agentPtr);
	
    RepastHPCAgent * createAgent(RepastHPCAgentPackage package);
	
    void updateAgent(RepastHPCAgentPackage package);
	
};


class RepastHPCModel{
	int stopAt;
    int gridBuffer;
	int countOfAgents;
    double minX,minY,maxX,maxY;
	repast::Properties* props;
	repast::SharedContext<RepastHPCAgent> context;
	
	RepastHPCAgentPackageProvider* provider;
	RepastHPCAgentPackageReceiver* receiver;

	repast::SVDataSet* agentValues;
    modelSpaceType* discreteSpace;
	
public:
	RepastHPCModel(std::string propsFile, int argc, char** argv, boost::mpi::communicator* comm);
	~RepastHPCModel();
	void init();
	void moveAgents();
	void step();
    void sync();
	void initSchedule(repast::ScheduleRunner& runner);
	void recordResults();
    modelSpaceType* space(){return discreteSpace;}
};
#endif
