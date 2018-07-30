/* Model.h */

#ifndef MODEL
#define MODEL

#include <boost/mpi.hpp>
#include "repast_hpc/Schedule.h"
#include "repast_hpc/Properties.h"
#include "repast_hpc/AgentRequest.h"
#include "repast_hpc/TDataSource.h"
#include "repast_hpc/SVDataSet.h"
#include "repast_hpc/SVDataSetBuilder.h"
#include "repast_hpc/SharedContext.h"
#include "repast_hpc/SharedDiscreteSpace.h"
#include "repast_hpc/GridComponents.h"
#include "Environment.h"
#include "agent.h"


class MadModel;
class AgentPackage;
//This is an alias for the background discrete grid
typedef repast::SharedDiscreteSpace<MadAgent, repast::WrapAroundBorders, repast::SimpleAdder<MadAgent> > modelSpaceType;


//------------------------------------------------------------------------------------------
//Structs for moving agents across threads
//------------------------------------------------------------------------------------------
class MadAgentPackageProvider {
	
private:
    repast::SharedContext<MadAgent>* agents;
	
public:
	
   MadAgentPackageProvider(repast::SharedContext<MadAgent>* agentPtr);

    void providePackage(MadAgent * agent, std::vector<AgentPackage>& out);

    void provideContent(repast::AgentRequest req, std::vector<AgentPackage>& out);
	
};
//------------------------------------------------------------------------------------------
class MadAgentPackageReceiver {
	
private:
    repast::SharedContext<MadAgent>* agents;
	
public:
	
    MadAgentPackageReceiver(repast::SharedContext<MadAgent>* agentPtr);
	
    MadAgent * createAgent(AgentPackage package);
	
    void updateAgent(AgentPackage package);
	
};
//------------------------------------------------------------------------------------------
//Actual model class
//------------------------------------------------------------------------------------------
class MadModel{
    
	int _stopAt;

	repast::Properties* _props;
	repast::SharedContext<MadAgent> _context;
    std::vector<repast::DataSet*> dataSets;
	
	MadAgentPackageProvider* provider;
	MadAgentPackageReceiver* receiver;

    modelSpaceType* discreteSpace;
    int _totalCohorts,_totalStocks;
    double _totalCohortAbundance,_totalCohortBiomass,_totalStockBiomass,_totalOrganciPool,_totalRespiratoryCO2Pool;
    int    _totalMerged,_totalReproductions,_totalDeaths,_totalMoved;
    void dataSetClose();
    void addDataSet(repast::DataSet*) ;

public:
    int _minX,_minY,_maxX,_maxY;
    int _xlo,_xhi,_ylo,_yhi;
    vector<Environment*> _Env;
	MadModel(std::string propsFile, int argc, char** argv, boost::mpi::communicator* comm);
	~MadModel();
	void init();
	void moveAgents();
	void step();
    void sync();
	void initSchedule(repast::ScheduleRunner& runner);
	void recordResults();
    modelSpaceType* space(){return discreteSpace;}
    static int _stockType, _cohortType;
    //outputs
    int CohortCount() const {
		return _totalCohorts;
	}
    int StockCount() const {
		return _totalStocks;
	}
	double CohortAbundance() const {
		return _totalCohortAbundance;
	}
	double CohortBiomass() const {
		return _totalCohortBiomass;
	}
	double StockBiomass() const {
		return _totalStockBiomass;
	}
	double OrganicPool() const {
		return _totalOrganciPool;
	}
	double RespiratoryCO2Pool() const {
		return _totalRespiratoryCO2Pool;
	}
	int Merged() const {
       return _totalMerged;
    }
    int Reproductions() const{
        return _totalReproductions;
    }
    int Deaths(){
        return _totalDeaths;
    }
    int Moved(){
        return _totalMoved;
    }

};
#endif
