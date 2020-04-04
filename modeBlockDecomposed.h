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
#include "EnvironmentCell.h"
#include "agent.h"


class MadModel;
class AgentPackage;
//This is an alias for the background discrete grid
typedef repast::SharedDiscreteSpace<MadAgent, repast::WrapAroundBorders, repast::SimpleAdder<MadAgent> > wrappedSpaceType;


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
    
    unsigned _startingStep;
	int _stopAt;
    
    bool _verbose;
    unsigned _randomSeed;
	repast::Properties* _props;
	repast::SharedContext<MadAgent> _context;
    std::vector<repast::DataSet*> dataSets;
	
	MadAgentPackageProvider* provider;
	MadAgentPackageReceiver* receiver;

    wrappedSpaceType* discreteSpace;
    int _totalCohorts,_totalStocks;
    double _totalCohortAbundance,_totalCohortBiomass,_totalStockBiomass,_totalOrganciPool,_totalRespiratoryCO2Pool;
    int    _totalMerged,_totalReproductions,_totalDeaths,_totalMoved;
    
    vector<int> _FinalCohortBreakdown;
      
    map< string,vector<double> > outputMaps;
    map<string,string> outputUnits;
    vector<string> outputNames;

    unsigned _restartInterval;
    unsigned _restartStep;
    std::string _restartDirectory;
    std::string _archiveFormat;
    
    std::string _filePrefix, _filePostfix;
    void dataSetClose();
    void addDataSet(repast::DataSet*) ;
    void setupNcOutput();
    void asciiOutput( unsigned step );
    void netcdfOutput( unsigned step );
    void setNcGridFile(std::string,std::string );
    void writeNcGridFile(unsigned,vector<double>&,std::string);
    std::vector<AgentPackage>_packages;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & _packages;
    }

public:
    bool _eating       ;
    bool _metabolism   ;
    bool _reproduction ;
    bool _death        ;
    bool _dispersal    ;
    bool _mergers      ;
    bool _output       ;
    
    int _minX,_minY,_maxX,_maxY,_dimX,_dimY;
    int _xlo,_xhi,_ylo,_yhi;
    int _noLongitudeWrap; //1 if domain does *not* span the global longitude range
    string _dispersalSelection;
    Environment _Env;
    
	MadModel(repast::Properties& ,  boost::mpi::communicator* comm);
	~MadModel();
	void init(unsigned);
	void moveAgents();
	void step();
    void sync();
	void initSchedule(unsigned, repast::ScheduleRunner& runner);
	void recordResults();
    wrappedSpaceType* space(){return discreteSpace;}

    static int _stockType, _cohortType, _humanType;
    //outputs
    void read_restart(unsigned);
    void write_restart();
    void setupOutputs();
    void tests();
    void setupCohortTestValues(Cohort*);
    void checkCohortTestValues(Cohort*);
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
