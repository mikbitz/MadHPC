/* Model.cpp */

#include <stdio.h>
#include <vector>
#include <sstream>
#include <boost/mpi.hpp>
#include "repast_hpc/AgentId.h"
#include "repast_hpc/RepastProcess.h"
#include "repast_hpc/Utilities.h"
#include "repast_hpc/Properties.h"
#include "repast_hpc/SharedBaseGrid.h"                      // VN2D broken without this!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "repast_hpc/VN2DGridQuery.h"
#include "repast_hpc/initialize_random.h"
#include "repast_hpc/Point.h"
#include "repast_hpc/AgentRequest.h"
#include "repast_hpc/SVDataSet.h"
#include "repast_hpc/SVDataSetBuilder.h"

#ifndef _WIN32
#include "repast_hpc/NCDataSetBuilder.h"
#endif

#include "model.h"
#include "agent.h"
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

//------------------------------------------------------------------------------------------------------------
//Constructor and destructor
//------------------------------------------------------------------------------------------------------------
MadModel::MadModel(std::string propsFile, int argc, char** argv, boost::mpi::communicator* comm): _context(comm){
    //-----------------
    //Pull in parameter from the model.props file
    _props = new repast::Properties(propsFile, argc, argv, comm);
    //Number of timesteps
	_stopAt = repast::strToInt(_props->getProperty("stop.at"));
    //extent of buffer zones in grid units - this many grid cells are sahred at the boundary between cores
    int gridBuffer = repast::strToInt(_props->getProperty("grid.buffer"));
    //Grid extent
    _minX=repast::strToInt(_props->getProperty("min.x"));
    _minY=repast::strToInt(_props->getProperty("min.y"));
    _maxX=repast::strToInt(_props->getProperty("max.x"));
    _maxY=repast::strToInt(_props->getProperty("max.y"));
    int dimX=repast::strToInt(_props->getProperty("proc.per.x"));
    int dimY=repast::strToInt(_props->getProperty("proc.per.y"));

    //-----------------
	//create the model grid
    repast::Point<double> origin(_minX,_minY);
    repast::Point<double> extent(_maxX-_minX+1, _maxY-_minY+1);
    
    repast::GridDimensions gd(origin, extent);
    
    std::vector<int> processDims;
    processDims.push_back(dimX);
    processDims.push_back(dimY);
  
    discreteSpace = new repast::SharedDiscreteSpace<MadAgent, repast::WrapAroundBorders, repast::SimpleAdder<MadAgent> >("AgentDiscreteSpace", gd, processDims, gridBuffer, comm);
	
   
    //The agent container is a context. Add the grid to it.
   	_context.addProjection(discreteSpace);
    
   //-----------------
    //Set up the cross-thread data transferers
//	provider = new MadAgentPackageProvider(&_context);
//	receiver = new MadAgentPackageReceiver(&_context);
    
}
//------------------------------------------------------------------------------------------------------------
MadModel::~MadModel(){
	delete _props;
//	delete provider;
//	delete receiver;

}
//------------------------------------------------------------------------------------------------------------
//Model Initialisation
//------------------------------------------------------------------------------------------------------------
void MadModel::initSchedule(repast::ScheduleRunner& runner){
	runner.scheduleEvent(1, 1, repast::Schedule::FunctorPtr(new repast::MethodFunctor<MadModel> (this, &MadModel::step)));
	//runner.scheduleEndEvent(repast::Schedule::FunctorPtr(new repast::MethodFunctor<MadModel> (this, &MadModel::recordResults)));
	runner.scheduleStop(_stopAt);
}
//------------------------------------------------------------------------------------------------------------
void MadModel::init(){
    //time the initialisation
    repast::Timer initTimer;
    initTimer.start();
    //rank (i.e. the number of this thread) will be needed to make agents unique *between* threads
    int rank = repast::RepastProcess::instance()->rank();
    
    //get the environmental data - this is stored in the background as a DataLayerSet
    FileReader F;
    F.ReadFiles();
    
    //now set up the environmental cells - note at present this uses the full grid, not just local to this thread
    //so that off-thread environment can be easily queried. Currently some duplication here, but it is not a huge amount of data.
    _Env.resize( (_maxX-_minX+1) * (_maxY-_minY+1) );
    for(int x=_minX;x<=_maxX;x++) {
        for(int y=_minY;y<=_maxY;y++) {
            Environment* E=new Environment(x,y);
            _Env[x-_minX+(_maxX-_minX+1)*(y-_minY)]=E;
        }
    }
 
    //get the definitions of stocks and cohorts
    StockDefinitions::Initialise(Constants::cStockDefinitionsFileName);
    CohortDefinitions::Initialise(Constants::cCohortDefinitionsFileName);
  
    unsigned numCohortGroups=CohortDefinitions::Get()->size();
    unsigned cohortCount = strToInt(_props->getProperty("cohort.count"));
    unsigned numStockGroups = StockDefinitions::Get()->size();
    
    //arbitrary numbers to distiguish the agents types
    _cohortType = 0;
    _stockType  = 1;
    
    //explicitly use the local bounds of the grid on this thread to create countOfAgents per cell.
    //Not doing this can lead to problems with agents in dsitant cells not within the local thread neighbourhood
    //see SharedBaseGrid.h moveTo method
    int xmin=        discreteSpace->dimensions().origin().getX() ;
    int xmax= xmin + discreteSpace->dimensions().extents().getX();
    int ymin=        discreteSpace->dimensions().origin().getY() ;
    int ymax= ymin + discreteSpace->dimensions().extents().getY();
     unsigned totalCohorts=0,totalStocks=0;

  unsigned cNum=0,sNum=0;

    unsigned totalStocksThisCell=0;

  
    int s=0;
    for (int x = xmin; x < xmax; x++){
        for (int y = ymin; y < ymax; y++){
             Environment* E=_Env[x-_minX+(_maxX-_minX+1)*(y-_minY)];
             unsigned totalCohortsThisCell=0;
             for (unsigned i=0;i<numCohortGroups;i++) if (E->_Realm==CohortDefinitions::Get()->Trait(i,"realm"))totalCohortsThisCell+=cohortCount;
             repast::Point<int> initialLocation(x,y);
             for (unsigned i=0;i<numCohortGroups;i++){
                 if (E->_Realm==CohortDefinitions::Get()->Trait(i,"realm")){
                     for (unsigned j=0;j<cohortCount;j++){
                         //make sure the agentId is unique on this thread!!
                         // values are int id, int startProc, int agentType, 
                         repast::AgentId id(Cohort::_NextID, rank, _cohortType);
                         id.currentRank(rank);
                         Cohort* c = new Cohort(id);
                         c->setup(i,totalCohortsThisCell, E);
                         _context.addAgent(c);
                         discreteSpace->moveTo(id, initialLocation);
                    }
                }
             }
             for (unsigned i=0;i<numStockGroups;i++){
                 //one stock per functional group - use s to make sure agentId are different
                 if (E->_Realm==StockDefinitions::Get()->Trait(i,"realm")){
                         repast::AgentId id(s, rank, _stockType);
                         id.currentRank(rank);
                         Stock* agent = new Stock(id,E);
                         agent->setup(i, E);
                         _context.addAgent(agent);
                         discreteSpace->moveTo(id, initialLocation);
                         s++;
                }
             }
            }
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
	_props->putProperty("init.time", ss.str());
}
//------------------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------------------
//Run the model
//------------------------------------------------------------------------------------------------------------
void MadModel::step(){
    unsigned CurrentTimeStep=RepastProcess :: instance ()->getScheduleRunner ().currentTick () - 1;
	if(repast::RepastProcess::instance()->rank() == 0) 
        std::cout << " TICK " << CurrentTimeStep << std::endl;

    //make sure all data is synced across threads
    sync();

	//for(auto& a:agents){
	//	a->step(this);
    //    }
/*	for (int i=1;i<3;i++){
    for (auto& a:agents){
 	    a->move(this,-1,-1);
        }
    }
    for(auto& a:agents){
        if (a->getId().startingRank()==0 && a->getId().id()==25)std::cout << " AGENT " << " ON Thread " << repast::RepastProcess::instance()->rank() << std::endl;

		a->reportLocation(this);
        }*/
    int xmin=        discreteSpace->dimensions().origin().getX() ;
    int xmax= xmin + discreteSpace->dimensions().extents().getX();
    int ymin=        discreteSpace->dimensions().origin().getY() ;
    int ymax= ymin + discreteSpace->dimensions().extents().getY();
    for(int x = xmin; x < xmax; x++){
        for(int y = ymin; y < ymax; y++){

            Environment* E=_Env[x-_minX+(_maxX-_minX+1)*(y-_minY)];
            repast::Point<int> location(x,y);
            std::vector<MadAgent*> agentsInCell;
            //query four neighbouring cells, distance 0 (i.e. just the centre cell) - "true" keeps the centre cell.
            repast::VN2DGridQuery<MadAgent> VN2DQuery(space());
            VN2DQuery.query(location, 0, true, agentsInCell);
            std::vector<Stock*> stocks;
            std::vector<Cohort*> cohorts;
            for (auto a:agentsInCell){
                if (a->getId().agentType()==_stockType) stocks.push_back( (Stock*) a); else cohorts.push_back( (Cohort*) a);
            }
            
            double allBiomass=0;
            for (auto s:stocks)allBiomass+=s->_TotalBiomass;
            for (auto s:stocks)s->step(allBiomass,E, CurrentTimeStep);
            for (auto c:cohorts)c->step(E, cohorts, stocks, CurrentTimeStep);
        }
    }
        
        
}
//------------------------------------------------------------------------------------------------------------
void MadModel::sync(){
    //These lines synchronize the agents across all threads - if there is more than one...
/*	discreteSpace->balance();
    repast::RepastProcess::instance()->synchronizeAgentStatus<MadAgent, MadAgentPackage, 
             MadAgentPackageProvider, MadAgentPackageReceiver>(_context, *provider, *receiver, *receiver);
    
    repast::RepastProcess::instance()->synchronizeProjectionInfo<MadAgent, MadAgentPackage, 
             MadAgentPackageProvider, MadAgentPackageReceiver>(_context, *provider, *receiver, *receiver);

	repast::RepastProcess::instance()->synchronizeAgentStates<MadAgentPackage, 
             MadAgentPackageProvider, MadAgentPackageReceiver>(*provider, *receiver);
             */
}
//------------------------------------------------------------------------------------------------------------
// Packages for exchanging agents across threads
//------------------------------------------------------------------------------------------------------------
/*

MadAgentPackageProvider::MadAgentPackageProvider(repast::SharedContext<MadAgent>* agentPtr): agents(agentPtr){ }
//------------------------------------------------------------------------------------------------------------

void MadAgentPackageProvider::providePackage(MadAgent * agent, std::vector<MadAgentPackage>& out){
    repast::AgentId id = agent->getId();
    MadAgentPackage package(id.id(), id.startingRank(), id.agentType(), id.currentRank(), agent->getC(), agent->getTotal());
    out.push_back(package);
}
//------------------------------------------------------------------------------------------------------------

void MadAgentPackageProvider::provideContent(repast::AgentRequest req, std::vector<MadAgentPackage>& out){
    std::vector<repast::AgentId> ids = req.requestedAgents();
    for(size_t i = 0; i < ids.size(); i++){
        providePackage(agents->getAgent(ids[i]), out);
    }
}
//------------------------------------------------------------------------------------------------------------


MadAgentPackageReceiver::MadAgentPackageReceiver(repast::SharedContext<MadAgent>* agentPtr): agents(agentPtr){}
//------------------------------------------------------------------------------------------------------------

MadAgent * MadAgentPackageReceiver::createAgent(MadAgentPackage package){
    repast::AgentId id(package.id, package.rank, package.type, package.currentRank);
    return new MadAgent(id, package.c, package.total);
}
//------------------------------------------------------------------------------------------------------------

void MadAgentPackageReceiver::updateAgent(MadAgentPackage package){
    repast::AgentId id(package.id, package.rank, package.type);
    MadAgent * agent = agents->getAgent(id);
    agent->set(package.currentRank, package.c, package.total);
}

//------------------------------------------------------------------------------------------------------------

*/
