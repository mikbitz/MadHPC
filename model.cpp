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
#include "randomizer.h"
#include "AgentPackage.h"

using namespace std;
using namespace repast;
//arbitrary numbers to distiguish the agents types

int MadModel::_stockType=0, MadModel::_cohortType=1;
//------------------------------------------------------------------------------------------------------------
//Constructor and destructor
//------------------------------------------------------------------------------------------------------------
MadModel::MadModel(std::string propsFile, int argc, char** argv, boost::mpi::communicator* comm): _context(comm){
    //-----------------
    //Pull in parameter from the model.props file
    _props = new repast::Properties(propsFile, argc, argv, comm);
    //Number of timesteps
	_stopAt = repast::strToInt(_props->getProperty("stop.at"));
    //extent of buffer zones in grid units - this many grid cells are shared at the boundary between cores
    int gridBuffer = repast::strToInt(_props->getProperty("grid.buffer"));
    //Grid extent
    _minX=repast::strToInt(_props->getProperty("min.x"));
    _minY=repast::strToInt(_props->getProperty("min.y"));
    _maxX=repast::strToInt(_props->getProperty("max.x"));
    _maxY=repast::strToInt(_props->getProperty("max.y"));
    _dimX=repast::strToInt(_props->getProperty("proc.per.x"));
    _dimY=repast::strToInt(_props->getProperty("proc.per.y"));

    //-----------------
	//create the model grid
    repast::Point<double> origin(_minX,_minY);
    repast::Point<double> extent(_maxX-_minX+1, _maxY-_minY+1);
    
    repast::GridDimensions gd(origin, extent);
    
    std::vector<int> processDims;
    processDims.push_back(_dimX);
    processDims.push_back(_dimY);
  
    discreteSpace = new repast::SharedDiscreteSpace<MadAgent, repast::WrapAroundBorders, repast::SimpleAdder<MadAgent> >("AgentDiscreteSpace", gd, processDims, gridBuffer, comm);
	
   
    //The agent container is a context. Add the grid to it.
   	_context.addProjection(discreteSpace);
    
   //-----------------
    //Set up the cross-thread data transferers
	provider = new MadAgentPackageProvider(&_context);
	receiver = new MadAgentPackageReceiver(&_context);
    //---------------
    //variables to hold totals across all cells
    _totalCohorts=0;
    _totalStocks=0;
    _totalCohortAbundance=0;
    _totalCohortBiomass=0;
    _totalStockBiomass=0;
    _totalOrganciPool=0;
    _totalRespiratoryCO2Pool=0;
    //----------------
    //local grid extents on this thread
    _xlo=        discreteSpace->dimensions().origin().getX() ;
    _xhi= _xlo + discreteSpace->dimensions().extents().getX();
    _ylo=        discreteSpace->dimensions().origin().getY() ;
    _yhi= _ylo + discreteSpace->dimensions().extents().getY();

}
//------------------------------------------------------------------------------------------------------------
MadModel::~MadModel(){
	delete _props;
	delete provider;
	delete receiver;
    for (size_t i = 0; i < dataSets.size(); ++i) {
		delete dataSets[i];
	}

}
//------------------------------------------------------------------------------------------------------------
//Model Initialisation
//------------------------------------------------------------------------------------------------------------
void MadModel::initSchedule(repast::ScheduleRunner& runner){
	runner.scheduleEvent(1, 1, repast::Schedule::FunctorPtr(new repast::MethodFunctor<MadModel> (this, &MadModel::step)));
	runner.scheduleStop(_stopAt);
    runner.scheduleEndEvent(Schedule::FunctorPtr(new MethodFunctor<MadModel> (this, &MadModel::dataSetClose)));

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
    
    randomizer random;
    random.SetSeed(100);
    
    //explicitly use the local bounds of the grid on this thread to create countOfAgents per cell.
    //Not doing this can lead to problems with agents in distant cells not within the local thread neighbourhood
    //see SharedBaseGrid.h moveTo method

    unsigned totalCohorts=0,totalStocks=0;

    unsigned cNum=0,sNum=0;

    unsigned totalStocksThisCell=0;

  
    int s=0;

    for (int x = _xlo; x < _xhi; x++){
        for (int y = _ylo; y < _yhi; y++){
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
                         //agent also needs id of its current thread
                         id.currentRank(rank);
                         Cohort* c = new Cohort(id);
                         c->setup(i,totalCohortsThisCell, E,random);
                         _context.addAgent(c);
                         discreteSpace->moveTo(id, initialLocation);
                         //to get movement right agent needs its own copy of location
                         c->setLocation(x,y);
                         _totalCohorts++;
                         _totalCohortAbundance += c->_CohortAbundance;
                         _totalCohortBiomass += ( c->_IndividualBodyMass + c->_IndividualReproductivePotentialMass ) * c->_CohortAbundance / 1000.;//g to kg
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
                         _totalStockBiomass+=agent->_TotalBiomass/1000; //g to kg
                         s++;
                }
             }
            }
        }
    cout<<"rank "<<rank<<" totalCohorts "<<_totalCohorts<<" totalStocks "<<s<<endl;
    setupOutputs();


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
    _totalCohorts=0;
    _totalStocks=0;
    _totalCohortAbundance=0;
    _totalCohortBiomass=0;
    _totalStockBiomass=0;
    _totalMerged=0;
    _totalReproductions=0;
    _totalDeaths=0;
    _totalMoved=0;
    
    unsigned CurrentTimeStep=RepastProcess :: instance ()->getScheduleRunner ().currentTick () - 1;
    //needed to advance the datalayers to the current timestep
    TimeStep::Get( )->SetMonthly( CurrentTimeStep);
    
	if(repast::RepastProcess::instance()->rank() == 0) 
        std::cout << " TICK " << CurrentTimeStep << std::endl;

    //make sure all data is synced across threads
    sync();
    vector<int> cohortBreakdown,finalCohortBreakdown;
    cohortBreakdown.resize(19);
    finalCohortBreakdown.resize(19);
    for (auto& k:cohortBreakdown){k=0;}
    for (auto& k:finalCohortBreakdown){k=0;}
    //loop over local grid cells

    for(int x = _xlo; x < _xhi; x++){
        for(int y = _ylo; y < _yhi; y++){

            Environment* E=_Env[x-_minX+(_maxX-_minX+1)*(y-_minY)];
            E->zeroPools();

            //store current location in a repast structure for later use
            repast::Point<int> location(x,y);
            
            std::vector<MadAgent*> agentsInCell;
            //query four neighbouring cells, distance 0 (i.e. just the centre cell) - "true" keeps the centre cell.
            repast::VN2DGridQuery<MadAgent> VN2DQuery(space());
            VN2DQuery.query(location, 0, true, agentsInCell);
            std::vector<Stock*> stocks;
            std::vector<Cohort*> cohorts;
            for (auto a:agentsInCell){
               if (a->getId().currentRank()==repast::RepastProcess::instance()->rank() && a->_alive){//agents must be local and living!
                if (a->getId().agentType()==_stockType) stocks.push_back( (Stock*) a); else cohorts.push_back( (Cohort*) a);
               }
            }
            //need to random shuffle cohorts - this uses the version in repast::Random.h
            
            shuffleList<Cohort>(cohorts);
            double allBiomass=0;
            //advance the stocks and cohorts by one timestep
            for (auto s:stocks)allBiomass+=s->_TotalBiomass;
            for (auto s:stocks){s->step(allBiomass,E, CurrentTimeStep);}
            for (auto c:cohorts)c->step(E, cohorts, stocks, CurrentTimeStep);
            for (auto s:stocks){_totalStockBiomass+=s->_TotalBiomass/1000;_totalStocks++;}

            //cohorts can have one offspring per timestep - add the offspring to the model
            vector<Cohort*> newCohorts;
            for (auto c:cohorts)if(c->_newH!=NULL){
                _context.addAgent(c->_newH);
                discreteSpace->moveTo(c->_newH->getId(), location);
                newCohorts.push_back(c->_newH);
                _totalReproductions++;
            }
            //make sure new cohorts can get marked for death/merged - new cohorts are guaranteed to be local
            for (auto n:newCohorts){cohorts.push_back(n);}
            newCohorts.clear();
            for (auto c:cohorts){c->markForDeath();if (!c->_alive)_totalDeaths++;}
            
            _totalMerged+=CohortMerger::MergeToReachThresholdFast(cohorts);

            for (auto c:cohorts){
                if (c->_alive){
                    _totalCohorts++;
                    _totalCohortAbundance += c->_CohortAbundance;
                    _totalCohortBiomass += ( c->_IndividualBodyMass + c->_IndividualReproductivePotentialMass ) * c->_CohortAbundance / 1000.;
                    cohortBreakdown[c->_FunctionalGroupIndex]++;
                }
            }
            //care with sync() here - need to get rid of non-local not-alive agents:currently this is a lazy delete for non-local agents (they get removed one timestep late)?
            for (auto a:agentsInCell)if (!a->_alive)_context.removeAgent(a->getId());//does this delete the agent storage? - yes if Boost:shared_ptr works OK
            _totalOrganciPool+=E->organicPool()/1000;
            _totalRespiratoryCO2Pool+=E->respiratoryPool()/1000;
        }
    }
    //numbers per functional group - here's how to add up across a vector over threads.
    MPI_Reduce(cohortBreakdown.data(), finalCohortBreakdown.data(), 19, MPI::INT, MPI::SUM, 0, MPI_COMM_WORLD);
    
   // if(repast::RepastProcess::instance()->rank() == 0){
   //     for (auto& k:finalCohortBreakdown){cout<<k<<" ";}cout<<endl;cout.flush();}

    //find out which agents need to move
    //_moved has been set to false for new agents
    vector<Cohort*> movers;
    for(int x = _xlo; x < _xhi; x++){
        for(int y = _ylo; y < _yhi; y++){
            Environment* E=_Env[x-_minX+(_maxX-_minX+1)*(y-_minY)];
            repast::Point<int> location(x,y);
            std::vector<MadAgent*> agentsInCell;
            //query four neighbouring cells, distance 0 (i.e. just the centre cell) - "true" keeps the centre cell.
            repast::VN2DGridQuery<MadAgent> VN2DQuery(space());
            VN2DQuery.query(location, 0, true, agentsInCell);

            for (auto a:agentsInCell){
                 if (a->getId().currentRank()==repast::RepastProcess::instance()->rank()){
                     if (a->getId().agentType()==MadModel::_cohortType) {
                         ((Cohort*) a)->moveIt(E,this);
                         if (a->_moved){
                             movers.push_back((Cohort *) a);
                        }
                     }
                }
            }
        }
    }
    _totalMoved=movers.size();

    //RPHC default limited in that it uses a cartesian grid of threads to map onto the model grid
    //and it can't cope if an agent tried to move more than one cartesian grid cell in one go.
    //For some arrangments of cores this causes a crash if maxMoveDist below (number of model grid cells moved this timstep)
    //is set too high. To work around this, some movers have to be displaced in multiple steps. moving across threads multiple times
    //Each time they move we need a sync() to copy their properties (including how far still left to go) and then we need to re-acquire
    //the agents on each thread to take account of tne changes. This is unnecessarily expensive and needs a change to RHPC defaults.
    //The main problem is in SharedBaseGrid::balance() which assumes moves are only to cartesian grid neighbours 
    //see sync() below where this is specified
    //Even using a large buffer will not fix this - an edit of RHPC code is required.
    int maxMoveDist=min((_maxX-_minX)/_dimX,(_maxY-_minY)/_dimY);
    assert(maxMoveDist>1);
    //Each thread has to find out locally if it has moved all agents required
    //globallyFinished is needed to check across all threads to see if any are still active
    //sync() needs to be called by ALL threads if any are still going - otherwise things can hang
    bool globallyFinished=false;
    vector<int>displacement={0,0};
    while (!globallyFinished){
     //are we finished locally?
     bool finished=true;
     for (auto& m:movers){
         displacement[0]=m->_destination[0]-m->_location[0];
         displacement[1]=m->_destination[1]-m->_location[1];
         //move things with less than maxMovedDist - these are then settled and _moved becomes false
         if ( m->_moved &&
             displacement[0] <  maxMoveDist + 1 && displacement[1] <  maxMoveDist + 1  &&
             displacement[0] > -maxMoveDist - 1 && displacement[1] > -maxMoveDist - 1
         ){
         space()->moveTo(m,m->_destination);m->_moved=false;m->_location=m->_destination;
         //now deal with fast movers
         }else{
             //not finished as these cohorts will need to move again
             finished=false;
             //move them by maxMoveDist and then adjust their remaining displacement
             vector<int> d={maxMoveDist,maxMoveDist};
             d[0]=d[0]*( (displacement[0] > 0) - (displacement[0] < 0));
             d[1]=d[1]*( (displacement[1] > 0) - (displacement[1] < 0));
             m->_location[0]=m->_location[0] + d[0];
             m->_location[1]=m->_location[1] + d[1];
             space()->moveByDisplacement(m,d);
         }
     }
     //check across all threads to see if we are finished globally
     MPI_Allreduce(&finished, &globallyFinished, 1, MPI::BOOL, MPI::LAND,MPI_COMM_WORLD);
     if (!globallyFinished){
         //some moves still needed - synchronize the temporary moves
         sync();
         movers.clear();
         //now check on each thread for moved agents that still need to be moved
     	 std::vector<MadAgent*> agents;
         _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
 		 for(auto& a:agents){
             if (a->getId().agentType()==MadModel::_cohortType && a->_moved) { 
               movers.push_back((Cohort*) a);
             }
         }
     }
    }
}


//------------------------------------------------------------------------------------------------------------
void MadModel::sync(){
    //These lines synchronize the agents across all threads - if there is more than one...
    //Question - are threads guaranteed to be in sync?? (i.e. are we sure that all threads are on the same timestep?)
    //Possibilites for sync in terms of where to send agents are POLL, USE_CURRENT, USE_LAST_OR_CURRENT, USE_LAST_OR_POLL
    //USE_CURRENT assumes agents do no move beyond the neighbours in the cartesian grid of threads. This fails
    //for some long distance moves - POLL seems the safest (not sure if guaranteed to work though...also maybe slower)
	discreteSpace->balance();
    repast::RepastProcess::instance()->synchronizeAgentStatus<MadAgent, AgentPackage, 
             MadAgentPackageProvider, MadAgentPackageReceiver>(_context, *provider, *receiver, *receiver,RepastProcess::POLL);
    
    repast::RepastProcess::instance()->synchronizeProjectionInfo<MadAgent, AgentPackage, 
             MadAgentPackageProvider, MadAgentPackageReceiver>(_context, *provider, *receiver, *receiver,RepastProcess::POLL);

	repast::RepastProcess::instance()->synchronizeAgentStates<AgentPackage, 
             MadAgentPackageProvider, MadAgentPackageReceiver>(*provider, *receiver);
             
}
//------------------------------------------------------------------------------------------------------------
// Packages for exchanging agents across threads
//------------------------------------------------------------------------------------------------------------


MadAgentPackageProvider::MadAgentPackageProvider(repast::SharedContext<MadAgent>* agentPtr): agents(agentPtr){ }
//------------------------------------------------------------------------------------------------------------

void MadAgentPackageProvider::providePackage(MadAgent* agent, std::vector<AgentPackage>& out){
    repast::AgentId id = agent->getId();
    if (id.agentType() == MadModel::_cohortType){
     AgentPackage package(id.id(), id.startingRank(), id.agentType(), id.currentRank());
     ((Cohort*)agent)->PushThingsIntoPackage(package);
     out.push_back(package);
    }
}

//------------------------------------------------------------------------------------------------------------

void MadAgentPackageProvider::provideContent(repast::AgentRequest req, std::vector<AgentPackage>& out){
    std::vector<repast::AgentId> ids = req.requestedAgents();
    for(size_t i = 0; i < ids.size(); i++){
        providePackage(agents->getAgent(ids[i]), out);
    }
}
//------------------------------------------------------------------------------------------------------------


MadAgentPackageReceiver::MadAgentPackageReceiver(repast::SharedContext<MadAgent>* agentPtr): agents(agentPtr){}
//------------------------------------------------------------------------------------------------------------

MadAgent * MadAgentPackageReceiver::createAgent(AgentPackage package){
    repast::AgentId id(package._id, package._rank, package._type, package._currentRank);
    if (id.agentType() == MadModel::_cohortType){
        Cohort* c=new Cohort(id,package);
        return c;
    } else {
        return NULL;
    }
}
//------------------------------------------------------------------------------------------------------------
//This function is needed if buffers are being used so that agents can interact across cellSize
//At present it does nothing
void MadAgentPackageReceiver::updateAgent(AgentPackage package){
    repast::AgentId id(package._id, package._rank, package._type);
    if (id.agentType() == MadModel::_cohortType){
      Cohort* agent = (Cohort*)(agents->getAgent(id));//I think this matches irrespective of the value of currentRank (AgentId== operator doesn't use it)
      agent->PullThingsOutofPackage(package);
      //agent->set(id.currentRank(),package;// Do not use !! this line is incorrect!!
    }
}
//------------------------------------------------------------------------------------------------------------
void MadModel::setupOutputs(){
    
	//The things added to the datasetbuilder will be accumulated over cores each timestep and output to data.csv
	SVDataSetBuilder svbuilder("./output/data.csv", ",", repast::RepastProcess::instance()->getScheduleRunner().schedule());

    DispersalSum* DSum = new DispersalSum(this);
    svbuilder.addDataSource(repast::createSVDataSource("Dispersals", DSum, std::plus<int>()));
    
    ExtinctionSum* ESum = new ExtinctionSum(this);
    svbuilder.addDataSource(repast::createSVDataSource("Extinctions", ESum, std::plus<int>()));
    
    ProductionSum* PSum = new ProductionSum(this);
    svbuilder.addDataSource(repast::createSVDataSource("Productions", PSum, std::plus<int>()));
    
    CombinationSum* coSum = new CombinationSum(this);
    svbuilder.addDataSource(repast::createSVDataSource("Combinations", coSum, std::plus<int>()));
    
    CohortSum* cSum = new CohortSum(this);
	svbuilder.addDataSource(repast::createSVDataSource("Total Cohorts", cSum, std::plus<int>()));
    
    StockSum* sSum = new StockSum(this);
	svbuilder.addDataSource(repast::createSVDataSource("Total Stocks", sSum, std::plus<int>()));
    
    CohortAbundanceSum* caSum = new CohortAbundanceSum(this);
  	svbuilder.addDataSource(repast::createSVDataSource("Total Cohort Abundance", caSum, std::plus<double>()));
        
    CohortOrganicPool* cpool = new CohortOrganicPool(this);
  	svbuilder.addDataSource(repast::createSVDataSource("OrganicPool", cpool, std::plus<double>()));
    
    CohortResp* rpool = new CohortResp(this);
  	svbuilder.addDataSource(repast::createSVDataSource("RespiratoryC02Pool", rpool, std::plus<double>()));
    
    StockBiomassSum* sbSum = new StockBiomassSum(this);
  	svbuilder.addDataSource(repast::createSVDataSource("Total Stock Biomass", sbSum, std::plus<double>()));
    
    CohortBiomassSum* cbSum = new CohortBiomassSum(this);
  	svbuilder.addDataSource(repast::createSVDataSource("Total Cohort Biomass", cbSum, std::plus<double>()));
    

	addDataSet(svbuilder.createDataSet());

//  Netcdf is supposed to act in a similar way - currently not functional - example is left over from the RHPC zombie demo code
//	NCDataSetBuilder builder("./output/data.ncf", RepastProcess::instance()->getScheduleRunner().schedule());
//	InfectionSum* infectionSum = new InfectionSum(this);
//	builder.addDataSource(repast::createNCDataSource("number_infected", infectionSum, std::plus<int>()));
//	addDataSet(builder.createDataSet());

}
//------------------------------------------------------------------------------------------------------------

void MadModel::dataSetClose() {
	for (size_t i = 0; i < dataSets.size(); ++i) {
		(dataSets[i])->write();
		(dataSets[i])->close();
	}
}
//---------------------------------------------------------------------------------------------------------------------------
void MadModel::addDataSet(repast::DataSet* dataSet) {
	dataSets.push_back(dataSet);
	ScheduleRunner& runner = RepastProcess::instance()->getScheduleRunner();
	runner.scheduleEvent(0.1, 1, Schedule::FunctorPtr(new MethodFunctor<repast::DataSet> (dataSet,&repast::DataSet::record)));
	Schedule::FunctorPtr dsWrite = Schedule::FunctorPtr(new MethodFunctor<repast::DataSet> (dataSet,&repast::DataSet::write));
    //output every 100 steps
	runner.scheduleEvent(100.2, 100, dsWrite);
}
//---------------------------------------------------------------------------------------------------------------------------
void MadModel::tests(){
    int rank = repast::RepastProcess::instance()->rank();
    randomizer random;
    random.SetSeed(100);
    int nranks=_dimX*_dimY;
    
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
    int x=_xlo,y=_ylo;
    repast::Point<int> initialLocation(x,y);
    repast::Point<int> origin(_minX,_minY);
    Environment* E=_Env[x-_minX+(_maxX-_minX+1)*(y-_minY)];
    
    //---------------------------------------------------
    //***-------------------TEST 1-------------------***//
    //---------------------------------------------------
    //initially all agents are just on rank 0
   int localTotals=0,globalTotals=0,n=1000;
    //create some agents on rank 0
    if (rank==0){
     cout<<"Test1: create n="<<n<<" agents on rank 0"<<endl;
     for (int i=0;i<n;i++){
      repast::AgentId id(Cohort::_NextID, rank, _cohortType);
      id.currentRank(rank);
      Cohort* c = new Cohort(id);
      c->setup(0,1, E,random);
      _context.addAgent(c);
      discreteSpace->moveTo(id, initialLocation);
     }
    }

    std::vector<MadAgent*> agents;
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    localTotals=agents.size();
    
    //check across all threads to see if agents all still exist
    
    MPI_Allreduce(&localTotals, &globalTotals, 1, MPI::INT, MPI::SUM,MPI_COMM_WORLD);
    assert(globalTotals==n);
    if (rank==0)cout<<"Test1 succeeded"<<endl;
    
    //---------------------------------------------------
    //***-------------------TEST 2-------------------***//
    //---------------------------------------------------
    //remove agents, but again only on rank 0
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    if (rank==0){
     cout<<"Test2: delete all "<<n<<" agents on rank 0"<<endl;
     for (auto a:agents)_context.removeAgent(a->getId());
    }
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    localTotals=agents.size();
    MPI_Allreduce(&localTotals, &globalTotals, 1, MPI::INT, MPI::SUM,MPI_COMM_WORLD);
    assert(globalTotals==0);
    if (rank==0)cout<<"Test2 succeeded "<<endl;

    //---------------------------------------------------
    //***-------------------TEST 3-------------------***//
    //---------------------------------------------------
    //add agents back to rank 0
    if (rank==0){
     cout<<"Test3: re-create "<<n<<" agents on rank 0"<<endl;
     for (int i=0;i<n;i++){
      repast::AgentId id(Cohort::_NextID, rank, _cohortType);
      id.currentRank(rank);
      Cohort* c = new Cohort(id);
      c->setup(0,1, E,random);
      _context.addAgent(c);
      discreteSpace->moveTo(id, initialLocation);
     }
    }
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    localTotals=agents.size();
    MPI_Allreduce(&localTotals, &globalTotals, 1, MPI::INT, MPI::SUM,MPI_COMM_WORLD);
    assert(globalTotals==n);
    if (rank==0)cout<<"Test3 succeeded "<<endl;
    
    //---------------------------------------------------
    //***-------------------TEST 4-------------------***//
    //---------------------------------------------------
    //all agents still exist: allow any thread to move them
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    if (rank==0)cout<<"Test4: move agents by 1 unit in x "<<endl;
    std::vector<int> displ{1,0};
    for (auto a:agents)space()->moveByDisplacement(a,displ);
    sync();
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    localTotals=agents.size();
    MPI_Allreduce(&localTotals, &globalTotals, 1, MPI::INT, MPI::SUM,MPI_COMM_WORLD);
    assert(globalTotals==n);
    vector<int> location;
    for (auto a:agents){
      space()->getLocation(a->getId(), location);
      assert(location[0]==_xlo+1 && location[1]==_ylo+0);
    }
    if (rank==0)cout<<"Test4 succeeded "<<endl;
    //---------------------------------------------------
    //***-------------------TEST 5-------------------***//
    //---------------------------------------------------
    //all agents still exist: allow any thread to move them
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    if (rank==0)cout<<"Test5: move agents by 1 unit in y "<<endl;
    std::vector<int> displ2{0,1};
    for (auto a:agents)space()->moveByDisplacement(a,displ2);
    sync();
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    localTotals=agents.size();
    MPI_Allreduce(&localTotals, &globalTotals, 1, MPI::INT, MPI::SUM,MPI_COMM_WORLD);
    assert(globalTotals==n);
    
    for (auto a:agents){
      space()->getLocation(a->getId(), location);
      assert(location[0]==_xlo+1 && location[1]==_ylo+1);
    }
    if (rank==0)cout<<"Test5 succeeded "<<endl;
    //---------------------------------------------------
    //***-------------------TEST 6-------------------***//
    //---------------------------------------------------
    //all agents still exist: allow any thread to move them
    //note grid is wrapped, so moves of any distance should work.
    //all agents are still co-located, so should end up on the same thread
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    int xm1=200,ym1=-15;
    if (rank==0)cout<<"Test6: move agents by "<<xm1<<" units in x and "<<ym1<<" units in y "<<endl;
    std::vector<int> displ3{xm1,ym1};
    for (auto a:agents)space()->moveByDisplacement(a,displ3);
    sync();
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    localTotals=agents.size();
    MPI_Allreduce(&localTotals, &globalTotals, 1, MPI::INT, MPI::SUM,MPI_COMM_WORLD);
    assert(globalTotals==n);
    for (auto a:agents){
      space()->getLocation(a->getId(), location);
      //calculate the wrapped position including previous move by +1,+1 
      int xpos=xm1+1+_minX,ypos=ym1+1+_minY;
      while (xpos<_minX)xpos+=(_maxX - _minX + 1);
      while (ypos<_minY)ypos+=(_maxY - _minY + 1);
      xpos=xpos % (_maxX - _minX + 1);
      ypos=ypos % (_maxY - _minY + 1);
      assert(location[0]==xpos && location[1]==ypos);
    }
    if (rank==0)cout<<"Test6 succeeded "<<endl;
    //---------------------------------------------------
    //***-------------------TEST 7-------------------***//
    //---------------------------------------------------
    //all agents still exist: allow any thread to move them
    //note grid is wrapped, so moves of any distance should work.
    //all agents are still co-located, so should end up on the same thread
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    int xm2=-79,ym2=3355;
    if (rank==0)cout<<"Test7: move agents by "<<xm2<<" units in x and "<<ym2<<" units in y "<<endl;
    std::vector<int> displ4{xm2,ym2};
    for (auto a:agents)space()->moveByDisplacement(a,displ4);
    sync();
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    localTotals=agents.size();
    MPI_Allreduce(&localTotals, &globalTotals, 1, MPI::INT, MPI::SUM,MPI_COMM_WORLD);
    assert(globalTotals==n);
    for (auto a:agents){
      space()->getLocation(a->getId(), location);
      //calculate the wrapped position including previous move by xm1+1,ym1+1 
      int xpos=xm1+xm2+1+_minX,ypos=ym1+ym2+1+_minY;
      while (xpos<_minX)xpos+=(_maxX - _minX + 1);
      while (ypos<_minY)ypos+=(_maxY - _minY + 1);
      xpos=xpos % (_maxX - _minX + 1);
      ypos=ypos % (_maxY - _minY + 1);
      assert(location[0]==xpos && location[1]==ypos);
    }
    //---------------------------------------------------
    //***-------------------TEST 8-------------------***//
    //---------------------------------------------------
    //move all agents back to origin (thread 0)
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);

    if (rank==0)cout<<"Test8: move agents back to origin: "<<_minX<<" "<<_minY<<endl;
    for (auto a:agents)space()->moveTo(a,origin);
    sync();
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    localTotals=agents.size();
    MPI_Allreduce(&localTotals, &globalTotals, 1, MPI::INT, MPI::SUM,MPI_COMM_WORLD);
    assert(globalTotals==n);
    for (auto a:agents){
      space()->getLocation(a->getId(), location);
      assert(location[0]==_minX && location[1]==_minY);
    }
    if (rank==0)cout<<"Test8 succeeded "<<endl;
    //---------------------------------------------------
    //***-------------------TEST 9-------------------***//
    //---------------------------------------------------
    //move all agents to random locations
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);

    if (rank==0)cout<<"Test9: move agents to random location: "<<endl;
    
    for (auto a:agents){
        vector<int> loc{int(random.GetUniform()*1000.),int(random.GetUniform()*1000.)};
        space()->moveTo(a,loc);
    }
    sync();
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    localTotals=agents.size();
    cout<<"Total on rank:"<<rank<<" "<<localTotals<<endl;
    MPI_Allreduce(&localTotals, &globalTotals, 1, MPI::INT, MPI::SUM,MPI_COMM_WORLD);
    assert(globalTotals==n);
    for (auto a:agents){
      space()->getLocation(a->getId(), location);
      assert(location[0]>=_xlo && location[1]>=_ylo && location[0]<_xhi && location[1]<_yhi);
    }
    if (rank==0)cout<<"Test9 succeeded "<<endl;
    //---------------------------------------------------
    //***-------------------TEST 10-------------------***//
    //---------------------------------------------------
    //test add and delete on all threads, not just rank 0
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);

    if (rank==0)cout<<"Test10: delete up to 10 agents per thread, add up to 100 new, move all randomly: "<<endl;
    int nr=0,globalAdded,globalRemoved;
    for (auto a:agents){
        if (nr<int(random.GetUniform()*10.)){a->_alive=false;_context.removeAgent(a->getId());nr++;}
        if (a->_alive){
            vector<int> loc{int(random.GetUniform()*781.),int(random.GetUniform()*500.-250)};
            space()->moveTo(a,loc);
        }
    }
    int nnew=int(random.GetUniform()*100.);
    for (int i=0;i<nnew;i++){
      repast::AgentId id(Cohort::_NextID, rank, _cohortType);
      id.currentRank(rank);
      Cohort* c = new Cohort(id);
      c->setup(0,1, E,random);
      _context.addAgent(c);
      vector<int> loc{int(random.GetUniform()*37.-99),int(random.GetUniform()*188.)};
      space()->moveTo(c,loc);
    }
    
    sync();
    MPI_Allreduce(&nnew, &globalAdded, 1, MPI::INT, MPI::SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&nr, &globalRemoved, 1, MPI::INT, MPI::SUM,MPI_COMM_WORLD);

    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    localTotals=agents.size();
    MPI_Allreduce(&localTotals, &globalTotals, 1, MPI::INT, MPI::SUM,MPI_COMM_WORLD);
    assert(globalTotals==n+globalAdded-globalRemoved);
    for (auto a:agents){
      space()->getLocation(a->getId(), location);
      assert(location[0]>=_xlo && location[1]>=_ylo && location[0]<_xhi && location[1]<_yhi);
      assert(a->getId().currentRank()==rank);
    }
    if (rank==0)cout<<"Test10 succeeded: "<<"removed:"<<globalRemoved<<" added:"<<globalAdded<<endl;
    cout.flush();
    sync();
}
