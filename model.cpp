/* Model.cpp */

#include <stdio.h>
#include <vector>
#include <sstream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/mpi.hpp>
#include "repast_hpc/AgentId.h"
#include "repast_hpc/RepastProcess.h"
#include "repast_hpc/Utilities.h"
#include "repast_hpc/Properties.h"
#include "repast_hpc/SharedBaseGrid.h"                      // VN2D broken without this!
#include "repast_hpc/VN2DGridQuery.h"
#include "repast_hpc/Moore2DGridQuery.h"
#include "repast_hpc/initialize_random.h"
#include "repast_hpc/Point.h"
#include "repast_hpc/AgentRequest.h"
#include "repast_hpc/SVDataSet.h"
#include "repast_hpc/SVDataSetBuilder.h"
#include "repast_hpc/Random.h"

#ifndef _WIN32
#include "repast_hpc/NCDataSetBuilder.h"
#endif

#include "model.h"
#include "agent.h"
#include "CohortMerger.h"
#include "Environment.h"
#include "Parameters.h"
#include "Groups.h"
#include "Stock.h"
#include "Cohort.h"
#include "FileReader.h"
#include "CohortSum.h"
#include "TimeStep.h"
#include "Constants.h"
#include "randomizer.h"
#include "RandomRepast.h"
#include "AgentPackage.h"
#include <netcdf>

//repast shuffleList only works on pointers
//this version is for vectors of unsigned

void shuffleList(std::vector<unsigned>& elementList){
  if(elementList.size() <= 1) return;
  repast::IntUniformGenerator rnd = repast::Random::instance()->createUniIntGenerator(0, elementList.size() - 1);
  unsigned swap;
  for(size_t i = 0, sz = elementList.size(); i < sz; i++){
    int other = rnd.next();
    swap = elementList[i];
    elementList[i] = elementList[other];
    elementList[other] = swap;
  }
}

using namespace std;
using namespace repast;
//arbitrary numbers to distiguish the agents types

int MadModel::_stockType=0, MadModel::_cohortType=1, MadModel::_humanType=2;
//------------------------------------------------------------------------------------------------------------
//Constructor and destructor
//------------------------------------------------------------------------------------------------------------
MadModel::MadModel(repast::Properties& props,  boost::mpi::communicator* comm): _context(comm){
    //-----------------
    //Pull in parameter from the model.props file
    _props = &props;
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
    _noLongitudeWrap=repast::strToInt(_props->getProperty("noLongitudeWrap"));
    _dispersalSelection=_props->getProperty("simulation.DispersalSelection");
    _randomSeed=repast::strToUInt(_props->getProperty("global.random.seed"));
    //-----------------
    //naming convention for output files
    if(repast::RepastProcess::instance()->rank() == 0){
       _filePrefix=                  _props->getProperty("experiment.output.directory")+
                      "/experiment."+_props->getProperty("experiment.name");
       if (!boost::filesystem::exists(_filePrefix))boost::filesystem::create_directories(_filePrefix);
       std::string runNumber= _props->getProperty("run.number");
       std::string m00="/run_";
       if (runNumber!=""){
           m00=m00+runNumber;
       }else{
         //auto-increment run number if run.number is not set
         int i=0;
         m00="/run_000",runNumber="000";
         std::string pfx="00";
         while(boost::filesystem::exists(_filePrefix+m00)){    // Find a new directory name
           i++;
           std::stringstream ss;
           if (i>9 ) pfx="0";
           if (i>99) pfx="";
           ss<<pfx<<i;
           runNumber=ss.str();
           m00="/run_"+runNumber;
         }
       }
       if (!boost::filesystem::exists(_filePrefix+m00))boost::filesystem::create_directories(_filePrefix+m00);
       _props->putProperty ("run.number",runNumber);
       _filePrefix= _filePrefix+m00+"/";
       _filePostfix="";
       cout<<"Outputfiles will be named "<<_filePrefix<<"<Data Name>"<<_filePostfix<<".<filenameExtension>"<<endl;
    }
    //-----------------
	//create the model grid
    repast::Point<double> origin(_minX,_minY);
    repast::Point<double> extent(_maxX-_minX+1, _maxY-_minY+1);
    
    repast::GridDimensions gd(origin, extent);
    
    std::vector<int> processDims;
    processDims.push_back(_dimX);
    processDims.push_back(_dimY);
  
    //because RHPC uses templates for grid wrapping, if you want any wrapping at all , you have to 
    //use a space that is wrapped in both x and y - to use a non-wrapped space implies templating all occurrences of
    //"model" - far too much like hard work; either that or we need a wrapper class for spaces, or different models (with classes per space)
    discreteSpace = new wrappedSpaceType("AgentDiscreteSpace", gd, processDims, gridBuffer, comm);
	
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
    
    _totalCohorts=0;
    _totalStocks=0;
    _totalCohortAbundance=0;
    _totalCohortBiomass=0;
    _totalStockBiomass=0;
    _totalMerged=0;
    _totalReproductions=0;
    _totalDeaths=0;
    _totalMoved=0;

    //time the initialisation
    repast::Timer initTimer;
    initTimer.start();
    //rank (i.e. the number of this thread) will be needed to make agents unique *between* threads
    int rank = repast::RepastProcess::instance()->rank();
    
    //get the environmental data - this is stored in the background as a DataLayerSet
    //FileReader F;
    //F.ReadFiles();
    
    //now set up the environmental cells - note at present this uses the full grid, not just local to this thread
    //so that off-thread environment can be easily queried. Currently some duplication here, but it is not a huge amount of data.
    _Env.resize( (_maxX-_minX+1) * (_maxY-_minY+1) );
    for(int x=_minX;x<=_maxX;x++) {
        for(int y=_minY;y<=_maxY;y++) {
            Environment* E=new Environment(x,y);
            _Env[x-_minX+(_maxX-_minX+1)*(y-_minY)]=E;
        }
    }
    //set up the static (i.e. shared) parameters for the Cohorts
    Cohort::setParameters(_props);

    //get the definitions of stocks and cohorts
    StockDefinitions::Initialise(_props->getProperty("input.DataDirectory")+"/"+_props->getProperty("input.StockDefinitionsFileName"));
    CohortDefinitions::Initialise(_props->getProperty("input.DataDirectory")+"/"+_props->getProperty("input.CohortDefinitionsFileName"));
  
    unsigned numCohortGroups=CohortDefinitions::Get()->size();
    unsigned cohortCount = strToInt(_props->getProperty("cohort.count"));
    unsigned numStockGroups = StockDefinitions::Get()->size();
    
         //Values only reduced on thread 0
    if(repast::RepastProcess::instance()->rank() == 0){
     _FinalCohortBiomassMap.resize( (_maxX-_minX+1) * (_maxY-_minY+1) );
     _FinalCohortAbundanceMap.resize( (_maxX-_minX+1) * (_maxY-_minY+1) );
     _FinalStockBiomassMap.resize( (_maxX-_minX+1) * (_maxY-_minY+1) );
     _FinalCohortBreakdown.resize(numCohortGroups);
    }
    
    
    randomizer* random=new RandomRepast;
    //random->SetSeed(100); seed is set from model.props file - see main.cpp

 
    //explicitly use the local bounds of the grid on this thread to create countOfAgents per cell.
    //Not doing this can lead to problems with agents in distant cells not within the local thread neighbourhood
    //see SharedBaseGrid.h moveTo method
    //although latest updates I have made to RHPC should have fixed this...
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
                         //stock needs a location to get herbivory right when cross cell interaction is allowed
                         agent->setLocation(x,y);
                         _totalStockBiomass+=agent->_TotalBiomass/1000; //g to kg
                         s++;
                }
             }
            }
        }
    if(_props->getProperty("verbose")=="true")cout<<"rank "<<rank<<" totalCohorts "<<_totalCohorts<<" totalStocks "<<s<<endl;
    setupOutputs();
    if (rank==0)setupNcOutput();

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
    int rank=repast::RepastProcess::instance()->rank();
    unsigned CurrentTimeStep=RepastProcess :: instance ()->getScheduleRunner ().currentTick () - 1;
    //needed to advance the datalayers to the current timestep
    TimeStep::Get( )->SetMonthly( CurrentTimeStep);
    
	if(rank == 0) 
        std::cout << " TICK " << CurrentTimeStep << std::endl;

    //make sure all data is synced across threads
    sync();
    //vectors length CohortDefinitions::Get()->size() initialized to 0
    vector<int> cohortBreakdown(CohortDefinitions::Get()->size(),0);
    //spatial distributions - for efficiency these should just the the local part, but easier initially to just keep the whole thing on each thread
    //with zeros for non-local, and then reduce onto thread 0
    vector<double> cohortBiomassMap( (_maxX-_minX+1) * (_maxY-_minY+1),0.0 ),cohortAbundanceMap( (_maxX-_minX+1) * (_maxY-_minY+1),0.0 );
    vector<double> stockBiomassMap( (_maxX-_minX+1) * (_maxY-_minY+1),0.0 );


    //loop over local grid cells
    //For cell interaction, we need to make sure copies are updated appropriately.
    //Since we use a moore neighbourhood, blocks of nine cells can be treated independently
    //when the centre cell of the nine is updated - so split the grid into 9 sets of
    //3x3 cells, spearate by a 3 cell stride, update date these, then sync to update the shared data in overlaps
    //and then move to the next block. Randomize the starting point on each time step to ensure ceratin cells do not
    //get privileged.
    vector<unsigned> initialCells(9,0);
    for(unsigned i=0;i<9;i++){
        initialCells[i]=i%3;
        //use cell coords 0,0 0,1 0,2; 1,0 1,1 1,2 etc.
        if (i>2)initialCells[i]+=(_maxX-_minX+1);
        if (i>5)initialCells[i]+=(_maxX-_minX+1);
    }
    
    //cells need to have the same random sequence on each thread (otherwise cross-thread interactions won't preserve sequences), 
    //so restart the rng: Use the timestep to alter the random sequence on successive timesteps
    repast::Random::initialize(_randomSeed + CurrentTimeStep);
    shuffleList(initialCells);//use the version defined above
    
    int buffer=repast::strToInt(_props->getProperty("grid.buffer"));
    int range=0;
    if (buffer==2) range=1;
    
    for (unsigned i=0;i<9;i++){
        int initialCell=initialCells[i];
        int x0=initialCell%(_maxX-_minX+1)+_minX; 
        int y0=initialCell/(_maxX-_minX+1)+_minY;
        while (x0<_xlo)x0+=3;
        while (y0<_ylo)y0+=3;
    for(int y = y0-buffer; y < _yhi+buffer; y+=3){  
     for(int xt = x0-buffer; xt < _xhi+buffer; xt+=3){
        //pick out the cells that are on this thread plus those in the overlap region
        //overlap cells have to run updates in the same order on each thread so that dynamics is
        //reproduced on independent threads
        //NB poles are currently a potential problem as Repast thinks we are on a torus! Must exclude overlaps at max and min latitudes.
        //In practice the internal locations of Cohorts should make this happen as distance measures will be large (x distances get wrapped, but not y, see Cohort.cpp). Stocks need locations though.
        int x=xt;

        if (!_noLongitudeWrap){
           if (xt  < _minX){x = x + (_maxX - _minX + 1);}
           if (xt  > _maxX){x = x - (_maxX - _minX + 1);}
        }

        if (y >= _minY && y <= _maxY){
         if (x >=_minX && x <= _maxX) {
          int cellIndex=x-_minX+(_maxX-_minX+1)*(y-_minY);

 //       if (  ( (x >= (_xlo - buffer) && x < (_xhi + buffer) ) || ( (x>=_xlo - buffer + _maxX - _minX || x < _xhi+buffer - _maxX + _minX) && !_noLongitudeWrap)  ) && 
 //                y >= (_ylo - buffer) && y < (_yhi + buffer) ){

            Environment* E=_Env[cellIndex];
            E->zeroPools();

            //store current location in a repast structure for later use
            repast::Point<int> location(x,y);
            
            std::vector<MadAgent*> agentsInCell,thingsToEat;
            //query four neighbouring cells, distance 0 (i.e. just the centre cell) - "true" keeps the centre cell.
            //these will be the cohorts that act in the current cell
            repast::VN2DGridQuery<MadAgent> VN2DQuery(space());
            VN2DQuery.query(location, 0, true, agentsInCell);
            //things that can be eaten by the agents above - they can be in any of the 
            //eight neighbouring cells (range =1 - requires grid.buffer=2) plus the current cell, or just the current cell (range=0, grid,buffer=0)
            repast::Moore2DGridQuery<MadAgent> Moore2DQuery(space());
            Moore2DQuery.query(location, range, true, thingsToEat);
            std::vector<Stock*> stocks,stocksToEat;
            std::vector<Cohort*> cohorts,cohortsToEat;

            for (auto a:agentsInCell){
               if ( a->_alive){//agents must be living but might not be local!
                //separate out stocks (currently all plants) from cohorts(animals)
                if (a->getId().agentType()==_stockType) stocks.push_back( (Stock*) a); else cohorts.push_back( (Cohort*) a);
               }
            }
            for (auto a:thingsToEat){
               if ( a->_alive){//agents must be living but might not be local!
                //separate out stocks (currently all plants) from cohorts(animals)
                if (a->getId().agentType()==_stockType) stocksToEat.push_back( (Stock*) a); else cohortsToEat.push_back( (Cohort*) a);
               }
            }
            //sort the lists of cohorts using their unique agent number -easiest to do this in place with a lambda
            //this ensures that the subsequent random shuffle can be preserved across threads
            sort(cohorts.begin(), cohorts.end(), [](const Cohort* c1, const Cohort* c2) -> bool { return c1->getId().id() > c2->getId().id();});
            //restart the rng yet again to be sure the shuffling will be the same on different threads
            //this is only here in order to make sure copies of cohorts from remote threads are handled the same in grid overlap regions. 
            repast::Random::initialize(_randomSeed + CurrentTimeStep + cellIndex);

            //need to random shuffle cohorts - so that the same cohort doesn't get to go first on every time step
            //this uses the version in repast::Random.h
            shuffleList<Cohort>(cohorts);
            double allBiomass=0;
            //advance the stocks and cohorts by one timestep
            for (auto s:stocks)allBiomass+=s->_TotalBiomass;
            for (auto s:stocks){s->step(allBiomass,E, CurrentTimeStep);}
            for (auto c:cohorts)c->step(E, cohortsToEat, stocksToEat, CurrentTimeStep,this);
        }
      }
     }
    }
    //each of the independent blocks is complete - now sync before moving to the next set of blocks, if required
    //it may be that a full sync is not really needed here - just synchronizeAgentStates, maybe?
     //if (buffer==2)sync();
     if (buffer==2)repast::RepastProcess::instance()->synchronizeAgentStates<AgentPackage, 
             MadAgentPackageProvider, MadAgentPackageReceiver>(*provider, *receiver);
    }

    //updates of offspring and mergers/death happen after all cells have updated
    //need to keep this separate if there is cross-cell interaction
    for(int y = _ylo; y < _yhi; y++){  
     for(int x = _xlo; x < _xhi; x++){
      int cellIndex=x-_minX+(_maxX-_minX+1)*(y-_minY);
            
            Environment* E=_Env[cellIndex];

            //store current location in a repast structure for later use
            repast::Point<int> location(x,y);
            
            std::vector<MadAgent*> agentsInCell;
            //query four neighbouring cells, distance 0 (i.e. just the centre cell) - "true" keeps the centre cell.
            repast::VN2DGridQuery<MadAgent> VN2DQuery(space());
            VN2DQuery.query(location, 0, true, agentsInCell);
            std::vector<Stock*> stocks;
            std::vector<Cohort*> cohorts;

            for (auto a:agentsInCell){
               if (a->getId().currentRank()==rank && a->_alive){//agents must be local and living!
                //separate out stocks (currently all plants) from cohorts(animals)
                if (a->getId().agentType()==_stockType) stocks.push_back( (Stock*) a); else cohorts.push_back( (Cohort*) a);
               }
            }
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
            
            //set up output stock map
            for (auto s:stocks){_totalStockBiomass+=s->_TotalBiomass/1000;stockBiomassMap[cellIndex]+=s->_TotalBiomass/1000/E->Area();_totalStocks++;}
            //acumulate other totals and spatial maps
            for (auto c:cohorts){
                if (c->_alive){
                    _totalCohorts++;
                    cohortAbundanceMap[cellIndex]+= c->_CohortAbundance/E->Area();
                    _totalCohortAbundance += c->_CohortAbundance;
                    cohortBiomassMap[cellIndex]+=( c->_IndividualBodyMass + c->_IndividualReproductivePotentialMass ) * c->_CohortAbundance / 1000.;//convert to kg
                    _totalCohortBiomass += cohortBiomassMap[cellIndex];
                    cohortBreakdown[c->_FunctionalGroupIndex]++;
                }
            }
            cohortBiomassMap[cellIndex]=cohortBiomassMap[cellIndex]/E->Area();//per square kilometre

            //care with sync() here - need to get rid of not-alive agents:currently this is a lazy delete for new/non-local agents (they get removed one timestep late)?
            for (auto a:agentsInCell)if (!a->_alive)_context.removeAgent(a->getId());//does this delete the agent storage? - yes if Boost:shared_ptr works OK
            _totalOrganciPool+=E->organicPool()/1000;
            _totalRespiratoryCO2Pool+=E->respiratoryPool()/1000;
        }
    }
    
    //find out which agents need to move
    //_moved has been set to false for new agents
    //NB this has to happen after above updates to individual Cohorts (otherwise some cells could get mixed before other have updated, so some cohorts could get updated twice)
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
                 //dispersers must be local and alive!
                 if (a->getId().currentRank()==repast::RepastProcess::instance()->rank() && a->_alive){
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

    vector<int>newPlace={0,0};
    for (auto& m:movers){
        newPlace[0]=int(m->_destination[0]);
        newPlace[1]=int(m->_destination[1]);
        //move things - these are then settled and _moved becomes false
        if ( m->_moved  ){
            space()->moveTo(m,newPlace);m->_moved=false;m->_location=m->_destination;
         }
    }
    //map outputs/vectors not yet available via svbuilder or netcdf builder.  
    //numbers per functional group - here's how to add up across a vector over threads.
    MPI_Reduce(cohortBreakdown.data(), _FinalCohortBreakdown.data(), CohortDefinitions::Get()->size(), MPI::INT, MPI::SUM, 0, MPI_COMM_WORLD);

    //also get the biomass maps
    MPI_Reduce(cohortBiomassMap.data(), _FinalCohortBiomassMap.data(), (_maxX-_minX+1) * (_maxY-_minY+1), MPI::DOUBLE, MPI::SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(cohortAbundanceMap.data(), _FinalCohortAbundanceMap.data(), (_maxX-_minX+1) * (_maxY-_minY+1), MPI::DOUBLE, MPI::SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(stockBiomassMap.data(), _FinalStockBiomassMap.data(), (_maxX-_minX+1) * (_maxY-_minY+1), MPI::DOUBLE, MPI::SUM, 0, MPI_COMM_WORLD);

    //if(repast::RepastProcess::instance()->rank() == 0){asciiOutput(CurrentTimeStep);}
    if(repast::RepastProcess::instance()->rank() == 0){netcdfOutput( CurrentTimeStep );}
}


//------------------------------------------------------------------------------------------------------------
void MadModel::sync(){
    //These lines synchronize the agents across all threads - if there is more than one...
    //Question - are threads guaranteed to be in sync?? (i.e. are we sure that all threads are on the same timestep?)
    //Possibilites for sync in terms of where to send agents are POLL, USE_CURRENT, USE_LAST_OR_CURRENT, USE_LAST_OR_POLL
    //USE_CURRENT assumes agents do not move beyond the neighbours in the cartesian grid of threads. This fails
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
     AgentPackage package(id);
     ((Cohort*)agent)->PushThingsIntoPackage(package);
     out.push_back(package);
    }
    if (id.agentType() == MadModel::_stockType){
     AgentPackage package(id);
     ((Stock*)agent)->PushThingsIntoPackage(package);
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
    repast::AgentId id=package.getId();
    if (id.agentType() == MadModel::_cohortType){
        Cohort* c=new Cohort(id,package);
        return c;
    } else {
        Stock* s=new Stock(id,package);
        return s;
    }
}
//------------------------------------------------------------------------------------------------------------
//This function is needed if buffers are being used so that agents can interact across cells
void MadAgentPackageReceiver::updateAgent(AgentPackage package){
    repast::AgentId id=package.getId();
    if (id.agentType() == MadModel::_cohortType){
      Cohort* agent = (Cohort*)(agents->getAgent(id));//I think this matches irrespective of the value of currentRank (AgentId== operator doesn't use it)
      agent->PullThingsOutofPackage(package);
    }
    if (id.agentType() == MadModel::_stockType){
      Stock* agent = (Stock*)(agents->getAgent(id));
      agent->PullThingsOutofPackage(package);
    }
}
//------------------------------------------------------------------------------------------------------------
void MadModel::setupOutputs(){
    
	//The things added to the datasetbuilder will be accumulated over cores each timestep and output to file filename

        
    std::string filename = _filePrefix+"global.outputs"+_filePostfix+".csv";        
                
	SVDataSetBuilder svbuilder(filename, ",", repast::RepastProcess::instance()->getScheduleRunner().schedule());

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
//----------------------------------------------------------------------------------------------
void MadModel::asciiOutput( unsigned step ) {
    stringstream A,B,C;
    string fAname,fBname,fCname;
    A<<"./output/cohortAbundance_SqKm_"<<(step+1)<<".asc";
    B<<"./output/cohortBiomassKg_SqKm_"<<(step+1)<<".asc";
    C<<"./output/stockBiomassKg_SqKm_"<<(step+1)<<".asc";
    A>>fAname;
    B>>fBname;
    C>>fCname;
    ofstream Abundout(fAname.c_str());
    Abundout.precision(15);
    ofstream Biomaout(fBname.c_str());
    Biomaout.precision(15);
    ofstream Stockout(fCname.c_str());
    Stockout.precision(15);
    const unsigned int NumLon = _maxX-_minX+1;
    const unsigned int NumLat = _maxY-_minY+1;
    Abundout<<"ncols "<<NumLon<<endl;
    Abundout<<"nrows "<<NumLat<<endl;
    Abundout<<"xllcorner     -180.0"<<endl;
    Abundout<<"yllcorner     -65.0"<<endl;
    Abundout<<"cellsize "<< Parameters::Get()->GetGridCellSize( )<<endl;
    Abundout<<"NODATA_value  -9999"<<endl;
    Biomaout<<"ncols "<<NumLon<<endl;
    Biomaout<<"nrows "<<NumLat<<endl;
    Biomaout<<"xllcorner     -180.0"<<endl;
    Biomaout<<"yllcorner     -65.0"<<endl;
    Biomaout<<"cellsize "<< Parameters::Get()->GetGridCellSize( )<<endl;
    Biomaout<<"NODATA_value  -9999"<<endl;
    Stockout<<"ncols "<<NumLon<<endl;
    Stockout<<"nrows "<<NumLat<<endl;
    Stockout<<"xllcorner     -180.0"<<endl;
    Stockout<<"yllcorner     -65.0"<<endl;
    Stockout<<"cellsize "<< Parameters::Get()->GetGridCellSize( )<<endl;
    Stockout<<"NODATA_value  -9999"<<endl;
    for (int la=NumLat-1;la>=0;la--){
     for (int lo=0;lo<NumLon;lo++){
       int cellIndex=lo-_minX+(_maxX-_minX+1)*(la-_minY);
       Abundout<<_FinalCohortAbundanceMap[cellIndex]<<" ";
       Biomaout<<_FinalCohortBiomassMap[cellIndex]<<" ";
       Stockout<<_FinalStockBiomassMap[cellIndex]<<" ";
     }
     Abundout<<endl;
     Biomaout<<endl;
     Stockout<<endl;
     }
}
//------------------------------------------------------------------------------------------------------------
void MadModel::setupNcOutput(){
        
        std::string filePath = _filePrefix+"totalCohortBreakdown"+_filePostfix+".nc";
        
        netCDF::NcFile cohortBreakdownFile( filePath.c_str(), netCDF::NcFile::replace ); // Creates file
        netCDF::NcDim TimeNcDim = cohortBreakdownFile.addDim( "time", _stopAt ); // Creates dimension
        netCDF::NcVar TimeNcVar = cohortBreakdownFile.addVar( "time", netCDF::ncUint, TimeNcDim ); // Creates variable
        TimeNcVar.putVar( Parameters::Get( )->GetMonthlyTimeStepArray( ) );
        TimeNcVar.putAtt( "units", "month" );
                
        netCDF::NcDim FGroupDim = cohortBreakdownFile.addDim("functionalGroupNumber" , _FinalCohortBreakdown.size() );
        netCDF::NcVar FGNcVar = cohortBreakdownFile.addVar( "functionalGroupNumber", netCDF::ncInt, FGroupDim );
        int* FGnums=new int[_FinalCohortBreakdown.size()];
        for (int i=0;i<_FinalCohortBreakdown.size();i++)FGnums[i]=i+1;
        FGNcVar.putVar( FGnums );
        FGNcVar.putAtt( "units", "number" );
            
        std::vector< netCDF::NcDim > dataDimensions={TimeNcDim,FGroupDim};

        netCDF::NcVar FGNumNcVar = cohortBreakdownFile.addVar(  "numberOfCohortsInFunctionalGroup", netCDF::ncInt, dataDimensions );
        FGNumNcVar.putAtt("units", "number" );

        //***//
        setNcGridFile("totalCohortBiomass", "kg/sq. km.");
        setNcGridFile("totalStockBiomass", "kg/sq. km.");
        setNcGridFile("totalCohortAbundance", "number/sq. km.");


}
//------------------------------------------------------------------------------------------------------------
void MadModel::setNcGridFile(std::string GridName, std::string units){

        std::string filePath = _filePrefix+GridName+_filePostfix+".nc";
        
        netCDF::NcFile gridFile( filePath.c_str(), netCDF::NcFile::replace ); // Creates file

        netCDF::NcDim gTimeNcDim = gridFile.addDim( "time", _stopAt ); // Creates dimension
        netCDF::NcVar gTimeNcVar = gridFile.addVar( "time", netCDF::ncUint, gTimeNcDim ); // Creates variable
        gTimeNcVar.putVar( Parameters::Get( )->GetMonthlyTimeStepArray( ) );
        gTimeNcVar.putAtt( "units", "month" );
                
        netCDF::NcDim longitudeDim =   gridFile.addDim( "Longitude", Parameters::Get( )->GetLengthUserLongitudeArray( ) );
        netCDF::NcVar longitudeNcVar = gridFile.addVar( "Longitude", netCDF::ncFloat, longitudeDim );
        longitudeNcVar.putVar( Parameters::Get( )->GetUserLongitudeArray( ) );
        longitudeNcVar.putAtt( "units", "degrees" );

        netCDF::NcDim latitudeDim =   gridFile.addDim( "Latitude", Parameters::Get( )->GetLengthUserLatitudeArray( ) );
        netCDF::NcVar latitudeNcVar = gridFile.addVar( "Latitude", netCDF::ncFloat, latitudeDim );
        latitudeNcVar.putVar( Parameters::Get( )->GetUserLatitudeArray( ) );
        latitudeNcVar.putAtt( "units", "degrees" );
                
        std::vector< netCDF::NcDim > gridDimensions={gTimeNcDim,latitudeDim,longitudeDim};

        netCDF::NcVar gridVar = gridFile.addVar(  GridName, netCDF::ncDouble, gridDimensions );
        gridVar.putAtt("units", units );
}
//------------------------------------------------------------------------------------------------------------
void MadModel::netcdfOutput( unsigned step ){

        std::string filePath = _filePrefix+"totalCohortBreakdown"+_filePostfix+".nc";

        try {

            netCDF::NcFile cohortBreakdownFile( filePath.c_str(), netCDF::NcFile::write );
            netCDF::NcVar FGNumNcVar=cohortBreakdownFile.getVar( "numberOfCohortsInFunctionalGroup" );

            vector<size_t> pos={step,0};vector<size_t> num={1,_FinalCohortBreakdown.size()};
            FGNumNcVar.putVar(pos, num,_FinalCohortBreakdown.data() );


        } catch( netCDF::exceptions::NcException& e ) {
                e.what( );
                std::cout << "ERROR> Write to \"" << filePath << "\" failed." << std::endl;
        }
        writeNcGridFile(step,_FinalCohortBiomassMap,"totalCohortBiomass");
        writeNcGridFile(step,_FinalStockBiomassMap,"totalStockBiomass");
        writeNcGridFile(step,_FinalCohortAbundanceMap,"totalCohortAbundance");


}
//------------------------------------------------------------------------------------------------------------

void MadModel::writeNcGridFile(unsigned step, vector<double>& GridDoubleVector,std::string GridName){
        
        std::string filePath = _filePrefix+GridName+_filePostfix+".nc";
        try {

            netCDF::NcFile gridFile( filePath.c_str(), netCDF::NcFile::write );
            netCDF::NcVar gridVar=gridFile.getVar( GridName );

            vector<size_t> pos={step,0,0};vector<size_t> num={1,Parameters::Get( )->GetLengthUserLatitudeArray( ),Parameters::Get( )->GetLengthUserLongitudeArray( )};
            gridVar.putVar(pos, num,GridDoubleVector.data() );
                    

        } catch( netCDF::exceptions::NcException& e ) {
                e.what( );
                std::cout << "ERROR> Write to \"" << filePath << "\" failed." << std::endl;
        }
                
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
//---------------------------------------------------------------------------------------------------------------------------
//***------------------------------------------------TESTING Section----------------------------------------------------***//
//---------------------------------------------------------------------------------------------------------------------------
void MadModel::tests(){
    int rank = repast::RepastProcess::instance()->rank();
    randomizer* random=new RandomRepast;
    //random->SetSeed(100); seed is set from model.props file - see main.cpp
    int nranks=_dimX*_dimY;
    
    //MB main has already got the environmental data - this is stored in the background as a DataLayerSet

    //check values from Parameters have got correctly placed in class grid extents
    assert(Parameters::Get()->GetLengthUserLongitudeArray( )==_maxX-_minX+1);
    assert(Parameters::Get()->GetLengthUserLatitudeArray( )==_maxY-_minY+1);
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
    StockDefinitions::Initialise(_props->getProperty("input.DataDirectory")+"/"+_props->getProperty("input.StockDefinitionsFileName"));
    CohortDefinitions::Initialise(_props->getProperty("input.DataDirectory")+"/"+_props->getProperty("input.CohortDefinitionsFileName"));
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
      //agent id number should be increasing
      assert(c->getId().id()==i);
      c->setup(0,1, E,random);
      _context.addAgent(c);
      discreteSpace->moveTo(id, initialLocation);
      c->setLocation(x,y);
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
      c->setLocation(x,y);
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
        vector<int> loc{int(random->GetUniform()*1000.),int(random->GetUniform()*1000.)};
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
    //remove up to 10 agents locally: move other agents around randomly
    for (auto a:agents){
        if (nr<int(random->GetUniform()*10.)){a->_alive=false;_context.removeAgent(a->getId());nr++;}
        if (a->_alive){
            vector<int> loc{int(random->GetUniform()*781.),int(random->GetUniform()*500.-250)};
            space()->moveTo(a,loc);
        }
    }
    //addup to 100 new agents on this thread and then move these at random 
    int nnew=int(random->GetUniform()*100.);
    for (int i=0;i<nnew;i++){
      repast::AgentId id(Cohort::_NextID, rank, _cohortType);
      id.currentRank(rank);
      Cohort* c = new Cohort(id);
      c->setup(0,1, E,random);
      _context.addAgent(c);
      vector<int> loc{int(random->GetUniform()*37.-99),int(random->GetUniform()*188.)};
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
    //---------------------------------------------------
    //***-------------------TEST 11-------------------***//
    //---------------------------------------------------
    //test whether data in cohorts move correctly across threads
    //first remove all existing agents
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    if (rank==0)cout<<"Test11: check agent data moves correctly across threads"<<endl;

    for (auto a:agents){
        _context.removeAgent(a->getId());
    }
    sync();
    //now add a new agent on each thread and set its properties to known values
    for (int i=0;i<n;i++){
      repast::AgentId id(Cohort::_NextID, rank, _cohortType);
      id.currentRank(rank);
      Cohort* c = new Cohort(id);
      c->setup(0,1, E,random);
      _context.addAgent(c);
      discreteSpace->moveTo(id, initialLocation);
      c->setLocation(x,y);
      setupCohortTestValues(c);
      checkCohortTestValues(c);
     }
    //move agents randomly 100 times
    for (int i=0;i<100;i++){
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    for (auto a:agents){
            vector<int> loc{int(random->GetUniform()*237.-125),int(random->GetUniform()*981.-250)};
            space()->moveTo(a,loc);
    }
    sync();}
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    for (auto a:agents){
      checkCohortTestValues((Cohort*) a);
    }
    //change some values and check they are still OK
    for (auto a:agents){
      ((Cohort*)a)->setLocation(112,-98);
      ((Cohort*)a)->_AdultMass=10.;
    }

    for (auto a:agents){
            vector<int> loc{int(random->GetUniform()*237.-125),int(random->GetUniform()*981.-250)};
            space()->moveTo(a,loc);
    }
    sync();
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    for (auto a:agents){
      assert(((Cohort*)a)->_location[0]==112);
      assert(((Cohort*)a)->_location[1]==-98);
      assert(((Cohort*)a)->_AdultMass==10.);
    }
    if (rank==0)cout<<"Test11 succeeded: values copied across threads unchanged"<<endl;

    cout.flush();
    sync();
    //---------------------------------------------------
    //***-------------------TEST 12-------------------***//
    //---------------------------------------------------
    //check out the RNG speed
    typedef std::chrono::high_resolution_clock Clock;
    auto t0 = Clock::now();
    //slightly faster to use a pointer to the singleton instance than re-create it each time
    repast::Random* R=repast::Random::instance();
    repast::DoubleUniformGenerator gen = R->createUniDoubleGenerator(0, 1);
    for (int i=0;i<1000000;i++){

        R->nextDouble();
        //gen.next( );
        repast::NormalGenerator NJ= R->createNormalGenerator(0,1);
        NJ.next();
        repast::LogNormalGenerator LNJ= R->createLogNormalGenerator(1,2.7);
        LNJ.next();
        //Unfortunately using the wrapper class RandomRepast - as in
        //random->GetUniform();
        //random->GetNormal();
        //random->GetLogNormal(1,2.7);
        //currently seems to slow the RNG calls down by more than a factor of two
        //using g++ -O3 in stead of -O2 also seems to have a negative effect (!)
    }


    
    auto t1 = Clock::now();
    //assert(std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count()<160000000);
    if (rank==0)cout<<"Test 12 : elapsed time"<<std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count()<<"ns: RNG alternative in subdirectory randomizers get about 140000000"<<endl;

    //---------------------------------------------------
    //***-------------------TEST 13-------------------***//
    //---------------------------------------------------
    //values from model.props.tests read into Cohort should match those in the parameter file
    Cohort::setParameters(_props);
    if (rank==0)cout<<"Test 13 : check values read into Cohort parmeters match those in model.props.tests"<<endl;

    assert(1== Cohort::_edibleFractionMarine);
    assert(2== Cohort::_AttackRateExponentMarine);
    assert(0.7== Cohort::_HandlingTimeExponentMarine);
    assert(0.7== Cohort::_HandlingTimeScalarMarine);
    assert(0.1== Cohort::_edibleFractionTerrestrial);
    assert(2== Cohort::_AttackRateExponentTerrestrial);
    assert(0.7== Cohort::_HandlingTimeExponentTerrestrial);
    assert(0.7== Cohort::_HandlingTimeScalarTerrestrial);
    assert(1== Cohort::_HerbivoryRateMassExponent);
    assert(1e-11== Cohort::_HerbivoryRateConstant);
    assert(1== Cohort::_ReferenceMass);
    assert(0.5== Cohort::_HandlingTimeScalar_C);
    assert(0.7== Cohort::_HandlingTimeExponent_C);
    assert(1e-06== Cohort::_SearchRateConstant);
    assert(0.7== Cohort::_FeedingPreferenceStandardDeviation);
    assert(12== Cohort::_NumberOfBins);
    assert(0.0278== Cohort::_DispersalSpeedBodyMassScalar);
    assert(0.48== Cohort::_DispersalSpeedBodyMassExponent);
    assert(100== Cohort::_HorizontalDiffusivity);
    assert(18== Cohort::_AdvectiveModelTimeStepLengthHours);
    assert(abs(6.48- Cohort::_HorizontalDiffusivityKmSqPerADTimeStep)<1.e-14);
    assert(40== Cohort::_AdvectionTimeStepsPerModelTimeStep);
    assert(2592== Cohort::_VelocityUnitConversion);
    assert(50000== Cohort::_DensityThresholdScaling);
    assert(0.8== Cohort::_StarvationDispersalBodyMassThreshold);
    assert(6.61== Cohort::_TerrestrialWarmingToleranceIntercept);
    assert(1.6== Cohort::_TerrestrialWarmingToleranceSlope);
    assert(1.51== Cohort::_TerrestrialTSMIntercept);
    assert(1.53== Cohort::_TerrestrialTSMSlope);
    assert(abs(3.14159265358979- Cohort::_Pi)<1.e-14);
    assert(100== Cohort::_CellAreaToHectares);
    assert(0.88== Cohort::_MetabolismMassExponentEcto);
    assert(148984000000== Cohort::_NormalizationConstantEcto);
    assert(0.69== Cohort::_ActivationEnergyEcto);
    assert(41918272883== Cohort::_NormalizationConstantBMR);
    assert(0.69== Cohort::_BasalMetabolismMassExponent);
    assert(0.0366972477064== Cohort::_EnergyScalarEcto);
    assert(0.7== Cohort::_MetabolismMassExponentEndo);
    assert(908090839730== Cohort::_NormalizationConstantEndo);
    assert(0.69== Cohort::_ActivationEnergyEndo);
    assert(8.617e-05== Cohort::_BoltzmannConstant);
    assert(0.0366972477064== Cohort::_EnergyScalarEndo);
    assert(273== Cohort::_TemperatureUnitsConvert);
    assert(1.5== Cohort::_MassRatioThreshold);
    assert(0.95== Cohort::_MassEvolutionProbabilityThreshold);
    assert(0.05== Cohort::_MassEvolutionStandardDeviation);
    assert(0.5== Cohort::_SemelparityAdultMassAllocation);
    assert(0.001== Cohort::_MortalityRateBackground);
    assert(0.003== Cohort::_MortalityRateMature);
    assert(0.6== Cohort::_LogisticInflectionPoint);
    assert(0.05== Cohort::_LogisticScalingParameter);
    assert(1== Cohort::_MaximumStarvationRate);
    if (rank==0)cout<<"Test 13 : succeeded"<<endl;
    //---------------------------------------------------
    //***-------------------TEST 14-------------------***//
    //---------------------------------------------------
    //test whether data in cohorts can be saved out to a boost archive file (intended for model restarts from a saved state -although randoms will make this non-reproducible)
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    if (rank==0)cout<<"Test14: check agent data saves correctly to a file and can be restored"<<endl;

    for (auto a:agents){
        _context.removeAgent(a->getId());
    }
    sync();
    //now add a new agent on each thread and set its properties to known values
    for (int i=0;i<n;i++){
      repast::AgentId id(Cohort::_NextID, rank, _cohortType);
      id.currentRank(rank);
      Cohort* c = new Cohort(id);
      c->setup(0,1, E,random);
      _context.addAgent(c);
      discreteSpace->moveTo(id, initialLocation);
      c->setLocation(x,y);
      setupCohortTestValues(c);
      checkCohortTestValues(c);
     }
     sync();
     agents.clear();
     _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
     std::stringstream s;
     s<<rank;
     std::ofstream ofs("TestAgentSerialization_rank_"+s.str());
     {
      boost::archive::text_oarchive oa(ofs);
      for (auto a:agents){
         AgentPackage package;
         package.setId(a->getId());
         ((Cohort*)a)->PushThingsIntoPackage( package );
         _packages.push_back(package);
      }
      oa<<_packages;
     }
     _packages.clear();
     std::ifstream ifs("TestAgentSerialization_rank_"+s.str());

     {
      boost::archive::text_iarchive ia(ifs);
      ia>>_packages;
      for (auto& p:_packages){
         repast::AgentId id=p.getId();
         int n=Cohort::_NextID;
         //Cohort::_NextID should be incremented here so set increase flag to true
         MadAgent* a=new Cohort( id,p,true );
         assert(Cohort::_NextID==n+1);
         assert(a->getId()==p._id);
         assert(a->getId().currentRank()==rank);
         assert(a->getId().agentType()==_cohortType);
         checkCohortTestValues((Cohort*)a);
      }
     }
     if (rank==0)cout<<"Test 14 : succeeded"<<endl;

    //---------------------------------------------------
    //***-------------------TEST 15-------------------***//
    //---------------------------------------------------
    //test whether dispersal using fractional cell co-ordinates works.
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    if (rank==0)cout<<"Test15: check direct agent dispersal works (movement by fraction of a cell permitted)"<<endl;

    for (auto a:agents){
        _context.removeAgent(a->getId());
    }
    sync();
    _noLongitudeWrap=true;
    //now add a new agent on each thread and set its properties to known values
    for (int i=0;i<n;i++){
      repast::AgentId id(Cohort::_NextID, rank, _cohortType);
      id.currentRank(rank);
      Cohort* c = new Cohort(id);
      E=_Env[x - _minX + (_maxX - _minX + 1)*(y - _minY)];
      c->setup(0,1, E,random);
      _context.addAgent(c);
      discreteSpace->moveTo(id, initialLocation);
      c->setLocation(x,y);c->setDestination(x,y);
      double uSpeed=0.9*E->Width(),vSpeed=0.5*E->Height();
      c->TryToDisperse(uSpeed, vSpeed,E, this);
      if (c->_moved){
        assert(c->_destination[0]-x-0.9<1.e-14);
        assert(c->_destination[1]-y-0.5<1.e-14);
        vector<int>newPlace={0,0};
        newPlace[0]=int(c->_destination[0]);
        newPlace[1]=int(c->_destination[1]);
        space()->moveTo(c,newPlace);
        c->_location=c->_destination;
      }
          
     }
     sync();
     agents.clear();
     _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
     //nothing should have moved rank at this point
     for (auto a:agents){
         assert(a->getId().currentRank()==a->getId().startingRank());
         assert(((Cohort*)a)->_location[0]-x-0.9<1.e-14);
         assert(((Cohort*)a)->_location[1]-y-0.5<1.e-14);
         vector<int>newPlace={_minX,(_maxY + _minY)/2};
         space()->moveTo(a,newPlace);
         ((Cohort*)a)->setDestination(newPlace[0],newPlace[1]);
         ((Cohort*)a)->setLocation(newPlace[0],newPlace[1]);
     }
     sync();
     agents.clear();
     _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
     //all agents now on same thread
     _noLongitudeWrap=false;
     _dispersalSelection="direct";
     double mean=0;int c=0;
     for (auto a:agents){
       c++;
       double uSpeed=0.02*E->Width(),vSpeed=0;
       ((Cohort*)a)->_Realm="all";//allow dispersal on land or sea
       for (int i=0; i<10000;i++){
        double r=repast::Random::instance()->nextDouble();
        double k;
        if (r>=0.5)k=-1; else k=1;
        ((Cohort*)a)->TryToDisperse(k*uSpeed, vSpeed,E, this);
     }
     //agents dispersing at random from a point should end up gaussian distributed with
     //mean 0 (around the initial point) and stdev stepsize^2*number of steps
     mean=mean+((Cohort*)a)->_destination[0];
     if (((Cohort*)a)->_destination[0] < (_maxX + _minX)/2) mean=mean+(_maxX - _minX);
    }
    if (agents.size()!=0){
      mean=mean/c;double stdev=0;
      for (auto a:agents){
          double u=((Cohort*)a)->_destination[0];
          if (u < (_maxX + _minX)/2) u+=(_maxX - _minX);//account for grid wrap - agents start from _minX
          stdev=stdev+pow((u-mean),2);
      }
      stdev=sqrt(stdev/c);
      cout<<"mean and stddev of position after dispersal "<<mean<<" "<<stdev<<endl;
      assert(mean-(_maxX - _minX)<0.02);
      assert(stdev-2<0.01);
      cout<<"Test 15 : succeeded"<<endl;
    }
}    
//---------------------------------------------------------------------------------------------------------------------------
//define some data values for the Cohort to check whether they are preserved on moving across threads
//see Cohort::setup
void MadModel::setupCohortTestValues(Cohort* c){
    

    c->_FunctionalGroupIndex=12;//functionalGroup;
    c->_Merged                      = false;
    c->_alive                       = true;

	c->_Heterotroph=false;//(CohortDefinitions::Get()->Trait(functionalGroup     , "heterotroph/autotroph")  =="heterotroph");   
    c->_Autotroph  =!c->_Heterotroph;
    c->_Endotherm  =true;//(CohortDefinitions::Get()->Trait(functionalGroup     , "endo/ectotherm")         =="endotherm");
    c->_Ectotherm  =!c->_Endotherm;
    c->_Realm      ="Astring";//CohortDefinitions::Get()->Trait(functionalGroup     , "realm");

    c->_Iteroparous=true;//(CohortDefinitions::Get()->Trait(functionalGroup     , "reproductive strategy")  =="iteroparity");
    c->_Semelparous=!c->_Iteroparous;
    c->_Herbivore=false;//(CohortDefinitions::Get()->Trait(functionalGroup       , "nutrition source")       =="herbivore");
    c->_Carnivore=true;//(CohortDefinitions::Get()->Trait(functionalGroup       , "nutrition source")       =="carnivore");
    c->_Omnivore= false;//(CohortDefinitions::Get()->Trait(functionalGroup       , "nutrition source")       =="omnivore");
    c->_IsPlanktonic= false;//(CohortDefinitions::Get()->Trait(functionalGroup   , "mobility")               =="planktonic");
    c->_IsFilterFeeder=false;//(CohortDefinitions::Get()->Trait(functionalGroup   , "diet")                   =="allspecial");
    
    c->_ProportionSuitableTimeActive= 0.67;//CohortDefinitions::Get()->Property(functionalGroup   ,"proportion suitable time active");
    
    c->_IsMature=false;
    c->_IndividualReproductivePotentialMass=0;
    
    c->_AssimilationEfficiency_H=0.15;//CohortDefinitions::Get()->Property(functionalGroup   ,"herbivory assimilation");
    c->_AssimilationEfficiency_C=0.321;//CohortDefinitions::Get()->Property(functionalGroup   ,"carnivory assimilation");
    c->_BirthTimeStep=0;
    c->_MaturityTimeStep=std::numeric_limits<unsigned>::max( );
    c->_MinimumMass=0.01;//CohortDefinitions::Get()->Property(functionalGroup   ,"minimum mass");
    c->_MaximumMass=75000;//CohortDefinitions::Get()->Property(functionalGroup   ,"maximum mass");

    //repast::DoubleUniformGenerator gen = repast::Random::instance()->createUniDoubleGenerator(0, 1);
    
    c->_AdultMass = 48.77;//pow( 10, ( gen.next( ) * ( log10( _MaximumMass ) - log10( 50 * _MinimumMass ) ) + log10( 50 * _MinimumMass ) ) );
    

    //NormalGenerator NJ= repast::Random::instance()->createNormalGenerator(0.1,0.02);
    
    c->_LogOptimalPreyBodySizeRatio = 0.05;//log(std::max( 0.01, NJ.next() ));
    
    //double expectedLnAdultMassRatio = 2.24 + 0.13 * log( _AdultMass );
    //in the original code the mean and sd are those of the underlying normal distribution
    //in the boost library they refer to the log distibution - see
    //https://www.boost.org/doc/libs/1_43_0/libs/math/doc/sf_and_dist/html/math_toolkit/dist/dist_ref/dists/lognormal_dist.html
    //LogNormalGenerator LNJ= repast::Random::instance()->createLogNormalGenerator(exp(expectedLnAdultMassRatio+0.5*0.5/2.), (exp(0.5*0.5)-1)*exp(2*expectedLnAdultMassRatio+0.5*0.5));
    
    /*if( _Realm=="terrestrial" ) {
          do {
            _JuvenileMass = _AdultMass  / (1.0 + LNJ.next());
          } while( _AdultMass <= _JuvenileMass || _JuvenileMass < _MinimumMass );
    } else {
          do {
            _JuvenileMass = _AdultMass  / (1.0 + 10 *LNJ.next());
          } while( _AdultMass <= _JuvenileMass || _JuvenileMass < _MinimumMass );
    }

    double NewBiomass = ( 3300. / numCohortsThisCell ) * 100 * 3000 * pow( 0.6, ( log10( _JuvenileMass ) ) ) * ( e->Area() );
    */
    
    c->_CohortAbundance = 100000;//NewBiomass / _JuvenileMass;
    c->_JuvenileMass=12.5;
    
    c->_MaximumAchievedBodyMass=99;//_JuvenileMass;
    c->_IndividualBodyMass=12567;//_JuvenileMass;
    c->_moved=false;
    c->_location={-12,75};
}
//---------------------------------------------------------------------------------------------------------------------------
//check that the values set in the above function are still maintained
void MadModel::checkCohortTestValues(Cohort* c){
    assert(c->_FunctionalGroupIndex==12);
    assert(c->_Merged                      == false);
    assert(c->_alive                       == true);

	assert(c->_Heterotroph==false);
    assert(c->_Autotroph  ==!c->_Heterotroph);
    assert(c->_Endotherm  ==true);
    assert(c->_Ectotherm  ==!c->_Endotherm);
    assert(c->_Realm      =="Astring");

    assert(c->_Iteroparous==true);
    assert(c->_Semelparous==!c->_Iteroparous);
    assert(c->_Herbivore==false);
    assert(c->_Carnivore==true);
    assert(c->_Omnivore== false);
    assert(c->_IsPlanktonic== false);
    assert(c->_IsFilterFeeder==false);
    
    assert(c->_ProportionSuitableTimeActive== 0.67);
    
    assert(c->_IsMature==false);
    assert(c->_IndividualReproductivePotentialMass==0);
    
    assert(c->_AssimilationEfficiency_H==0.15);
    assert(c->_AssimilationEfficiency_C==0.321);
    assert(c->_BirthTimeStep==0);
    assert(c->_MaturityTimeStep==std::numeric_limits<unsigned>::max( ));
    assert(c->_MinimumMass==0.01);
    assert(c->_MaximumMass==75000);

    assert(c->_AdultMass == 48.77);
    
    assert(c->_LogOptimalPreyBodySizeRatio == 0.05);

    
    assert(c->_CohortAbundance == 100000);
    assert(c->_JuvenileMass==12.5);
    assert(c->_MaximumAchievedBodyMass==99);
    assert(c->_IndividualBodyMass==12567);
    assert(c->_moved==false);
    assert(c->_location[0]==-12);
    assert(c->_location[1]==75);
}
//---------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------
//***------------------------------------------------END TESTING Section----------------------------------------------------***//
//---------------------------------------------------------------------------------------------------------------------------
