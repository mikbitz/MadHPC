/* Model.cpp */

#include <stdio.h>
#include <vector>
#include <sstream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/archive_exception.hpp>
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
#include "EnvironmentCell.h"
#include "Parameters.h"
#include "Groups.h"
#include "Stock.h"
#include "Cohort.h"
#include "Human.h"
#include "FileReader.h"
#include "CohortSum.h"
#include "TimeStep.h"
#include "Constants.h"
#include "randomizer.h"
#include "RandomRepast.h"
#include "AgentPackage.h"
#include "UtilityFunctions.h"

#include <netcdf>

//repast shuffleList only works on pointers
//this version is for vectors of int

void shuffleList(std::vector<int>& elementList){
  if(elementList.size() <= 1) return;
  repast::IntUniformGenerator rnd = repast::Random::instance()->createUniIntGenerator(0, elementList.size() - 1);
  int swap;
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
    //switch on all model aspects - these might need to be switched off for test purposes.
    _eating       = props.getProperty("simulation.IncludeEating")      !="false";
    _metabolism   = props.getProperty("simulation.IncludeMetabolism")  !="false";
    _reproduction = props.getProperty("simulation.IncludeReproduction")!="false";
    _death        = props.getProperty("simulation.IncludeDeath")       !="false";
    _dispersal    = props.getProperty("simulation.IncludeDispersal")   !="false";
    _mergers      = props.getProperty("simulation.IncludeMergers")     !="false";
    _output       = props.getProperty("simulation.IncludeOutput")      !="false";
    _verbose      = props.getProperty("verbose")=="true";
    //-----------------
    //Pull in parameter from the model.props file
    _props = &props;
    //Number of timesteps
	_stopAt = repast::strToInt(_props->getProperty("stop.at"));
    
    std::string rstrt=props.getProperty("simulation.RestartEvery");
    if (rstrt!="")_restartInterval=repast::strToInt(rstrt); else _restartInterval=0;

    rstrt=props.getProperty("simulation.RestartStep");
    if (rstrt!="")_restartStep=repast::strToInt(rstrt); else _restartStep=0;
    
    rstrt=props.getProperty("simulation.RestartDirectory");
    if (rstrt!="")_restartDirectory=rstrt+"/"; else _restartDirectory="./";

    _archiveFormat ="binary";
    if (props.getProperty("simulation.RestartFormat")=="text")_archiveFormat="text";

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
    if(repast::RepastProcess::instance()->rank() == 0 && _output){
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
    //slightly tricky to get the file prefix to all other threads (needed, for example, for restart names)
    //only thread 0 can know the name since it needs to create a new name based on existing directory names on disk 
    int prefix_size = _filePrefix.size();
    MPI_Bcast(&prefix_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (repast::RepastProcess::instance()->rank() != 0)_filePrefix.resize(prefix_size);
    MPI_Bcast(const_cast<char*>(_filePrefix.data()), prefix_size, MPI_CHAR, 0, MPI_COMM_WORLD);
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
void MadModel::initSchedule(unsigned startingStep,repast::ScheduleRunner& runner){
	runner.scheduleEvent(startingStep, 1, repast::Schedule::FunctorPtr(new repast::MethodFunctor<MadModel> (this, &MadModel::step)));
	runner.scheduleStop(_stopAt);
    runner.scheduleEndEvent(Schedule::FunctorPtr(new MethodFunctor<MadModel> (this, &MadModel::dataSetClose)));

}
//------------------------------------------------------------------------------------------------------------
void MadModel::init(unsigned startingStep){
    _startingStep=startingStep;//starts at 1 to correspond to timestep 0 - confusing...needs a fix
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
    
    //get the environmental data - this is stored in the background as a DataLayerSet - NB this is currently done in main.cpp because of some strange problems with legacy fileReading code
    //FileReader F;
    //F.ReadFiles();
    
    //now set up the environmental cells - note at present this uses the full grid, not just local to this thread
    //so that off-thread environment can be easily queried. Currently some duplication here, but it is not a huge amount of data.
    _Env=Environment(_minX,_maxX,_minY,_maxY);

    //set up the static (i.e. shared) parameters for the Cohorts
    Cohort::setParameters(_props);
    //get the definitions of stocks and cohorts
    StockDefinitions::Initialise(_props->getProperty("input.DataDirectory")+"/"+_props->getProperty("input.StockDefinitionsFileName"));
    CohortDefinitions::Initialise(_props->getProperty("input.DataDirectory")+"/"+_props->getProperty("input.CohortDefinitionsFileName"));
  
    unsigned numCohortGroups=CohortDefinitions::Get()->size();
    unsigned humanCount = strToInt(_props->getProperty("human.count"));
    unsigned cohortCount = strToInt(_props->getProperty("cohort.count"));
    unsigned numStockGroups = StockDefinitions::Get()->size();

    outputNames.push_back("totalCohortBiomass");
    outputNames.push_back("totalStockBiomass");
    outputNames.push_back("totalCohortAbundance");
    //Values only reduced on thread 0
    if(repast::RepastProcess::instance()->rank() == 0){
      outputUnits["totalCohortBiomass"]   ="kg/sq. km.";
      outputUnits["totalStockBiomass"]    ="kg/sq. km.";
      outputUnits["totalCohortAbundance"] ="number/sq. km.";
      for (auto name: outputNames) outputMaps[name]    =  vector<double> ( (_maxX-_minX+1) * (_maxY-_minY+1),0.0 );

     _FinalCohortBreakdown.resize(numCohortGroups);
    }
    
    
    randomizer* random=new RandomRepast;
    //seed is set from model.props file - see main.cpp

 
    //explicitly use the local bounds of the grid on this thread to create countOfAgents per cell.
    //Not doing this can lead to problems with agents in distant cells not within the local thread neighbourhood
    //see SharedBaseGrid.h moveTo method
    //although latest updates I have made to RHPC should have fixed this...
    unsigned totalCohorts=0,totalStocks=0;

    unsigned cNum=0,sNum=0;

    unsigned totalStocksThisCell=0;

  
    int s=0;
    unsigned hF=10000;//functionalgroup ID for humans
    if (_restartStep==0){
    for (int x = _xlo; x < _xhi; x++){
        for (int y = _ylo; y < _yhi; y++){
             EnvironmentCell* E=_Env[x][y];
             repast::Point<int> initialLocation(x,y);


             if (E->_Realm=="terrestrial"){
               for (unsigned j=0;j<humanCount;j++){
                    //make sure the agentId is unique on this thread!!
                    // values are int id, int startProc, int agentType, 
                    repast::AgentId id(Human::_NextID, rank, _humanType);
                    //agent also needs id of its current thread
                    id.currentRank(rank);
                    Human* h = new Human(id);
                    h->setup(hF,humanCount, E,random);
                    _context.addAgent(h);
                    discreteSpace->moveTo(id, initialLocation);
                    //to get movement right agent needs its own copy of location
                    h->setLocation(x,y);
               }
             }
             
             unsigned totalCohortsThisCell=0;
             for (unsigned i=0;i<numCohortGroups;i++) if (E->_Realm==CohortDefinitions::Get()->Trait(i,"realm"))totalCohortsThisCell+=cohortCount;
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
                         double xr=0,yr=0;
                         if (_dispersalSelection=="direct"){
                         //allow cohorst to be at locations other than cell centres initially - note using fractional cell co-ordinates
                             xr=(1-2*random->GetUniform())*0.5;
                             yr=(1-2*random->GetUniform())*0.5;
                         }
                         c->setLocation(x+xr,y+yr);
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

        }
        if (_restartStep>0){
            read_restart(_restartStep);
            vector<MadAgent*>agents;
            _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
            for (auto a:agents){
                if (a->getId().agentType()==_cohortType){
                    Cohort* c=(Cohort *)a;
                    _totalCohorts++;
                    _totalCohortAbundance += c->_CohortAbundance;
                    _totalCohortBiomass += ( c->_IndividualBodyMass + c->_IndividualReproductivePotentialMass ) * c->_CohortAbundance / 1000.;//g to kg
                }
                if (a->getId().agentType()==_stockType){
                    Stock* s=(Stock* )a;
                    _totalStockBiomass+=s->_TotalBiomass/1000; //g to kg
                }
            }
        }
    if(_verbose)cout<<"rank "<<rank<<" totalCohorts "<<_totalCohorts<<" totalStocks "<<s<<endl;
    if (_output){
     setupOutputs();
     if (rank==0)setupNcOutput();
    }
    long double t = initTimer.stop();
	std::stringstream ss;
	ss << t;
	_props->putProperty("init.time", ss.str());
    sync();

}
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
    //needed to advance the environmental datalayers to the current timestep
   
    _Env.update(CurrentTimeStep);
    
	if(rank == 0) std::cout << " TICK " << CurrentTimeStep << std::endl;

    //vectors length CohortDefinitions::Get()->size() initialized to 0
    vector<int> cohortBreakdown(CohortDefinitions::Get()->size(),0);
    //spatial distributions - for efficiency these should just the the local part, but easier initially to just keep the whole thing on each thread
    //with zeros for non-local, and then reduce onto thread 0
    map< string,vector<double> > localMaps;
    for (auto name: outputNames) localMaps[name]    =  vector<double> ( (_maxX-_minX+1) * (_maxY-_minY+1),0.0 );

    int buffer=repast::strToInt(_props->getProperty("grid.buffer"));
    
    //loop over local grid cells
    //For cell interaction, we need to make sure copies are updated appropriately.
    //Since we use a moore neighbourhood, blocks of nine cells can be treated independently when the centre cell of the nine is updated
    //- so split the grid into 9 sets of 3x3 cells, 
    //separate by a 3 cell stride, update date these, then sync to update the shared data in overlaps and then move to the next block. 
    //Randomize the starting point on each time step to ensure certain cells do not get privileged.
    
    vector<int> initialCells(9,0);
    if (rank==0){
    for(int i=0;i<9;i++){
        initialCells[i]=i%3;
        //use cell coords 0,0 0,1 0,2; 1,0 1,1 1,2 etc.
        if (i>2)initialCells[i]+=(_maxX-_minX+1);
        if (i>5)initialCells[i]+=(_maxX-_minX+1);
    }
    
    shuffleList(initialCells);//use the version defined above
    }
    //cells need to have the same random sequence on each thread (otherwise cross-thread interactions won't preserve sequences), 
    //Do not try to do this by  restarting the rng: random sequences will no longer be random
    //Instead permute on thread 0 then broadcast to make sure all threads use the same cell permutation
    MPI_Bcast(initialCells.data(), 9, MPI_INT, 0, MPI_COMM_WORLD);


    int range=0;
    //cohorts need an interaction range in the case of interaction
    //this should be strictly less than one cell at present
    if (buffer==2) range=1;
    //create a random sequence over the cohorts in each cell. make sure this is shared across threads so that overlap regions get the same data.
    for(int y = _ylo; y < _yhi; y++){  
     for(int x = _xlo; x < _xhi; x++){
        repast::Point<int> location(x,y);
        std::vector<MadAgent*> agentsInCell;
        //query four neighbouring cells, distance 0 (i.e. just the centre cell) - "true" keeps the centre cell.
        //these will be the cohorts that act in the current cell
        repast::VN2DGridQuery<MadAgent> VN2DQuery(space());
        VN2DQuery.query(location, 0, true, agentsInCell);
        vector<Cohort*> cohorts;
        for (auto a:agentsInCell){
            if ( a->_alive && a->getId().currentRank()==rank){//agents must be living and local!
            //separate out stocks (currently all plants) from cohorts(animals)
               if (a->getId().agentType()==_cohortType)cohorts.push_back( (Cohort*) a);
            }
        }

        vector<int> seq(cohorts.size(),0);
        for (int i=0;i<cohorts.size();i++)seq[i]=i;
        shuffleList(seq);
        for (int i=0;i<cohorts.size();i++)cohorts[i]->_sequencer=seq[i];
     }
    }
    //synchronize the sequence data
    if (buffer==2)repast::RepastProcess::instance()->synchronizeAgentStates<AgentPackage, MadAgentPackageProvider, MadAgentPackageReceiver>(*provider, *receiver);
    
    //now dynamical loop over cells, in 9 chunks
    for (unsigned i=0;i<9;i++){
        int initialCell=initialCells[i];
        int x0=initialCell%(_maxX-_minX+1)+_minX; 
        int y0=initialCell/(_maxX-_minX+1)+_minY;
        while (x0<_xlo)x0+=3;
        while (y0<_ylo)y0+=3;
        for(int y = y0-range; y < _yhi+range; y+=3){  
         for(int xt = x0-range; xt < _xhi+range; xt+=3){
            //pick out the cells that are on this thread plus those in the overlap region (but only those in range need to actually do the dynamics)
            //overlap cells have to run updates in the same order on each thread so that dynamics is reproduced on independent threads
            //NB poles are currently a potential problem as Repast thinks we are on a torus! Must exclude overlaps at max and min latitudes.
            //In practice the internal locations of Cohorts should make this happen as distance measures will be large (x distances get wrapped, but not y, see Cohort.cpp). Stocks need locations though.
            int x=xt;

            if (!_noLongitudeWrap){
               if (xt  < _minX){x = x + (_maxX - _minX + 1);}
               if (xt  > _maxX){x = x - (_maxX - _minX + 1);}
            }

            if (y >= _minY && y <= _maxY){
             if (x >=_minX && x <= _maxX) {

                EnvironmentCell* E=_Env[x][y];
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
                //sort the lists of cohorts using their unique sequence number -easiest to do this in place with a lambda
                //this ensures that random cohort order computed above can be preserved across threads
                sort(cohorts.begin(), cohorts.end(), [](const Cohort* c1, const Cohort* c2) -> bool { return c1->_sequencer > c2->_sequencer;});

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
        //a full sync is not really needed here - just synchronizeAgentStates.

        if (buffer==2)repast::RepastProcess::instance()->synchronizeAgentStates<AgentPackage, 
               MadAgentPackageProvider, MadAgentPackageReceiver>(*provider, *receiver);
    }

    //updates of offspring and mergers/death happen after all cells have updated
    //need to keep this separate if there is cross-cell interaction
 
    for(int y = _ylo; y < _yhi; y++){  
     for(int x = _xlo; x < _xhi; x++){
      int cellIndex=x-_minX+(_maxX-_minX+1)*(y-_minY);
            
            EnvironmentCell* E=_Env[x][y];

            //store current location in a repast structure for later use
            repast::Point<int> location(x,y);
            
            std::vector<MadAgent*> agentsInCell,argents;
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
            for (auto c:cohorts)if(c->_newH!=NULL && _reproduction){
                _context.addAgent(c->_newH);
                discreteSpace->moveTo(c->_newH->getId(), location);
                newCohorts.push_back(c->_newH);
                _totalReproductions++;
            }
            //make sure new cohorts can get marked for death/merged - new cohorts are guaranteed to be local
            for (auto n:newCohorts){cohorts.push_back(n);}
            newCohorts.clear();
            for (auto c:cohorts){c->markForDeath();if (!c->_alive)_totalDeaths++;}

            
            if (_mergers) _totalMerged+=CohortMerger::MergeToReachThresholdFast(cohorts);
            
            //set up output stock map
            for (auto s:stocks){
                _totalStockBiomass+=s->_TotalBiomass/1000;
                localMaps["totalStockBiomass"][cellIndex]+=s->_TotalBiomass/1000/E->Area();
                _totalStocks++;
            }
            //acumulate other totals and spatial maps

            for (auto c:cohorts){
                if (c->_alive){
                    _totalCohorts++;
                    localMaps["totalCohortAbundance"][cellIndex]+= c->_CohortAbundance/E->Area();
                    _totalCohortAbundance += c->_CohortAbundance;
                    localMaps["totalCohortBiomass"][cellIndex]+=( c->_IndividualBodyMass + c->_IndividualReproductivePotentialMass ) * c->_CohortAbundance / 1000.;//convert to kg
                    _totalCohortBiomass += ( c->_IndividualBodyMass + c->_IndividualReproductivePotentialMass ) * c->_CohortAbundance / 1000.;
                    cohortBreakdown[c->_FunctionalGroupIndex]++;
                }
            }

            localMaps["totalCohortBiomass"][cellIndex]=localMaps["totalCohortBiomass"][cellIndex]/E->Area();//per square kilometre

            //care with sync() here - need to get rid of not-alive agents:currently this is a lazy delete for new/non-local agents (they get removed one timestep late)?
            for (auto a:agentsInCell )if (!a->_alive && _death)_context.removeAgent(a->getId());//does this delete the agent storage? - yes if Boost:shared_ptr works OK
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
            EnvironmentCell* E=_Env[x][y];
            repast::Point<int> location(x,y);
            std::vector<MadAgent*> agentsInCell;
            //query four neighbouring cells, distance 0 (i.e. just the centre cell) - "true" keeps the centre cell.
            repast::VN2DGridQuery<MadAgent> VN2DQuery(space());
            VN2DQuery.query(location, 0, true, agentsInCell);

            for (auto a:agentsInCell){
                 //dispersers must be local and alive!
                 if (a->getId().currentRank()==repast::RepastProcess::instance()->rank() && a->_alive){
                     if (a->getId().agentType()==MadModel::_cohortType && _dispersal) {
                         ((Cohort*) a)->moveIt(E,this);
                         if (a->_moved){ //cohorts that have changed cell
                             movers.push_back((Cohort *) a);a->_moved=false;
                        }
                     }
                }
            }
        }
    }
    _totalMoved=movers.size();

    //agent data may have changed locally - ensure this is synced before anything gets moved, otherwise values do not move across threads correctly when there are buffers.
    if (buffer==2)repast::RepastProcess::instance()->synchronizeAgentStates<AgentPackage, 
            MadAgentPackageProvider, MadAgentPackageReceiver>(*provider, *receiver);

    vector<int>newPlace={0,0};
    for (auto& m:movers){
        newPlace[0]=int(m->_destination[0]);
        newPlace[1]=int(m->_destination[1]);
        //move things - these are then settled
        space()->moveTo(m,newPlace);
    }
    // ***** state of the model will not be fully consistent until sync() *****
    sync();
    
    if (_output){
     //map outputs/vectors not yet available via svbuilder or netcdf builder.  
     //numbers per functional group - here's how to add up across a vector over threads.
     MPI_Reduce(cohortBreakdown.data(), _FinalCohortBreakdown.data(), CohortDefinitions::Get()->size(), MPI::INT, MPI::SUM, 0, MPI_COMM_WORLD);

     //also get the maps
     for (auto name:outputNames)MPI_Reduce(localMaps[name].data(), outputMaps[name].data(), (_maxX-_minX+1) * (_maxY-_minY+1), MPI::DOUBLE, MPI::SUM, 0, MPI_COMM_WORLD);

     //if(repast::RepastProcess::instance()->rank() == 0){asciiOutput(CurrentTimeStep);}
     if(repast::RepastProcess::instance()->rank() == 0){netcdfOutput( CurrentTimeStep - _startingStep + 1);}
    }

    if (_restartInterval>0 && CurrentTimeStep>_restartStep && (CurrentTimeStep+1-_restartStep)%_restartInterval==0)write_restart();

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
    Abundout<<"cellsize "<< Parameters::instance()->GetGridCellSize( )<<endl;
    Abundout<<"NODATA_value  -9999"<<endl;
    Biomaout<<"ncols "<<NumLon<<endl;
    Biomaout<<"nrows "<<NumLat<<endl;
    Biomaout<<"xllcorner     -180.0"<<endl;
    Biomaout<<"yllcorner     -65.0"<<endl;
    Biomaout<<"cellsize "<< Parameters::instance()->GetGridCellSize( )<<endl;
    Biomaout<<"NODATA_value  -9999"<<endl;
    Stockout<<"ncols "<<NumLon<<endl;
    Stockout<<"nrows "<<NumLat<<endl;
    Stockout<<"xllcorner     -180.0"<<endl;
    Stockout<<"yllcorner     -65.0"<<endl;
    Stockout<<"cellsize "<< Parameters::instance()->GetGridCellSize( )<<endl;
    Stockout<<"NODATA_value  -9999"<<endl;
    for (int la=NumLat-1;la>=0;la--){
     for (int lo=0;lo<NumLon;lo++){
       int cellIndex=lo-_minX+(_maxX-_minX+1)*(la-_minY);
       //Abundout<<_FinalCohortAbundanceMap[cellIndex]<<" ";
       //Biomaout<<_FinalCohortBiomassMap[cellIndex]<<" ";
       //Stockout<<_FinalStockBiomassMap[cellIndex]<<" ";
     }
     Abundout<<endl;
     Biomaout<<endl;
     Stockout<<endl;
     }
}
//------------------------------------------------------------------------------------------------------------
void MadModel::setupNcOutput(){
        
        std::string filePath = _filePrefix+"totalCohortBreakdown"+_filePostfix+".nc";
        auto times=TimeStep::instance()->TimeStepArray();
        netCDF::NcFile cohortBreakdownFile( filePath.c_str(), netCDF::NcFile::replace ); // Creates file
        netCDF::NcDim TimeNcDim = cohortBreakdownFile.addDim( "time", times.size() ); // Creates dimension
        netCDF::NcVar TimeNcVar = cohortBreakdownFile.addVar( "time", netCDF::ncUint, TimeNcDim ); // Creates variable
        TimeNcVar.putVar( times.data() );
        TimeNcVar.putAtt( "units", TimeStep::instance()->TimeStepUnits() );
                
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
        for (auto name:outputNames)setNcGridFile(name,outputUnits[name]);

}
//------------------------------------------------------------------------------------------------------------
void MadModel::setNcGridFile(std::string GridName, std::string units){

        std::string filePath = _filePrefix+GridName+_filePostfix+".nc";
        
        netCDF::NcFile gridFile( filePath.c_str(), netCDF::NcFile::replace );             // Creates file
        auto times=TimeStep::instance()->TimeStepArray();
        netCDF::NcDim gTimeNcDim = gridFile.addDim( "time", times.size() );                    // Creates dimension
        netCDF::NcVar gTimeNcVar = gridFile.addVar( "time", netCDF::ncUint, gTimeNcDim ); // Creates variable
        gTimeNcVar.putVar( times.data() );
        gTimeNcVar.putAtt( "units", TimeStep::instance()->TimeStepUnits() );
                
        netCDF::NcDim longitudeDim =   gridFile.addDim( "Longitude", Parameters::instance()->GetLengthLongitudeArray( ) );
        netCDF::NcVar longitudeNcVar = gridFile.addVar( "Longitude", netCDF::ncFloat, longitudeDim );
        longitudeNcVar.putVar( Parameters::instance()->GetLongitudeArray( ) );
        longitudeNcVar.putAtt( "units", "degrees" );

        netCDF::NcDim latitudeDim =   gridFile.addDim( "Latitude", Parameters::instance()->GetLengthLatitudeArray( ) );
        netCDF::NcVar latitudeNcVar = gridFile.addVar( "Latitude", netCDF::ncFloat, latitudeDim );
        latitudeNcVar.putVar( Parameters::instance()->GetLatitudeArray( ) );
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
        for (auto name:outputNames)writeNcGridFile(step,outputMaps[name],name);

}
//------------------------------------------------------------------------------------------------------------

void MadModel::writeNcGridFile(unsigned step, vector<double>& GridDoubleVector,std::string GridName){
        
        std::string filePath = _filePrefix+GridName+_filePostfix+".nc";
        try {

            netCDF::NcFile gridFile( filePath.c_str(), netCDF::NcFile::write );
            netCDF::NcVar gridVar=gridFile.getVar( GridName );

            vector<size_t> pos={step,0,0};vector<size_t> num={1,Parameters::instance()->GetLengthLatitudeArray( ),Parameters::instance()->GetLengthLongitudeArray( )};
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
	runner.scheduleEvent(0.1+_startingStep-1, 1, Schedule::FunctorPtr(new MethodFunctor<repast::DataSet> (dataSet,&repast::DataSet::record)));
	Schedule::FunctorPtr dsWrite = Schedule::FunctorPtr(new MethodFunctor<repast::DataSet> (dataSet,&repast::DataSet::write));
    //output every 100 steps
	runner.scheduleEvent(100.2, 100, dsWrite);
}
//---------------------------------------------------------------------------------------------------------------------------
// Restarts
//------------------------------------------------------------------------------------------------------------
 void MadModel::write_restart(){
     //each thread writes its own restart file. Read restart is able to deal with this later...
     unsigned step=RepastProcess :: instance ()->getScheduleRunner ().currentTick ();
     std::vector<MadAgent*> agents;
     _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
     std::stringstream s;
     s<<step<<"_"<<repast::RepastProcess::instance()->rank();
     std::ofstream ofs(_filePrefix+"Restart_step_rank_"+s.str());
     //archive saves when destructor called - this block should ensure this happens
     {

      for (auto a:agents){
          if (a->_alive){

              AgentPackage package(a->getId());
              if (a->getId().agentType()==_cohortType){
                  ((Cohort*)a)->PushThingsIntoPackage( package );
              }
              if (a->getId().agentType()==_stockType){
                  ((Stock*)a)->PushThingsIntoPackage( package );
              }
              _packages.push_back(package);
          }
      }
      
      if (_archiveFormat =="binary") {boost::archive::binary_oarchive oa(ofs);oa<<_packages;}
      if (_archiveFormat =="text"  ) {boost::archive::text_oarchive   oa(ofs);oa<<_packages;}
     
      
      if (_verbose) cout<<"Wrote "<<_packages.size()<<" objects to restart: "<<"Restart_step_rank_"<<s.str()<<endl;
     }
     _packages.clear();

     
 }
//------------------------------------------------------------------------------------------------------------
void MadModel::read_restart(unsigned step){
    
    unsigned r=0,error=0,rank=repast::RepastProcess::instance()->rank();
    unsigned numProcs=repast::RepastProcess::instance()->worldSize();

    std::stringstream s;
    s<<step<<"_"<<r;
    std::string filename=_restartDirectory+"Restart_step_rank_"+s.str();
    if (rank==0){
        if (!boost::filesystem::exists(filename)){
            cout<<"No restarts found for "<<filename<<endl;
            error=1;
        }
    }
    MPI_Bcast(&error, 1, MPI_INT, 0 , MPI_COMM_WORLD);
    if (error==0){
        int nxtID[numProcs];
        for (int i=0;i<numProcs;i++)nxtID[i]=0;
        while(boost::filesystem::exists(filename)){
            
            if (rank == 0){
                
                if (error==0)cout<<"Reading restarts from "<<filename<<endl;
                std::ifstream ifs(filename);
                {// this block ensures the archive gets closed on block exit
                    try {
                        if (_archiveFormat =="binary"){boost::archive::binary_iarchive ia(ifs);ia>>_packages;}
                        if (_archiveFormat =="text"  ){boost::archive::text_iarchive   ia(ifs);ia>>_packages;}
                        
                        cout<<"Read in "<<_packages.size()<<" mad agents..."<<endl;
                        for (auto& p:_packages){
                            repast::AgentId id=p.getId();
                            int n=Cohort::_NextID;
                            //Cohort::_NextID should be incremented here so set increase flag to true
                            MadAgent* a=NULL;
                            if(id.agentType()==_cohortType){
                                a=new Cohort( id,p,true );
                            }
                            if(id.agentType()==_stockType) {
                                a=new Stock( id,p );
                            }
                            if (a!=NULL){
                                //moveTo checks if we are local - agents may previously have been on another core, so make sure local now
                                a->getId().currentRank(0);
                                //agent unique Ids depend on id no. type and startingRank - to make sure any new agents are unique
                                //need to find the max id no. associated with a given previous rank and add 1
                                //number of threads may have decreased so check to be sure we don't overflow nxtID
                                if (a->getId().startingRank()<numProcs)
                                    nxtID[a->getId().startingRank()]=max<int>(a->getId().id()+1,nxtID[a->getId().startingRank()]);
                                _context.addAgent(a);
                                repast::Point<int> initialLocation(int(a->getLocation()[0]),int(a->getLocation()[1]));
                                space()->moveTo(id, initialLocation);
                            }else{ cout<<"Warning NULL agents on reading restart file "<<filename<<endl;}
                        }
                        _packages.clear();
                    } catch (exception e){
                        if (r==0)cout <<"************* Error attempting to open archive: wrong file format? *************";
                        error=1;
                    }
                }
            }
            
            //check for errors
            MPI_Bcast(&error, 1, MPI_INT, 0 , MPI_COMM_WORLD);
            //sync so thread 0 doesn't have to carry all the agents from a possibly multi-core previous run
            //remember every thread needs to do the sync, not just thread 0!
            sync();
            //share the nextID data - this may update progressively, as startingRanks of current live agents can be arbitrary
            //if there are more threads than previously, those with rank not previously present can have nextID=0
            //if there are fewer threads, those with startingRank>=numProcs are safe to ignore (as there will be no new agents with these startingRanks)
            MPI_Bcast(&nxtID, numProcs, MPI_INT, 0 , MPI_COMM_WORLD);
            Cohort::_NextID=nxtID[rank];
            
            r++;
            std::stringstream s;
            s<<step<<"_"<<r;
            filename=_restartDirectory+"Restart_step_rank_"+s.str();
        }
        
    }

    if (error!=0){
            MPI_Finalize();
            exit(error);
    }
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
    assert(Parameters::instance()->GetLengthLongitudeArray( )==_maxX-_minX+1);
    assert(Parameters::instance()->GetLengthLatitudeArray( )==_maxY-_minY+1);
    //now set up the environmental cells - note at present this uses the full grid, not just local to this thread
    //so that off-thread environment can be easily queried. Currently some duplication here, but it is not a huge amount of data.
    _Env=Environment(_minX,_maxX,_minY,_maxY);

    //get the definitions of stocks and cohorts
    StockDefinitions::Initialise(_props->getProperty("input.DataDirectory")+"/"+_props->getProperty("input.StockDefinitionsFileName"));
    CohortDefinitions::Initialise(_props->getProperty("input.DataDirectory")+"/"+_props->getProperty("input.CohortDefinitionsFileName"));
    int x=_xlo,y=_ylo;
    repast::Point<int> initialLocation(x,y);
    repast::Point<int> origin(_minX,_minY);
    EnvironmentCell* E=_Env[x][y];
    
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
         if (nr<int(random->GetUniform()*10.)){a->_alive=false;_context.removeAgent(a->getId());nr++;}//care here - after the removal, the pointer to a is no longer valid
         else{
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

   //if there's a buffer, and you change local agent data, seems you need a sync of states *before* you do anything like a move...
   //sync after a move leads to errors
   repast::RepastProcess::instance()->synchronizeAgentStates<AgentPackage,MadAgentPackageProvider, MadAgentPackageReceiver>(*provider, *receiver);
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
    //cout<<Cohort::_AdvectionTimeStepsPerModelTimeStep<<endl;
    //assert(40== Cohort::_AdvectionTimeStepsPerModelTimeStep);
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
         AgentPackage package(a->getId());
         ((Cohort*)a)->PushThingsIntoPackage( package );
         _packages.push_back(package);
      }
      oa<<_packages;
     }
     cout<<"Wrote "<<_packages.size()<<" cohorts"<<endl;
     _packages.clear();
     std::ifstream ifs("TestAgentSerialization_rank_"+s.str());

     {
      boost::archive::text_iarchive ia(ifs);
      ia>>_packages;
      cout<<"Read "<<_packages.size()<<" cohorts"<<endl;
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
/*    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    if (rank==0)cout<<"Test15: check direct agent dispersal subroutine works (movement by fraction of a cell permitted)"<<endl;

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
      E=_Env[x][y];
      c->setup(0,1, E,random);
      _context.addAgent(c);
      discreteSpace->moveTo(id, initialLocation);
      c->setLocation(x,y);c->setDestination(x,y);
      double uSpeed=0.9*E->Width(),vSpeed=0.5*E->Height();
      c->TryToDisperse(uSpeed, vSpeed,E, this);
      if (c->_moved){
//        assert(c->_destination[0]-x-0.9<1.e-14);
//        assert(c->_destination[1]-y-0.5<1.e-14);
        vector<int>newPlace={0,0};
        newPlace[0]=int(c->_destination[0]);
        newPlace[1]=int(c->_destination[1]);
        c->_location=c->_destination;
      }
          
     }
     //if you have a buffer you have to sync agent states when agent data has changed, before you can do a move.
     repast::RepastProcess::instance()->synchronizeAgentStates<AgentPackage,MadAgentPackageProvider, MadAgentPackageReceiver>(*provider, *receiver);
     agents.clear();
     _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
     for (auto a:agents){
      if (a->_moved){
        vector<int>newPlace={int(a->_destination[0]),int(a->_destination[1])};
        space()->moveTo(a,newPlace);
      }
     }
     sync();
     agents.clear();
     _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
     //nothing should have moved rank at this point
     for (auto a:agents){
         assert(a->getId().currentRank()==a->getId().startingRank());
//         assert(((Cohort*)a)->_location[0]-x-0.9<1.e-14);
//         assert(((Cohort*)a)->_location[1]-y-0.5<1.e-14);
         vector<int>newPlace={_minX,(_maxY + _minY)/2};

         ((Cohort*)a)->setDestination(newPlace[0],newPlace[1]);
         ((Cohort*)a)->setLocation(newPlace[0],newPlace[1]);

     }
     //if you have a buffer you have to sync agent states when agent data has changed, before you can do a move.
     repast::RepastProcess::instance()->synchronizeAgentStates<AgentPackage,MadAgentPackageProvider, MadAgentPackageReceiver>(*provider, *receiver);

     for (auto a:agents){
        vector<int>newPlace={int(a->_destination[0]),int(a->_destination[1])};
        space()->moveTo(a,newPlace);
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
       for (int i=0; i<1;i++){
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
//      assert(mean-(_maxX - _minX)<0.02);
//      assert(stdev-2<0.01);
      cout<<"Test 15 : succeeded"<<endl;
    }*/
    //---------------------------------------------------
    //***-------------------TEST 16-------------------***//
    //---------------------------------------------------
    //test whether advective dispersal works.
    agents.clear();
    _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
    if (rank==0)cout<<"Test16: test dispersal works for ocean plankton"<<endl;

    for (auto a:agents){
        _context.removeAgent(a->getId());
    }
    sync();
    _eating          = false;
    _metabolism      = false;
    _reproduction    = false;
    _death           = false;
    _mergers         = false;
    _dispersal       = true;
    _output          = false;
    _noLongitudeWrap = false;
    _dispersalSelection=_props->getProperty("simulation.DispersalSelection");
    //now add a new agent on each thread and set its properties to known values
    
    for (int i=0;i<1;i++){
      repast::AgentId id(Cohort::_NextID, rank, _cohortType);
      id.currentRank(rank);
      Cohort* c = new Cohort(id);
      E=_Env[x][y];
      c->setup(9,1, E,random);
      _context.addAgent(c);
      discreteSpace->moveTo(id, initialLocation);
      c->setLocation(_maxX,y);c->setDestination(_maxX,y);
      vector<int>newPlace={x,y};
      space()->moveTo(c,newPlace);
      assert(c->_IsPlanktonic);
    }
    sync();    
    for (int i=0;i<1200;i++){
         step();

         agents.clear();
         _context.selectAgents(repast::SharedContext<MadAgent>::LOCAL,agents);
         for (auto a:agents){
            assert(((Cohort *)a)->_IsPlanktonic);
            if (a->getId().id()==3015)cout<<rank<<" "<<" "<<a->getId().id()<<" "<<a->_location[0]<<" "<<a->_location[1]<<endl;
         }
    }

     // cout<<"Test 16 : succeeded"<<endl;
    
}    
//---------------------------------------------------------------------------------------------------------------------------
//define some data values for the Cohort to check whether they are preserved on moving across threads
//see Cohort::setup
void MadModel::setupCohortTestValues(Cohort* c){
    

    c->_FunctionalGroupIndex=floor(repast::Random::instance()->nextDouble()*CohortDefinitions::Get()->size());//functionalGroup;
    c->_Merged                      = false;
    c->_alive                       = true;
    
    c->setPropertiesFromCohortDefinitions(c->_FunctionalGroupIndex);

    c->_IsMature=false;
    c->_IndividualReproductivePotentialMass=0;

    c->_BirthTimeStep=0;
    c->_MaturityTimeStep=std::numeric_limits<unsigned>::max( );


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
    //assert(c->_FunctionalGroupIndex==12);
    assert(c->_Merged                      == false);
    assert(c->_alive                       == true);

    assert(c->_Heterotroph==(CohortDefinitions::Get()->Trait(c->_FunctionalGroupIndex     , "heterotroph/autotroph")  =="heterotroph"));
    assert(c->_Autotroph  ==!c->_Heterotroph);
    assert(c->_Endotherm  ==(CohortDefinitions::Get()->Trait(c->_FunctionalGroupIndex     , "endo/ectotherm")         =="endotherm"));
    assert(c->_Ectotherm  ==!c->_Endotherm);
    assert(c->_Realm      ==CohortDefinitions::Get()->Trait(c->_FunctionalGroupIndex     , "realm"));

    assert(c->_Iteroparous==(CohortDefinitions::Get()->Trait(c->_FunctionalGroupIndex     , "reproductive strategy")  =="iteroparity"));
    assert(c->_Semelparous==!c->_Iteroparous);
    assert(c->_Herbivore==(CohortDefinitions::Get()->Trait(c->_FunctionalGroupIndex       , "nutrition source")       =="herbivore"));
    assert(c->_Carnivore==(CohortDefinitions::Get()->Trait(c->_FunctionalGroupIndex       , "nutrition source")       =="carnivore"));
    assert(c->_Omnivore== (CohortDefinitions::Get()->Trait(c->_FunctionalGroupIndex       , "nutrition source")       =="omnivore"));
    assert(c->_IsPlanktonic== (CohortDefinitions::Get()->Trait(c->_FunctionalGroupIndex   , "mobility")               =="planktonic"));
    assert(c->_IsFilterFeeder==(CohortDefinitions::Get()->Trait(c->_FunctionalGroupIndex   , "diet")                   =="allspecial"));
    
    assert(c->_ProportionSuitableTimeActive== CohortDefinitions::Get()->Property(c->_FunctionalGroupIndex   ,"proportion suitable time active"));
    
    assert(c->_IsMature==false);
    assert(c->_IndividualReproductivePotentialMass==0);
    
    assert(c->_AssimilationEfficiency_H==CohortDefinitions::Get()->Property(c->_FunctionalGroupIndex   ,"herbivory assimilation"));
    assert(c->_AssimilationEfficiency_C==CohortDefinitions::Get()->Property(c->_FunctionalGroupIndex   ,"carnivory assimilation"));
    assert(c->_BirthTimeStep==0);
    assert(c->_MaturityTimeStep==std::numeric_limits<unsigned>::max( ));
    assert(c->_MinimumMass==CohortDefinitions::Get()->Property(c->_FunctionalGroupIndex   ,"minimum mass"));
    assert(c->_MaximumMass==CohortDefinitions::Get()->Property(c->_FunctionalGroupIndex   ,"maximum mass"));

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
