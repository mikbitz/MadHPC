/* Model.cpp */

#include <stdio.h>
#include <vector>
#include <boost/mpi.hpp>
#include "repast_hpc/AgentId.h"
#include "repast_hpc/RepastProcess.h"
#include "repast_hpc/Utilities.h"
#include "repast_hpc/Properties.h"


#include "repast_hpc/initialize_random.h"
#include "repast_hpc/Point.h"
#include "model.h"
#include "agent.h"
//------------------------------------------------------------------------------------------------------------
//Constructor and destructor
//------------------------------------------------------------------------------------------------------------
RepastHPCModel::RepastHPCModel(std::string propsFile, int argc, char** argv, boost::mpi::communicator* comm): context(comm){
    //-----------------
    //Pull in parameter from the model.props file
    props = new repast::Properties(propsFile, argc, argv, comm);
    //Number of timesteps
	stopAt = repast::strToInt(props->getProperty("stop.at"));
    //extent of buffer zones in grid units - this many grid cells are sahred at the boundary between cores
    gridBuffer = repast::strToInt(props->getProperty("grid.buffer"));
    //Total agents in the run
	countOfAgents = repast::strToInt(props->getProperty("count.of.agents"));
    //Grid extent
    minX=repast::strToInt(props->getProperty("min.x"));
    minY=repast::strToInt(props->getProperty("min.y"));
    maxX=repast::strToInt(props->getProperty("max.x"));
    maxY=repast::strToInt(props->getProperty("max.y"));
    int dimX=repast::strToInt(props->getProperty("proc.per.x"));
    int dimY=repast::strToInt(props->getProperty("proc.per.y"));
    //-----------------
    //initialize default random number generator
    initializeRandom(*props, comm);
    //-----------------
    //save out the properties for this run in the output
    if(repast::RepastProcess::instance()->rank() == 0) props->writeToSVFile("./output/record.csv");
    //-----------------
	//create the model grid
    repast::Point<double> origin(minX,minY);
    repast::Point<double> extent(maxX-minX+1, maxY-minY+1);
    
    repast::GridDimensions gd(origin, extent);
    
    std::vector<int> processDims;
    processDims.push_back(dimX);
    processDims.push_back(dimY);
  
    discreteSpace = new repast::SharedDiscreteSpace<RepastHPCAgent, repast::WrapAroundBorders, repast::SimpleAdder<RepastHPCAgent> >("AgentDiscreteSpace", gd, processDims, gridBuffer, comm);
	
   
    //The agent container is a context. Add the grid to it.
   	context.addProjection(discreteSpace);
    
   //-----------------
    //Set up the cross-thread data transferers
	provider = new RepastHPCAgentPackageProvider(&context);
	receiver = new RepastHPCAgentPackageReceiver(&context);
    
    //-----------------
	// Set up Data collection
	// Create the data set builder
	std::string fileOutputName("./output/agent_total_data.csv");
	repast::SVDataSetBuilder builder(fileOutputName.c_str(), ",", repast::RepastProcess::instance()->getScheduleRunner().schedule());
	
	// Create the individual data sets to be added to the builder
	DataSource_AgentTotals* agentTotals_DataSource = new DataSource_AgentTotals(&context);
	builder.addDataSource(createSVDataSource("Total", agentTotals_DataSource, std::plus<int>()));

	DataSource_AgentCTotals* agentCTotals_DataSource = new DataSource_AgentCTotals(&context);
	builder.addDataSource(createSVDataSource("C", agentCTotals_DataSource, std::plus<int>()));

	// Use the builder to create the data set
	agentValues = builder.createDataSet();
	
}
//------------------------------------------------------------------------------------------------------------
RepastHPCModel::~RepastHPCModel(){
	delete props;
	delete provider;
	delete receiver;

}
//------------------------------------------------------------------------------------------------------------
//Model Initialisation
//------------------------------------------------------------------------------------------------------------
void RepastHPCModel::initSchedule(repast::ScheduleRunner& runner){
	runner.scheduleEvent(1, 1, repast::Schedule::FunctorPtr(new repast::MethodFunctor<RepastHPCModel> (this, &RepastHPCModel::step)));
	runner.scheduleEndEvent(repast::Schedule::FunctorPtr(new repast::MethodFunctor<RepastHPCModel> (this, &RepastHPCModel::recordResults)));
	runner.scheduleStop(stopAt);
}
void RepastHPCModel::init(){
    int rank = repast::RepastProcess::instance()->rank();
    if (rank==0) cout<<"In init "<<discreteSpace->dimensions().extents().getX()<<" "<<discreteSpace->dimensions().extents().getY()<<endl;
    //explicitly use the local bounds of the grid on this thread to create countOfAgents per cell.
    //Not doing this can lead to problems with agents in dsitant cells not within the local thread neighbourhood
    //see SharedBaseGrid.h moveTo method
    int xmin=        discreteSpace->dimensions().origin().getX() ;
    int xmax= xmin + discreteSpace->dimensions().extents().getX();
    int ymin=        discreteSpace->dimensions().origin().getY() ;
    int ymax= ymin + discreteSpace->dimensions().extents().getY();
    int c=0;
	for(int i = 0; i < countOfAgents; i++){
        for(int x = xmin; x < xmax; x++){
            for(int y = ymin; y < ymax; y++){
             repast::Point<int> initialLocation(x,y);
             //make sure the agentId is unique on this thread!!
             // values are int id, int startProc, int agentType, 
		     repast::AgentId id(c, rank, 0);
		     id.currentRank(rank);
		     RepastHPCAgent* agent = new RepastHPCAgent(id);
		     context.addAgent(agent);
             discreteSpace->moveTo(id, initialLocation);
             c++;
            }
        }
	}
}
//------------------------------------------------------------------------------------------------------------
//Run the model
//------------------------------------------------------------------------------------------------------------
void RepastHPCModel::step(){
	if(repast::RepastProcess::instance()->rank() == 0) 
        std::cout << " TICK " << repast::RepastProcess::instance()->getScheduleRunner().currentTick() << std::endl;
    sync();

	std::vector<RepastHPCAgent*> agents;
    context.selectAgents(repast::SharedContext<RepastHPCAgent>::LOCAL,agents);

	//for(auto& a:agents){
	//	a->step(this);
    //    }
	for (int i=1;i<3;i++){
    for (auto& a:agents){
 	    a->move(this,-1,-1);
        }
    }
    for(auto& a:agents){
        if (a->getId().startingRank()==0 && a->getId().id()==25)std::cout << " AGENT " << " ON Thread " << repast::RepastProcess::instance()->rank() << std::endl;

		a->reportLocation(this);
        } 
   
}
//------------------------------------------------------------------------------------------------------------
void RepastHPCModel::sync(){
    //These lines synchronize the agents across all threads - if there is more than one...
	discreteSpace->balance();
    repast::RepastProcess::instance()->synchronizeAgentStatus<RepastHPCAgent, RepastHPCAgentPackage, 
             RepastHPCAgentPackageProvider, RepastHPCAgentPackageReceiver>(context, *provider, *receiver, *receiver);
    
    repast::RepastProcess::instance()->synchronizeProjectionInfo<RepastHPCAgent, RepastHPCAgentPackage, 
             RepastHPCAgentPackageProvider, RepastHPCAgentPackageReceiver>(context, *provider, *receiver, *receiver);

	repast::RepastProcess::instance()->synchronizeAgentStates<RepastHPCAgentPackage, 
             RepastHPCAgentPackageProvider, RepastHPCAgentPackageReceiver>(*provider, *receiver);
}
//------------------------------------------------------------------------------------------------------------
// Packages for exchanging agents across threads
//------------------------------------------------------------------------------------------------------------


RepastHPCAgentPackageProvider::RepastHPCAgentPackageProvider(repast::SharedContext<RepastHPCAgent>* agentPtr): agents(agentPtr){ }
//------------------------------------------------------------------------------------------------------------

void RepastHPCAgentPackageProvider::providePackage(RepastHPCAgent * agent, std::vector<RepastHPCAgentPackage>& out){
    repast::AgentId id = agent->getId();
    RepastHPCAgentPackage package(id.id(), id.startingRank(), id.agentType(), id.currentRank(), agent->getC(), agent->getTotal());
    out.push_back(package);
}
//------------------------------------------------------------------------------------------------------------

void RepastHPCAgentPackageProvider::provideContent(repast::AgentRequest req, std::vector<RepastHPCAgentPackage>& out){
    std::vector<repast::AgentId> ids = req.requestedAgents();
    for(size_t i = 0; i < ids.size(); i++){
        providePackage(agents->getAgent(ids[i]), out);
    }
}
//------------------------------------------------------------------------------------------------------------


RepastHPCAgentPackageReceiver::RepastHPCAgentPackageReceiver(repast::SharedContext<RepastHPCAgent>* agentPtr): agents(agentPtr){}
//------------------------------------------------------------------------------------------------------------

RepastHPCAgent * RepastHPCAgentPackageReceiver::createAgent(RepastHPCAgentPackage package){
    repast::AgentId id(package.id, package.rank, package.type, package.currentRank);
    return new RepastHPCAgent(id, package.c, package.total);
}
//------------------------------------------------------------------------------------------------------------

void RepastHPCAgentPackageReceiver::updateAgent(RepastHPCAgentPackage package){
    repast::AgentId id(package.id, package.rank, package.type);
    RepastHPCAgent * agent = agents->getAgent(id);
    agent->set(package.currentRank, package.c, package.total);
}

//------------------------------------------------------------------------------------------------------------
// DATA SAVING SECTION
//------------------------------------------------------------------------------------------------------------

void RepastHPCModel::recordResults(){
	if(repast::RepastProcess::instance()->rank() == 0){
		props->putProperty("Result","Passed");
		std::vector<std::string> keyOrder;
		keyOrder.push_back("RunNumber");
		keyOrder.push_back("stop.at");
		keyOrder.push_back("Result");
		props->writeToSVFile("./output/results.csv", keyOrder);
    }
}

//------------------------------------------------------------------------------------------------------------

DataSource_AgentTotals::DataSource_AgentTotals(repast::SharedContext<RepastHPCAgent>* c) : context(c){ }

int DataSource_AgentTotals::getData(){
	int sum = 0;
	repast::SharedContext<RepastHPCAgent>::const_local_iterator iter    = context->localBegin();
	repast::SharedContext<RepastHPCAgent>::const_local_iterator iterEnd = context->localEnd();
	while( iter != iterEnd) {
		sum+= (*iter)->getTotal();
		iter++;
	}
	return sum;
}
//------------------------------------------------------------------------------------------------------------


DataSource_AgentCTotals::DataSource_AgentCTotals(repast::SharedContext<RepastHPCAgent>* c) : context(c){ }

int DataSource_AgentCTotals::getData(){
	int sum = 0;
	repast::SharedContext<RepastHPCAgent>::const_local_iterator iter    = context->localBegin();
	repast::SharedContext<RepastHPCAgent>::const_local_iterator iterEnd = context->localEnd();
	while( iter != iterEnd) {
		sum+= (*iter)->getC();
		iter++;
	}
	return sum;
}
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

