#ifndef COLLECTIVES_H
#include agent.h
#include <string>
class Ateam{
//in MadAgent - memberOf, leaderOf
//designed to be able to be a hierarchy - each agent in the member can be a lead of a sub-team
//agents can create temporary teams with no leader (but then might vote in or choose a leader or spokesperson), create a team as leader with no memebers and then recruit, 
MadAgent* leader;
vector<MadAgent*>members;
};
class Mother{
//in MadAgent motherOf (+parentOf?fatherOf?). Index by mother as this is data the mother (and others) can be sure of (!)
//link to relationships for birth probabilties - NB a "family" could be all genetic relatives (but this is all people eventually! a family of strangers) or maybw know relatives?
//what about adtoption and legitmacy?
vector <MadAgent*> children;
vector <MadAgent*> carers;
MadAgent* mother,father;
};
class place{
 int name;//name is an index into lcation names; places can be points, polygons, rasters, buildings, or other fixed-location objects.
};
class country {
 place* location;//might be a point- but places can be various types - could be the country polygon?
};
class nationalGovernment
 country* home;//what about local government or city states?
};
class productionSystem{
//e.g., industries, hunting, gathering, farming, sleep, religion (outputs wellbeing? social capital?), search (for places, outputs could be place and route)...maybe this one is more complex
    int name;//index for what kind of production system this is - maybe only need one prototype instantiated that gets used by all agents? may need to store system states though (how far through production, how much produced?)
    vector<int> inputs;//index for the input type -can production systems be hierarchies/networks(production system indexes another production system)? what about technolgies and efficiency?
    vector<int> inputQuantityPerUnitOutput//one per input type - units need thinking about (e.g. land area for farming)
    int output;//index into output type. multiple outputs? what about something like a hospital (could be a prodcution system network?). index by output(s) for goal setting 
    int outputQuality;//index into some kind of quality structure?
    vector<int> waste;//index of waste types
    vector<int> wastePerUnitOutput;//quantity of each waste type. implies a requirement for storage/disposal?
    vector<int> energyInputs;//may need fuel+electricity, for example.
    vectro<int> energyRequirementPerUnitOutput;//
    int productionTimePerUnitOutput;//these next two might not be known ahead of time? (e.g. for hunting).
    int labourRequirementPerUnitOutput;//they also depend on efficiency - labour productivity somehow? 
    //need capital as well as inputs?
    //agents needs can be met by a search for suitable production system, then a team created to implement. Other existing production systems need to be known about for inputs.
};
#endif
