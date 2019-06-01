/*
 *
 * Human.cpp
 *
 *  Created on: Aug 1, 2018
 *      Author: Mike
 *  
 */


#include "repast_hpc/RepastProcess.h"
#include "repast_hpc/Random.h"

#include "Groups.h"
#include "Human.h"
#include "Stock.h"
#include "Parameters.h"
#include "model.h"

using namespace repast;

//needs to be class variable so that Cohorts can have unique numbers
unsigned Human::_NextID=0;
//------------------------------------------------------------------------------------------------------------
void Human::setup(unsigned functionalGroup,unsigned numHumansThisCell,EnvironmentCell* e,randomizer* r){
    //ResetAccounts( );

    _FunctionalGroupIndex=functionalGroup;
    _Merged                      = false;
    _alive                       = true;

    _Realm      ="terrestrial";

    _ProportionSuitableTimeActive= 0.5;
    
    _IsMature=false;
    _IndividualReproductivePotentialMass=0;
    
    _AssimilationEfficiency_H=0.4; 
    _AssimilationEfficiency_C=0.64;
    _BirthTimeStep=0;
    _MaturityTimeStep=std::numeric_limits<unsigned>::max( );
    _MinimumMass=184800;  //6st the minimum adult mass in this cohort
    _MaximumMass=1848000;//60st - possible...

  
    _AdultMass = pow( 10, ( repast::Random::instance()->nextDouble() * ( log10( _MaximumMass ) - log10( 50 * _MinimumMass ) ) + log10( 50 * _MinimumMass ) ) );
    
    NormalGenerator NJ= repast::Random::instance()->createNormalGenerator(0.1,0.02);
    
    _LogOptimalPreyBodySizeRatio = log(std::max( 0.01, NJ.next() ));
    
    double expectedLnAdultMassRatio = 2.24 + 0.13 * log( _AdultMass );
    //in the original code the mean and sd are those of the underlying normal distribution
    //in the boost library they refer to the log distibution - see
    //https://www.boost.org/doc/libs/1_43_0/libs/math/doc/sf_and_dist/html/math_toolkit/dist/dist_ref/dists/lognormal_dist.html
    //NB indepedent calls to create an RNG generator will follow the correct random sequence (i.e. they will not just re-create the initial random value) as this appears to use a singleton (tested -OK!)
    LogNormalGenerator LNJ= repast::Random::instance()->createLogNormalGenerator(exp(expectedLnAdultMassRatio+0.5*0.5/2.), (exp(0.5*0.5)-1)*exp(2*expectedLnAdultMassRatio+0.5*0.5));

    do {
        _JuvenileMass = _AdultMass  / (1.0 + LNJ.next());
    } while( _AdultMass <= _JuvenileMass || _JuvenileMass < _MinimumMass );


    double NewBiomass = ( 3300. / numHumansThisCell ) * 100 * 3000 * pow( 0.6, ( log10( _JuvenileMass ) ) ) * ( e->Area() );

    _CohortAbundance = NewBiomass / _JuvenileMass;
    _MaximumAchievedBodyMass=_JuvenileMass;
    _IndividualBodyMass=_JuvenileMass;
    _moved=false;
    _location={0,0};
}
//------------------------------------------------------------------------------------------------------------
void Human::setupOffspring( Human* actingHuman, double juvenileBodyMass, double adultBodyMass, double initialBodyMass, double initialAbundance, unsigned birthTimeStep ) {
     
     _newH=NULL;
    _FunctionalGroupIndex        = actingHuman->_FunctionalGroupIndex;
    _JuvenileMass                = juvenileBodyMass;
    _AdultMass                   = adultBodyMass;
    _IndividualBodyMass          = initialBodyMass;
    _CohortAbundance             = initialAbundance;
    _CurrentTimeStep               = birthTimeStep;
    _BirthTimeStep               = birthTimeStep;
    _MaturityTimeStep            = std::numeric_limits<unsigned>::max( );
    _LogOptimalPreyBodySizeRatio = actingHuman->_LogOptimalPreyBodySizeRatio;
    _MaximumAchievedBodyMass     = juvenileBodyMass;
    _Merged                      = false;
    _alive                       = true;
    _IndividualReproductivePotentialMass = 0;
  
    _Realm      =    actingHuman->_Realm;

    _MinimumMass=    actingHuman->_MinimumMass;
    _MaximumMass=    actingHuman->_MaximumMass;
    
    _ProportionSuitableTimeActive= actingHuman->_ProportionSuitableTimeActive;
    
    _IsMature=false;
    
    _AssimilationEfficiency_H=0.4; 
    _AssimilationEfficiency_C=0.64;
    
    _moved=false;
    _location=actingHuman->_location;
    _destination=_location;


}
//------------------------------------------------------------------------------------------------------------
//Required by RHPC for cross-core copy - NB "Accounts" do not need to be included as they are instantaneous within a timestep
void Human::PullThingsOutofPackage( const AgentPackage& package ) {

    _FunctionalGroupIndex        = package._contents._FunctionalGroupIndex;
    _JuvenileMass                = package._contents._JuvenileMass;
    _AdultMass                   = package._contents._AdultMass;
    _IndividualBodyMass          = package._contents._IndividualBodyMass;
    _CohortAbundance             = package._contents._CohortAbundance;
    _BirthTimeStep               = package._contents._BirthTimeStep;
    _MaturityTimeStep            = package._contents._MaturityTimeStep;
    _LogOptimalPreyBodySizeRatio = package._contents._LogOptimalPreyBodySizeRatio;
    _MaximumAchievedBodyMass     = package._contents._MaximumAchievedBodyMass;
    _Merged                      = package._contents._Merged;
    _alive                       = package._contents._alive;
    _IndividualReproductivePotentialMass = package._contents._IndividualReproductivePotentialMass ;
  
    _Realm      ="terrestrial";

    _ProportionSuitableTimeActive= 0.5;   
    _AssimilationEfficiency_H=0.4; 
    _AssimilationEfficiency_C=0.64;
    _MinimumMass=184800;  //6st the minimum adult mass in this cohort
    _MaximumMass=1848000;//60st - possible...
    
    _IsMature=package._contents._IsMature;
    
    
    _moved=package._contents._moved;
    _location=package._contents._location;
    _destination=package._contents._destination;
    
}
//------------------------------------------------------------------------------------------------------------
//Required by RHPC for cross-core copy
void Human::PushThingsIntoPackage( AgentPackage& package ) {

    package._contents._FunctionalGroupIndex        =  _FunctionalGroupIndex;
    package._contents._JuvenileMass                =  _JuvenileMass;
    package._contents._AdultMass                   =  _AdultMass;
    package._contents._IndividualBodyMass          =  _IndividualBodyMass;
    package._contents._CohortAbundance             =  _CohortAbundance;
    package._contents._BirthTimeStep               =  _BirthTimeStep;
    package._contents._MaturityTimeStep            =  _MaturityTimeStep;
    package._contents._LogOptimalPreyBodySizeRatio =  _LogOptimalPreyBodySizeRatio;
    package._contents._MaximumAchievedBodyMass     =  _MaximumAchievedBodyMass;
    package._contents._Merged                      =  _Merged;
    package._contents._alive                       =  _alive;
    package._contents._IndividualReproductivePotentialMass =  _IndividualReproductivePotentialMass ;
  
    
    package._contents._IsMature= _IsMature;
       
    package._contents._moved=_moved;
    package._contents._location=_location;
    package._contents._destination=_destination;
    
}

//------------------------------------------------------------------------------------------------------------
void Human::step(EnvironmentCell* e,vector<Human*>& preys,vector<Stock*>& stocks,const unsigned Timestep,MadModel* m) {
    _newH=NULL;//make sure the reproduction pointer has been zeroed out

    if (_CohortAbundance - Parameters::Get( )->GetExtinctionThreshold( ) <= 0)return;
    _CurrentTimeStep=Timestep;

    //note - passing in preys here from above ensures that new cohorts created later do not immediately get eaten.
    //since they *do* get added to the global cohort list straight away.
    //pass in of environment avoids having to find it with an expensive lookup
    //ResetAccounts( );
    //for (auto A: _Accounting)for (auto B: A.second)if(B.first!="mortality")assert(B.second==0);
    //assignTimeActive(e);
    //eat(e,preys,stocks,m);
    //metabolize(e);
    //reproduce(e);
    //mort();
    //applyEcology(e);

}











