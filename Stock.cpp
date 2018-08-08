/*

 *
 * Stock.cpp
 *
 *  Created on: May 7, 2018
 *      Author: Mike Bithell
 * 
 *  Derived from Origianl C# code by
 * Emergent Global Patterns of Ecosystem Structure and Function from a Mechanistic General Ecosystem Model , Michael B. J. Harfoot , Tim Newbold, Derek P. Tittensor,  Stephen Emmott, Jon Hutton, Vassily Lyutsarev, Matthew J. Smith, Jörn P. W. Scharlemann, Drew W. Purves PLOS Biology 12(4) 2014 e1001841 https://doi.org/10.1371/journal.pbio.1001841
 */

#include "Stock.h"
#include "Cohort.h"

#include "Environment.h"
#include "AutotrophProcessor.h"
#include "TerrestrialCarbon.h"
#include "HANPP.h"
#include "Groups.h"


//-------------------------------------------------------------------------------------------------------------------
 
void Stock::setup(unsigned functionalGroup,Environment* LocalEnvironment){
     _FunctionalGroupIndex = functionalGroup;
    // Get the individual body masses for organisms in each stock functional group

     _IndividualBodyMass = StockDefinitions::Get()->Property(functionalGroup, "individual mass" );
    
    _Marine              =(StockDefinitions::Get()->Trait(functionalGroup   , "realm")         =="marine");
    _Deciduous           =(StockDefinitions::Get()->Trait(functionalGroup   , "leaf strategy") =="deciduous");

    // If it is a functional group that corresponds to the current realm, then seed the stock
    if( ! _Marine && LocalEnvironment->Precipitation() != Constants::cMissingValue && LocalEnvironment->Temperature()!= Constants::cMissingValue ) {
            // An instance of the terrestrial carbon model class
            TerrestrialCarbon PlantModel;

            // Calculate predicted leaf mass at equilibrium for this stock
            _TotalBiomass = PlantModel.CalculateEquilibriumLeafMass(LocalEnvironment , _Deciduous );

            
    } else if( _Marine && LocalEnvironment->NPP()!= Constants::cMissingValue ) {
            _TotalBiomass = 1.e12;
    }

}
//-------------------------------------------------------------------------------------------------------------------
void Stock::step(double& AllBiomass,Environment* LocalEnvironment,const unsigned CurrentTimeStep) {


    if( _Marine ) {
        AutotrophProcessor A;
        double NPP=A.ConvertNPPToAutotroph(LocalEnvironment);

        _TotalBiomass += NPP*Constants::cDay;

        // If the biomass of the autotroph stock has been made less than zero (i.e. because of negative NPP) then reset to zero
        if( _TotalBiomass < 0.0 ) _TotalBiomass = 0.0;

        
    } else {
        TerrestrialCarbon PlantModel;
        HumanAutotrophMatterAppropriation HumanRemoval;
        // Run the dynamic plant model to update the leaf stock for this time step
        double NPPWetMatter = PlantModel.UpdateLeafStock( LocalEnvironment ,this, _Deciduous );//??could pass TotalBiomass by reference?

        // Apply human appropriation of NPP - note in the latest C# version this is changed to include the NPPWetMatter calculated above
        // AllBiomass is the current amount in this cell, so adjust for HANPP as well.
        double FracBiomass=_TotalBiomass/AllBiomass;
        AllBiomass-=_TotalBiomass;
        double fhanpp = HumanRemoval.RemoveHumanAppropriatedMatter(LocalEnvironment, NPPWetMatter, FracBiomass, CurrentTimeStep);
        _TotalBiomass += NPPWetMatter * ( 1 - fhanpp );
        AllBiomass+=_TotalBiomass;
    }

}
