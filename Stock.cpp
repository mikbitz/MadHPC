/*
 *   Repast for High Performance Computing (Repast HPC)
 *
 *   Copyright (c) 2010 Argonne National Laboratory
 *   All rights reserved.
 *
 *   Redistribution and use in source and binary forms, with
 *   or without modification, are permitted provided that the following
 *   conditions are met:
 *
 *     Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *
 *     Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *
 *     Neither the name of the Argonne National Laboratory nor the names of its
 *     contributors may be used to endorse or promote products derived from
 *     this software without specific prior written permission.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *   PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE TRUSTEES OR
 *   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 *   EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 * Stock.cpp
 *
 *  Created on: Sep 1, 2010
 *      Author: nick
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
