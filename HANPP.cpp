/*
 *
 * HANPP.cpp
 *
 *  Created on: May 7, 2018
 *      Author: Mike Bithell
 * 
 *  Derived from Origianl C# code by
 * Emergent Global Patterns of Ecosystem Structure and Function from a Mechanistic General Ecosystem Model , Michael B. J. Harfoot , Tim Newbold, Derek P. Tittensor,  Stephen Emmott, Jon Hutton, Vassily Lyutsarev, Matthew J. Smith, Jörn P. W. Scharlemann, Drew W. Purves PLOS Biology 12(4) 2014 e1001841 https://doi.org/10.1371/journal.pbio.1001841
 */
#include "HANPP.h"
#include "Parameters.h"

/** \file HANPP.cpp
 * \brief the HANPP implemetation file
 */

    HumanAutotrophMatterAppropriation::HumanAutotrophMatterAppropriation( ) {
    }
    //----------------------------------------------------------------------------------------------
    double HumanAutotrophMatterAppropriation::RemoveHumanAppropriatedMatter(Environment* LocalEnvironment, double NPPWetMatter,double fracBiomass,unsigned currentTimestep) {
        
        // Factor to convert NPP from units per m2 to units per km2
        double m2Tokm2Conversion = 1000000.0;

        //extraction parameters
        std::string humanNPPScenarioType= Parameters::Get( )->GetHumanNPPScenarioType();
        double humanNPPExtractionScale  = Parameters::Get( )->GetHumanNPPExtractionScale();
        double humanNPPScenarioDuration = Parameters::Get( )->GetHumanNPPScenarioDuration();
        unsigned burninSteps            = Parameters::Get( )->GetBurninSteps();
        unsigned impactSteps            = Parameters::Get( )->GetImpactSteps();
        unsigned recoverySteps          = Parameters::Get( )->GetRecoverySteps();
        
        //return value
        double RemovalRate=0;

        if( humanNPPScenarioType == "hanpp" ) {
            
          if (currentTimestep > burninSteps){

          
            // If HANPP value is missing, then assume zero
            
            // Get the total amount of NPP appropriated by humans from this cell
            double HANPP = LocalEnvironment->TerrestrialHANPP() ;

            if( HANPP == Constants::cMissingValue ) HANPP = 0.0;else HANPP = HANPP * LocalEnvironment->Seasonality()*fracBiomass;

            // Convert gC/m2/month to gC/km2/month
            HANPP *= m2Tokm2Conversion;

            // Multiply by cell area (in km2) to get g/cell/day
            HANPP *= LocalEnvironment->Area();

            // Convert from gC to g dry matter
            double DryMatterAppropriated = HANPP * 2;

            // Convert from g dry matter to g wet matter
            double WetMatterAppropriated = DryMatterAppropriated * 2;
            //Calculate the rate of HANPP offtake
            if (NPPWetMatter <= 0.0){
                    RemovalRate = 0.0;
            } else {
                    RemovalRate = min(1.0, WetMatterAppropriated / NPPWetMatter);
            }
          }

        } else if( humanNPPScenarioType == "constant" ){
            
          if (currentTimestep > burninSteps){

            // If the burn-in period has been completed, then remove the specified constant
            // fraction from the acting autotroph stock
                RemovalRate = humanNPPExtractionScale;
          }
          
        } else if (humanNPPScenarioType == "temporary") {
            // If the spin-up period has been completed and the period of impact has not elapsed,
            // then remove the specified constant fraction from the acting autotroph stock
            if ((currentTimestep > burninSteps) && (currentTimestep <= (burninSteps + impactSteps))){
                RemovalRate = humanNPPExtractionScale;
            }

        } else if (humanNPPScenarioType == "escalating") {
            // If the spin-up period has been completed, then remove a proportion of plant matter
            // according to the number of time-steps that have elapsed since the spin-up ended
            if (currentTimestep > burninSteps)  {
                RemovalRate = (min(1.0, (((currentTimestep - burninSteps) / 12.0) * humanNPPExtractionScale)));
            }

        } else if (humanNPPScenarioType == "temp-escalating") {
            // If the spin-up period has been completed and the period of impact has not elapsed, 
            // then remove a proportion of plant matter
            // according to the number of time-steps that have elapsed since the spin-up ended
            if ((currentTimestep > burninSteps) && (currentTimestep <= (burninSteps + impactSteps))) {
                        RemovalRate = (min(1.0, (((currentTimestep - burninSteps) / 12.0) * humanNPPExtractionScale)));
            }

        } else if (humanNPPScenarioType == "temp-escalating-const-rate")  {
            // If the spin-up period has been completed and the period of impact (specified by the third scenario element
            // has not elapsed, then remove a proportion of plant matter
            // according to the number of time-steps that have elapsed since the spin-up ended

            unsigned ConstImpactSteps = humanNPPScenarioDuration* Constants::cYear ;
            if ((currentTimestep > burninSteps) && (currentTimestep <= (burninSteps + ConstImpactSteps))) {
                RemovalRate = (min(1.0, (((currentTimestep - burninSteps) / 12.0) * humanNPPExtractionScale)));
            }

        } else if (humanNPPScenarioType == "temp-escalating-const-rate-duration") {
            // If the spin-up period has been completed and the period of impact (specified by the third scenario element)
            // has not elapsed, then remove a proportion of plant matter
            // according to the number of time-steps that have elapsed since the spin-up ended

            unsigned ConstImpactSteps = humanNPPScenarioDuration* Constants::cYear ;

            if ((currentTimestep > burninSteps) && (currentTimestep <= (burninSteps + impactSteps)))  {
                RemovalRate = (min(1.0,
                                   min(((ConstImpactSteps / 12.0) * humanNPPExtractionScale),
                                      (((currentTimestep - burninSteps) / 12.0) * humanNPPExtractionScale))));
            }
;
        } else if (humanNPPScenarioType == "temp-escalating-declining") {
            // If the spin-up period has been completed, then apply a level of harvesting
            // according to the number of time-steps that have elapsed since the spin-up ended
           if ((currentTimestep > burninSteps) && (currentTimestep <= (burninSteps + impactSteps)))  {
                RemovalRate = (min(1.0, (((currentTimestep - burninSteps) / 12.0) * humanNPPExtractionScale)));
            } else if ((currentTimestep > (burninSteps + impactSteps)) & (currentTimestep <= (burninSteps + impactSteps + recoverySteps))){
                RemovalRate = (min(1.0, (((burninSteps + impactSteps + recoverySteps - currentTimestep) / 12.0) * humanNPPExtractionScale)));
            }

        } else if( humanNPPScenarioType == "none" ) {
          //removal rate remains zero

        } else {
            cout<<"humanNPPScenarioType "<<humanNPPScenarioType<<endl;
            cout<<"humanNPPExtractionScale "<<humanNPPExtractionScale <<endl;
            cout<<"humanNPPScenarioDuration "<<humanNPPScenarioDuration <<endl;
            cout<<"burninSteps "<<burninSteps <<endl;           
            cout<<"impactSteps "<<impactSteps <<endl;           
            cout<<"recoverySteps "<<recoverySteps <<endl;          
            assert(humanNPPScenarioType =="none"  );
        }

        return RemovalRate;
    }
    //----------------------------------------------------------------------------------------------

