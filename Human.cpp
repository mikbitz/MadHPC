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
 * Human.cpp
 *
 *  Created on: Sep 1, 2010
 *      Author: nick
 */

#include "relogo/AgentSet.h"
#include "relogo/Patch.h"
#include "MadObserver.h"

#include "repast_hpc/RepastProcess.h"
#include "repast_hpc/Random.h"

#include "Human.h"
#include "Stock.h"
#include "Environment.h"

using namespace repast::relogo;
using namespace repast;

struct CountStocksOnPatch {
	double operator()(const Patch* patch) const {
		AgentSet<Stock> set;
		patch->turtlesOn(set);
		return set.size();
	}
};

//void Human::infect() {
//	if (!_infected) {
//		_infected = true;
		//((MadObserver*) _observer)->incrementInfectionCount();
//	}
//}

void Human::step() {
    _CurrentTimeStep=RepastProcess :: instance ()->getScheduleRunner ().currentTick ();
	// if human is now dead we can't move it because
	// it will be removed from the sim and the synchronization
	// mechanism cannot move it.
	bool alive = true;

	//if (_infected) {
	//	_infectionTime++;
	//	if (_infectionTime == 50) {
	//		// should be safe to die here as nothing else
	//		// will need this object, so OK to delete it
	//		_observer->hatch<Stock> (this);
	//		die();
	//		alive = false;
	//	}
	//}
	
    /*assignTimeActive();
    eat();
    metabolize();
    reproduce();
    mort();
    applyEcology();*/

}
void Human::expire(){
   _alive = true;
//    if (_CohortAbundance - Parameters::Get( )->GetExtinctionThreshold( ) <= 0 || _IndividualBodyMass <= 0){ 
//      die();
//      alive=false;
//    }
}
void Human::moveIt(){
        if (!_alive)return;
        //AgentSet<Patch> nghs = patchHere<Patch> ()->neighbors<Patch> ();
	    //Patch* winningPatch = nghs.minOneOf(CountStocksOnPatch());
		//face(winningPatch);
        //heading(5);
		//double distanceToMove = 1.5; // For non-toroidal worlds, need to check to make sure move is not out of bounds
		//while((_observer->patchAtOffset(location(), heading(), distanceToMove) == 0) && (distanceToMove > 0)) distanceToMove--;
      	//if(distanceToMove > 0) move(distanceToMove);
        // Calculate the scalar to convert from the time step units used by this implementation of dispersal to the global model time step units
        double DeltaT = Constants::cMonth;
        double latCellLength = patchHere<Environment>()->_Height;
        double lonCellLength = patchHere<Environment>()->_Width;
        double CellArea      = patchHere<Environment>()->_Area;

        double DispersalSpeedBodyMassScalar = 0.0278;
        double DispersalSpeedBodyMassExponent = 0.48;

        
        double dispersalSpeed=0.;

        if( _IsPlanktonic ) {
          // Advective dispersal
          //dispersalName = "advective";
          double HorizontalDiffusivity = 100;
          double AdvectiveModelTimeStepLengthHours = 18;
          double HorizontalDiffusivityKmSqPerADTimeStep = HorizontalDiffusivity / ( 1000 * 1000 ) * 60 * 60 * AdvectiveModelTimeStepLengthHours;

          // Initialise the advective dispersal temporal scaling to adjust between time steps appropriately
          double AdvectionTimeStepsPerModelTimeStep = Constants::cDay * 24 / AdvectiveModelTimeStepLengthHours;

          // Convert velocity from m/s to km/month. Note that if the _TimeUnitImplementation changes, this will also have to change.
          double VelocityUnitConversion = 60 * 60 * 24 * Constants::cDay * Constants::cMonth  / 1000;
          
          for( int mm = 0; mm < AdvectionTimeStepsPerModelTimeStep; mm++ ) {
            // Get the u speed and the v speed from the cell data
            double uAdvectiveSpeed = patchHere<Environment>()->_uVel;
            assert( uAdvectiveSpeed > -9999 );

            double vAdvectiveSpeed = patchHere<Environment>()->_vVel;
            assert( vAdvectiveSpeed > -9999 );

            // Note that this formulation drops the delta t because we set the horizontal diffusivity to be at the same temporal
            // scale as the time step
            NormalGenerator NJ= repast::Random::instance()->createNormalGenerator(0,1);
            // Calculate the distance travelled in this dispersal (not global) time step. both advective and diffusive speeds need to have been converted to km / advective model time step
            double uSpeed = uAdvectiveSpeed * VelocityUnitConversion / AdvectionTimeStepsPerModelTimeStep + NJ.next() * sqrt( ( 2.0 * HorizontalDiffusivityKmSqPerADTimeStep ) );
            double vSpeed = vAdvectiveSpeed * VelocityUnitConversion / AdvectionTimeStepsPerModelTimeStep + NJ.next() * sqrt( ( 2.0 * HorizontalDiffusivityKmSqPerADTimeStep ) );

             TryToDisperse( uSpeed,vSpeed );
          }
        }// Otherwise, if mature do responsive dispersal
        else if( _IsMature ) {
          //dispersalName = "responsive";
          double DensityThresholdScaling = 50000;
          double StarvationDispersalBodyMassThreshold = 0.8;
          dispersalSpeed=DispersalSpeedBodyMassScalar * pow( _AdultMass, DispersalSpeedBodyMassExponent);

          // Check for starvation-driven dispersal
          // A boolean to check whether a cohort has dispersed
          bool cohortHasDispersed = false;

          // Check for starvation driven dispersal. Note that at present we are just tracking starvation for adults
          // Assume a linear relationship between probability of dispersal and body mass loss, up to _StarvationDispersalBodyMassThreshold at which point the cohort will try to disperse every time step
          if( _IndividualBodyMass < _AdultMass ) {

             // If the body mass loss is greater than the starvation dispersal body mass threshold, then the cohort tries to disperse
             if( _IndividualBodyMass / _AdultMass < StarvationDispersalBodyMassThreshold ) {
                // Cohort tries to disperse
                TryToDisperse( dispersalSpeed );
                // Note that regardless of whether or not it succeeds,  it is counted as having dispersed for the purposes of not then allowing it to disperse based on its density.
                cohortHasDispersed = true;
                // Otherwise, the cohort has a chance of trying to disperse proportional to its mass lass
             } else {
               // Cohort tries to disperse with a particular probability
               if( ( ( 1.0 - _IndividualBodyMass / _AdultMass ) / ( 1.0 - StarvationDispersalBodyMassThreshold ) ) > repast::Random::instance()->nextDouble() ) {
                 TryToDisperse( dispersalSpeed );
                 cohortHasDispersed = true;
               }
             }
          }
          //not starvation - density driven?
          if( !cohortHasDispersed ) {
            // If below the density threshold
            if( ( _CohortAbundance / CellArea ) < DensityThresholdScaling / _AdultMass ) {
                TryToDisperse( dispersalSpeed );
           }
          }
        }// If the cohort is immature, run diffusive dispersal
        else {
           dispersalSpeed=DispersalSpeedBodyMassScalar * pow( _IndividualBodyMass, DispersalSpeedBodyMassExponent);

           TryToDisperse( dispersalSpeed );
        }
}
 void Human::TryToDisperse(double dispersalSpeed){
    double randomDirection = repast::Random::instance()->nextDouble()* 2 * acos( -1. );

    // Calculate the u and v components given the dispersal speed
    double uSpeed = dispersalSpeed * cos( randomDirection );
    double vSpeed = dispersalSpeed * sin( randomDirection );
    TryToDisperse(uSpeed, vSpeed);
 }
 void Human::TryToDisperse(double uSpeed, double vSpeed){
 // Pick a direction at random
    double latCellLength = patchHere<Environment>()->_Height;
    double lonCellLength = patchHere<Environment>()->_Width;
    double CellArea      = patchHere<Environment>()->_Area;

  
    double AreaOutsideBoth = abs( uSpeed * vSpeed );

    // Calculate the area of the grid cell that is now outside in the u direction (not including the diagonal)
    double AreaOutsideU = abs( uSpeed * latCellLength ) - AreaOutsideBoth;

    // Calculate the proportion of the grid cell that is outside in the v direction (not including the diagonal
    double AreaOutsideV = abs( vSpeed * lonCellLength ) - AreaOutsideBoth;

    // Convert areas to a probability
    double DispersalProbability = ( AreaOutsideU + AreaOutsideV + AreaOutsideBoth ) / CellArea;
    // Check that we don't have any issues
    if( DispersalProbability > 1 ){
       DispersalProbability=1.;
    }
    // Check to see in which axis the cohort disperses

   // Note that the values in the dispersal array are the proportional area moved outside the grid cell in each direction; we simply compare the random draw to this
   // to determine the direction in which the cohort moves probabilistically
   double RandomValue=repast::Random::instance()->nextDouble();
   if( DispersalProbability >= RandomValue ) {
      int signu = ( uSpeed > 0 ) - ( uSpeed < 0 );
      int signv = ( vSpeed > 0 ) - ( vSpeed < 0 );
      // Longitudinally
      if( RandomValue <= AreaOutsideU / CellArea ) {
       signv = 0;
      } else {
       //Latitudinally
       if( RandomValue <= ( AreaOutsideU / CellArea + AreaOutsideV / CellArea ) ) {
           signu = 0;
       }
     }
     facexy(signu,signv);
     if (patchRightAndAhead<Environment>(0,1)->_Marine==patchHere<Environment>()->_Marine)
        moveTo(patchRightAndAhead<Environment>(0,1));
   }
     
}
void Human::assignTimeActive(){

    double WarmingTolerance = 0;
    double ThermalSafetyMargin = 0;
    double Topt = 0;
    double CTmax = 0;
    double CTmin = 0;
    double AmbientTemp = 0;
    double DTR = 0;
    // Initialise ecological parameters for predation
    // Source: Deutsch et al (2008), Impacts of climate warming on terrestrial ectotherms across latitude, PNAS.
    double TerrestrialWarmingToleranceIntercept = 6.61;
    double TerrestrialWarmingToleranceSlope = 1.6;
    double TerrestrialTSMIntercept = 1.51;
    double TerrestrialTSMSlope = 1.53;
    double Pi = acos(-1);

    //Only work on heterotroph cohorts
    if( _Heterotroph ) {
        //Check if this is an endotherm or ectotherm
        if( _Endotherm ) {
            //Assumes the whole timestep is suitable for endotherms to be active - actual time active is therefore the proportion specified for this functional group.
            _ProportionTimeActive = _ProportionSuitableTimeActive;;
        } else {
            //If ectotherm then use realm specific function
            if( _Marine){
                 double ProportionTimeSuitableMarine = 1.0;
                _ProportionTimeActive = ProportionTimeSuitableMarine * _ProportionSuitableTimeActive;
            } else {
               AmbientTemp = patchHere<Environment> ()->_Temperature;
               DTR         = patchHere<Environment> ()->_DiurnalTemperatureRange;

               //Calculate the Warming tolerance and thermal safety margin given standard deviation of monthly temperature
               WarmingTolerance =    TerrestrialWarmingToleranceSlope * patchHere<Environment> ()->_SDTemperature + TerrestrialWarmingToleranceIntercept;
               ThermalSafetyMargin = TerrestrialTSMSlope *              patchHere<Environment> ()->_SDTemperature + TerrestrialTSMIntercept;

               Topt  = ThermalSafetyMargin + patchHere<Environment> ()->_AnnualTemperature;
               CTmax = WarmingTolerance    + patchHere<Environment> ()->_AnnualTemperature;

               double PerformanceStandardDeviation = ( CTmax - Topt ) / 12;

               CTmin = Topt - 4 * PerformanceStandardDeviation;

               //Calculate the diurnal maximum in the current month
               double DTmax = AmbientTemp + ( 0.5 * DTR );
               double DTmin = AmbientTemp - ( 0.5 * DTR );

               double temp = 2 * ( CTmax - AmbientTemp ) / DTR;
               if( CTmax - DTmax > 0.0 ) {
                  temp = 1.0;
               } else if( CTmax - DTmin < 0.0 ) {
                  temp = -1.0;
               }

               //Proportion of time for which ambient temperatures are greater than the critical upper temperature
               double POver = ( ( Pi / 2.0 ) - asin( temp ) ) / Pi;

               if( CTmin - DTmax > 0.0 ) {
                  temp = 1.0;
               } else if( CTmin - DTmin < 0.0 ) {
                  temp = -1.0;
               } else {
                  temp = 2 * ( CTmin - AmbientTemp ) / DTR;
               }
               //Proportion of time for which ambient temperatures are below the critical lower temperature

               double PBelow = 1 - ( ( Pi / 2.0 ) - asin( temp ) ) / Pi;

               double  ProportionDaySuitable = 1 - ( POver + PBelow );
               double ProportionTimeSuitableTerrestrial= ProportionDaySuitable;
               _ProportionTimeActive = ProportionTimeSuitableTerrestrial * _ProportionSuitableTimeActive;

            }
        }
    }
}
void Human::metabolize(){
    if (_Ectotherm){
    
    // Parameters from fitting to Nagy 1999 Field Metabolic Rates for reptiles - assumes that reptile FMR was measured with animals at their optimal temp of 30degC
    double MetabolismMassExponent = 0.88;
    double NormalizationConstant = 148984000000; //1.4898373851E+11;
    double ActivationEnergy = 0.69; // includes endotherms in hibernation and torpor
    double BoltzmannConstant = 8.617e-5;
    // BMR normalisation constant from Brown et al 2004 - original units of J/s so scale to kJ/d
    double NormalizationConstantBMR = 41918272883; //exp( 20 )*60 * 60 * 24 / 1000;
    double BasalMetabolismMassExponent = 0.69;
    // Currently a very rough estimate based on calorific values of fat, protein and carbohydrate - assumes organism is metabolising mass of 1/4 protein, 1/4 carbohydrate and 1/2 fat 
    double EnergyScalar = 0.036697248; //1 / 27.25;
    // Set the constant to convert temperature in degrees Celsius to Kelvin
    double TemperatureUnitsConvert = 273.0;
    // Calculate the scalar to convert from the time step units used by this implementation of metabolism to the global  model time step units
    double DeltaT = Constants::cDay;
    double temperature=patchHere<Environment> ()->_Temperature + TemperatureUnitsConvert;
    bool ProportionTimeActiveCalculatedThisTimestep = false;
    double FieldMetabolicLosskJ = NormalizationConstant    * pow( _IndividualBodyMass, MetabolismMassExponent )      * exp( -( ActivationEnergy / ( BoltzmannConstant * temperature ) ) );
    double BasalMetabolicLosskJ = NormalizationConstantBMR * pow( _IndividualBodyMass, BasalMetabolismMassExponent ) * exp( -( ActivationEnergy / ( BoltzmannConstant * temperature ) ) );

    // metabolic loss in grams
    double IndividualMetabolicRate = (( _ProportionTimeActive * FieldMetabolicLosskJ ) + ( ( 1 - _ProportionTimeActive ) * ( BasalMetabolicLosskJ ) ) ) * EnergyScalar;
    _MassAccounting["biomass"]["metabolism"] = -IndividualMetabolicRate * DeltaT;

    // If metabolic loss is greater than individual body mass after herbivory and predation, then set equal to individual body mass
    _MassAccounting["biomass"]["metabolism"] = std::max( _MassAccounting["biomass"]["metabolism"], -( _IndividualBodyMass + _MassAccounting["biomass"]["predation"] + _MassAccounting["biomass"]["herbivory"] ) );

    // Add total metabolic loss for all individuals in the cohort to delta biomass for metabolism in the respiratory CO2 pool
    _MassAccounting["respiratoryCO2pool"]["metabolism"] = -_MassAccounting["biomass"]["metabolism"] * _CohortAbundance;

    }
    if (_Endotherm){

    // Parameters from fitting to Nagy 1999 Field Metabolic Rates for mammals and birds, and assuming that these endotherms are metabolising with a body temperature of 310K (37C)
    double MetabolismMassExponent = 0.7;
    double NormalizationConstant = 9.0809083973E+11;
    double ActivationEnergy = 0.69; // includes endotherms in hibernation and torpor
    double BoltzmannConstant = 8.617e-5;
    // Currently a very rough estimate based on calorific values of fat, protein and carbohydrate - assumes organism is metabolising mass of 1/4 protein, 1/4 carbohydrate and 1/2 fat 
    double EnergyScalar = 1 / 27.25;
    // Set the constant to convert temperature in degrees Celsius to Kelvin
    double TemperatureUnitsConvert = 273.0;
    // Assume all endotherms have a constant body temperature of 37degC
    double EndothermBodyTemperature = 37.0 + TemperatureUnitsConvert;
    // Calculate the scalar to convert from the time step units used by this implementation of metabolism to the global  model time step units
    double DeltaT = Constants::cDay;
    // Calculate metabolic loss in kJ
    double metabolicLosskJ = NormalizationConstant * pow( _IndividualBodyMass, MetabolismMassExponent ) * exp( -( ActivationEnergy / ( BoltzmannConstant * EndothermBodyTemperature ) ) );
    // metabolic loss in grams
    double IndividualMetabolicRate= metabolicLosskJ * EnergyScalar;
    // Calculate metabolic loss for an individual and add the value to the delta biomass for metabolism
    _MassAccounting[ "biomass" ][ "metabolism" ] = -IndividualMetabolicRate * DeltaT;

    // If metabolic loss is greater than individual body mass after herbivory and predation, then set equal to individual body mass
    _MassAccounting[ "biomass" ][ "metabolism" ] = std::max( _MassAccounting[ "biomass" ][ "metabolism" ], -( _IndividualBodyMass + _MassAccounting[ "biomass" ][ "predation" ] + _MassAccounting[ "biomass" ][ "herbivory" ] ) );

    // Add total metabolic loss for all individuals in the cohort to delta biomass for metabolism in the respiratory CO2 pool
    _MassAccounting[ "respiratoryCO2pool" ][ "metabolism" ] = -_MassAccounting[ "biomass" ][ "metabolism" ] * _CohortAbundance;
    }
    //heterotroph - in the original model but not used */
}
void Human::reproduce(){
        
    // Biomass per individual in each cohort to be assigned to reproductive potential
    double BiomassToAssignToReproductivePotential;

    // Net biomass change from other ecological functions this time step
    double NetBiomassFromOtherEcologicalFunctionsThisTimeStep;

    // Reset variable holding net biomass change of individuals in this cohort as a result of other ecological processes
    NetBiomassFromOtherEcologicalFunctionsThisTimeStep = 0.0;

    // Loop over all items in the biomass deltas
    for( auto Biomass: _MassAccounting[ "biomass" ] ) {
        // Add the delta biomass to net biomass
        NetBiomassFromOtherEcologicalFunctionsThisTimeStep += Biomass.second;
    }

    // If individual body mass after the addition of the net biomass from processes this time step will yield a body mass 
    // greater than the adult body mass for this cohort, then assign the surplus to reproductive potential
    if( ( _IndividualBodyMass + NetBiomassFromOtherEcologicalFunctionsThisTimeStep ) > _AdultMass ) {
        // Calculate the biomass for each individual in this cohort to be assigned to reproductive potential
        BiomassToAssignToReproductivePotential = _IndividualBodyMass + NetBiomassFromOtherEcologicalFunctionsThisTimeStep - _AdultMass;
        // Check that a positive biomass is to be assigned to reproductive potential
        assert( BiomassToAssignToReproductivePotential >= 0.0 && "Assignment of negative reproductive potential mass" );

        // If this is the first time reproductive potential mass has been assigned for this cohort, 
        // then set the maturity time step for this cohort as the current model time step
        if( _MaturityTimeStep == std::numeric_limits< unsigned >::max( ) ) {
            _MaturityTimeStep = _CurrentTimeStep;
        }

        // Assign the specified mass to reproductive potential mass and remove it from individual biomass
        _MassAccounting[ "reproductivebiomass" ][ "reproduction" ] += BiomassToAssignToReproductivePotential;
        _MassAccounting[ "biomass" ][ "reproduction" ] -= BiomassToAssignToReproductivePotential;

    } else {
        // Cohort has not gained sufficient biomass to assign any to reproductive potential, so take no action
    }

    // Run reproductive events. Note that we can't skip juveniles here as they could conceivably grow to adulthood and get enough biomass to reproduce in a single time step
    // due to other ecological processes

    std:: string TimeUnitImplementation = "month";
    double MassRatioThreshold = 1.5;
    double MassEvolutionProbabilityThreshold = 0.95;
    double MassEvolutionStandardDeviation = 0.05;
    double SemelparityAdultMassAllocation = 0.5;

    // Calculate the scalar to convert from the time step units used by this implementation of dispersal to the global model time step units
    double DeltaT = Constants::cMonth;

    // Adult non-reproductive biomass lost by semelparous organisms
    double adultMassLost;

    // Offspring cohort abundance
    double offspringCohortAbundance;

    // Mass ratio of body mass + reproductive mass to adult body mass
    double currentMassRatio;

    // Individual body mass including change this time step as a result of other ecological processes
    double bodyMassIncludingChangeThisTimeStep;

    // Individual reproductive mass including change this time step as a result of other ecological processes
    double reproductiveMassIncludingChangeThisTimeStep;

    // Calculate the biomass of an individual in this cohort including changes this time step from other ecological processes  
    bodyMassIncludingChangeThisTimeStep = 0.0;

    for( auto& Biomass: _MassAccounting[ "biomass" ] ) {
        // Add the delta biomass to net biomass
        bodyMassIncludingChangeThisTimeStep += Biomass.second;
    }
    bodyMassIncludingChangeThisTimeStep += _IndividualBodyMass;

    // Calculate the reproductive biomass of an individual in this cohort including changes this time step from other ecological processes  
    reproductiveMassIncludingChangeThisTimeStep = 0.0;

    for( auto& ReproBiomass: _MassAccounting[ "reproductivebiomass" ] ) {
        // Add the delta reproductive biomass to net biomass
        reproductiveMassIncludingChangeThisTimeStep += ReproBiomass.second;
    }
    reproductiveMassIncludingChangeThisTimeStep += _IndividualReproductivePotentialMass;
    if ( _IndividualBodyMass > 1.e-200 ) {
        // Get the current ratio of total individual mass (including reproductive potential) to adult body mass
        currentMassRatio = ( bodyMassIncludingChangeThisTimeStep + reproductiveMassIncludingChangeThisTimeStep ) / _AdultMass;

        // Must have enough mass to hit reproduction threshold criterion, and either (1) be in breeding season, or (2) be a marine cell (no breeding season in marine cells)
        if ( ( currentMassRatio > MassRatioThreshold ) && (patchHere<Environment> ()->_Breeding_Season == 1.0  ||  _Marine  ) ) {
            // Iteroparous and semelparous organisms have different strategies
            if ( _Iteroparous ) {
                // Iteroparous organisms do not allocate any of their current non-reproductive biomass to reproduction
                adultMassLost = 0.0;

                // Calculate the number of offspring that could be produced given the reproductive potential mass of individuals
                offspringCohortAbundance = _CohortAbundance * reproductiveMassIncludingChangeThisTimeStep / _JuvenileMass;
            } else {
                // Semelparous organisms allocate a proportion of their current non-reproductive biomass (including the effects of other ecological processes) to reproduction
                adultMassLost = SemelparityAdultMassAllocation * bodyMassIncludingChangeThisTimeStep;

                // Calculate the number of offspring that could be produced given the reproductive potential mass of individuals
                offspringCohortAbundance = ( ( _CohortAbundance ) * ( adultMassLost + reproductiveMassIncludingChangeThisTimeStep ) ) / _JuvenileMass;
            }
            // Check that the abundance in the cohort to produce is greater than or equal to zero
            assert( offspringCohortAbundance >= 0.0 && "Offspring abundance < 0" );

            // Get the adult and juvenile masses of the offspring cohort

            // Determine whether offspring cohort 'evolves' in terms of adult and juvenile body masses
            // If not, it just gets the same values as the parent cohort
            double RandomValue = repast::Random::instance()->nextDouble();//mRandomNumber.GetUniform( );
            double newJuvenileMass = _JuvenileMass;
            double newAdultMass    = _AdultMass;
            if( RandomValue > MassEvolutionProbabilityThreshold ) {
             // Determine the new juvenile body mass
             NormalGenerator NJ= repast::Random::instance()->createNormalGenerator(_JuvenileMass,MassEvolutionStandardDeviation * _JuvenileMass);
             double RandomValueJ = NJ.next();
             newJuvenileMass = std::max( RandomValueJ, _MinimumMass );

             // Determine the new adult body mass
             NormalGenerator NA= repast::Random::instance()->createNormalGenerator(_AdultMass,MassEvolutionStandardDeviation * _AdultMass);

             double RandomValueA = NA.next();
             newAdultMass = std::min( RandomValueA, _MaximumMass );

             } 
             // Update cohort abundance in case juvenile mass has been altered through 'evolution'
             offspringCohortAbundance = offspringCohortAbundance * ( _JuvenileMass / newJuvenileMass );


            Human* newH = _observer->hatch<Human> (this);
            newH->setOffspring( this, newJuvenileMass, newAdultMass, newJuvenileMass, offspringCohortAbundance, _CurrentTimeStep);

            // Subtract all of the reproductive potential mass of the parent cohort, which has been used to generate the new
            // cohort, from the delta reproductive potential mass and delta adult body mass
            _MassAccounting[ "reproductivebiomass" ][ "reproduction" ] -= reproductiveMassIncludingChangeThisTimeStep;
            _MassAccounting[ "biomass" ][ "reproduction" ] -= adultMassLost;
        } else {
            // Organism is not large enough, or it is not the breeding season, so take no action
        }
    }
}
void Human::setOffspring( Human* actingCohort, double juvenileBodyMass, double adultBodyMass, double initialBodyMass, double initialAbundance, unsigned birthTimeStep ) {
    _FunctionalGroupIndex        = actingCohort->_FunctionalGroupIndex;
    _JuvenileMass                = juvenileBodyMass;
    _AdultMass                   = adultBodyMass;
    _IndividualBodyMass          = initialBodyMass;
    _CohortAbundance             = initialAbundance;
    _BirthTimeStep               = birthTimeStep;
    _MaturityTimeStep            = std::numeric_limits<unsigned>::max( );
    _LogOptimalPreyBodySizeRatio = actingCohort->_LogOptimalPreyBodySizeRatio;
    _MaximumAchievedBodyMass     = juvenileBodyMass;
    _Merged                      = false;
    _ProportionTimeActive        = actingCohort->_ProportionTimeActive;
    _IndividualReproductivePotentialMass = 0;

}
void Human::eat(){
    
    double  PotentialBiomassEaten=0.,HandlingTime=0.,HandlingTimeScaled=0.,BiomassesEaten=0.,IndividualHerbivoryRate=0.;
    double edibleFraction=0,AttackRateExponent=0, HandlingTimeExponent=0,HandlingTimeScalar=0;
    double _SearchRateConstant=0,_DaysInATimeStep=Constants::cDay;
    double PredatorAbundanceMultipliedByTimeEating=0;
    double _FeedingPreferenceStandardDeviation=0.;
    double CellAreaHectares=patchHere<Environment> ()->_Area/100;
    unsigned _NumberOfBins=12;
    AgentSet<Human> preys= turtlesHere<Human>();
    AgentSet<Stock> stocks = turtlesHere<Stock>();

    if (_Herbivore || _Omnivore){

       //acumulate handling time - these values will be for herbies - note the issue if these are stored in the class with omnis
       if (_Marine){
       edibleFraction       = _edibleFractionMarine;
       AttackRateExponent   = _AttackRateExponentMarine;
       HandlingTimeExponent = _HandlingTimeExponentMarine;
       HandlingTimeScalar   = _HandlingTimeScalarMarine;
       }
       if(_Terrestrial){
       edibleFraction       = _edibleFractionTerrestrial;
       AttackRateExponent   = _AttackRateExponentTerrestrial;
       HandlingTimeExponent = _HandlingTimeExponentTerrestrial;
       HandlingTimeScalar   = _HandlingTimeScalarTerrestrial;
       }
       IndividualHerbivoryRate               = _HerbivoryRateConstant    * pow( _IndividualBodyMass, ( _HerbivoryRateMassExponent ) );
       for (auto stock: stocks){
        if( stock->_TotalBiomass > 0.0 ) {
          PotentialBiomassEaten  = IndividualHerbivoryRate   * pow( (stock->_TotalBiomass*edibleFraction) / CellAreaHectares, AttackRateExponent ); // could store this to avoid second calc.?
          HandlingTimeScaled     = HandlingTimeScalar        * pow( ( _ReferenceMass / _IndividualBodyMass ), HandlingTimeExponent );
          HandlingTime          += PotentialBiomassEaten* HandlingTime;
        }        
       }

    }
    std::map<int, vector<double> > BinnedPreyDensities;
    if (_Carnivore || _Omnivore){
         if( _CohortAbundance > 0 ) {

        // Pre-calculate individual values for this predator to speed things up
        double _MaximumSearchRate = _SearchRateConstant * _IndividualBodyMass; // Ref: Predator-only component of eq. 35 from text S1.
        double _DaysEating = _DaysInATimeStep * _ProportionTimeActive;

        double _LogOptimalPreySize = _LogOptimalPreyBodySizeRatio;

        // If a filter feeder, then optimal body size is a value not a ratio: convert it to a ratio to ensure that all calculations work correctly
        if( _IsFilterFeeder ) {
            // Optimal body size is actually a value, not a ratio, so convert it to a ratio based on the present body size
            _LogOptimalPreySize = std::log( std::exp( _LogOptimalPreyBodySizeRatio ) / _IndividualBodyMass ); //LogOptimalPreyBodySizeRatio-log(actingCohort.IndividualBodyMass);
        }
        // Calculate the reference mass scaling ratio
        HandlingTimeScaled = HandlingTimeScalar * std::pow( _ReferenceMass / _IndividualBodyMass, HandlingTimeExponent );
        PredatorAbundanceMultipliedByTimeEating = _CohortAbundance * _DaysEating;

        // Calculate the abundance of prey in each of the prey mass bins

        // Loop through prey functional groups


        for (auto prey:preys){
            // Calculate the difference between the actual body size ratio and the optimal ratio, 
            // and then divide by the standard deviation in log ratio space to determine in 
            // which bin to assign the prey item.
            if (BinnedPreyDensities.count(prey->_FunctionalGroupIndex)==0){BinnedPreyDensities[prey->_FunctionalGroupIndex].resize(_NumberOfBins);
                  for( unsigned binIndex = 0; binIndex < _NumberOfBins; binIndex++ ) BinnedPreyDensities[ prey->_FunctionalGroupIndex][ binIndex ] = 0;}
            if( prey->_IndividualBodyMass > 0 ) {

                unsigned binIndex = ( unsigned )( ( ( std::log( prey->_IndividualBodyMass / _IndividualBodyMass ) - _LogOptimalPreySize ) / ( 0.5 * _FeedingPreferenceStandardDeviation ) ) + ( _NumberOfBins / 2 ) );

                if( ( 0 < binIndex ) && ( binIndex < _NumberOfBins ) ) {
                    BinnedPreyDensities[ prey->_FunctionalGroupIndex ][ binIndex ] += prey->_CohortAbundance / CellAreaHectares;
                }
            }
        }
        
        // Loop over potential prey functional groups
        for (auto prey: preys){
                prey->_PotentialAbundanceEaten = 0;
                //No Cannibalism
                if( prey->getId() != getId() &&  prey->_IndividualBodyMass > 0) {
                   unsigned binIndex = ( unsigned )( ( ( std::log( prey->_IndividualBodyMass / _IndividualBodyMass ) - _LogOptimalPreySize ) / ( 0.5 * _FeedingPreferenceStandardDeviation ) ) + ( _NumberOfBins / 2 ) );
                    if( (_IsFilterFeeder &&  prey->_IsPlanktonic ) || !_IsFilterFeeder) {

                        if(   ( 0 < binIndex ) && ( binIndex < _NumberOfBins ) ) {
                            // Calculate the potential abundance from this cohort eaten by the acting cohort
                            
                                // Calculate the relative feeding preference from a log-normal distribution with mean equal to the optimal prey to predator ratio and standard deviation as specified
    double preferenceForPrey = std::exp( -( std::pow( ( ( std::log( prey->_IndividualBodyMass / _IndividualBodyMass ) - _LogOptimalPreySize ) / _FeedingPreferenceStandardDeviation ), 2 ) ) ); // Ref: Text S1, eq 36.
    // Calculate the killing rate of an individual predator per unit prey density per hectare per time unit
    double searchRate = _MaximumSearchRate * preferenceForPrey; // Ref: Text S1, eq 35.
    // Calculate the potential number of prey killed given the number of prey detections
    double potentialAbundanceEaten = searchRate * ( prey->_CohortAbundance / CellAreaHectares ) * BinnedPreyDensities[ prey->_FunctionalGroupIndex ][ binIndex ]; // Ref: Text S1, eq. 34.
                            
                            
                            prey->_PotentialAbundanceEaten = potentialAbundanceEaten;

                            // Add the time required to handle the potential abundance eaten from this cohort to the cumulative total for all cohorts
                            HandlingTime += potentialAbundanceEaten * ( HandlingTimeScaled * prey->_IndividualBodyMass );
                       
                    } 
                }
            }
        }
    }
    }
    if (_Herbivore || _Omnivore){
       //accumulate eaten
       for (auto stock: stocks){
        double InstantFractionEaten = 0.0;
        double EdibleMass=0;
        if( stock->_TotalBiomass > 0.0 ) {
          PotentialBiomassEaten               = IndividualHerbivoryRate   * pow( (stock->_TotalBiomass*edibleFraction) / CellAreaHectares, AttackRateExponent );
          EdibleMass                          = stock->_TotalBiomass*edibleFraction;
          InstantFractionEaten                = _CohortAbundance * ( ( PotentialBiomassEaten / ( 1 + HandlingTime ) ) / EdibleMass );
          BiomassesEaten                      = EdibleMass * ( 1 - exp( -InstantFractionEaten * Constants::cDay*_ProportionTimeActive ) ); // should be min(stock._TotalBiomass,...)
          stock->_TotalBiomass               -= BiomassesEaten;
        } 
        // Return the total  biomass of the autotroph stock eaten
       
        if( _CohortAbundance > 0 )_MassAccounting["biomass"]["herbivory"] += BiomassesEaten * _AssimilationEfficiency / _CohortAbundance;
        _MassAccounting["organicpool"]["herbivory"] += BiomassesEaten * ( 1 - _AssimilationEfficiency );
       }
    }
    if (_Carnivore || _Omnivore){
    if( _CohortAbundance > 0 ) {
        double biomassConsumed = 0;

        _AssimilationEfficiency = 0;//FunctionalGroups::Get( )->GetCohorts( ).GetBiologicalPropertyOneFunctionalGroup( "carnivory assimilation", predator->mFunctionalGroupIndex );

        // Loop over potential prey functional groups
        for (auto prey: preys){
                // Calculate the actual abundance of prey eaten from this cohort
                double _AbundancesEaten = 0;
                if( prey->_CohortAbundance > 0 ) {
                    // Calculate the actual abundance of prey eaten from this cohort
                    _AbundancesEaten = prey->_CohortAbundance * ( 1 - std::exp( -( PredatorAbundanceMultipliedByTimeEating * ( ( prey->_PotentialAbundanceEaten / ( HandlingTime + 1 ) ) / prey->_CohortAbundance ) ) ) );
                }
                // Remove number of prey eaten from the prey cohort
                prey->_CohortAbundance -= _AbundancesEaten;

                biomassConsumed += ( prey->_IndividualBodyMass + prey->_IndividualReproductivePotentialMass ) * _AbundancesEaten / _CohortAbundance; //per capita
        }
        // Add the biomass eaten and assimilated by an individual to the delta biomass for the acting (predator) cohort
        _MassAccounting[ "biomass" ][ "carnivory" ] = biomassConsumed * _AssimilationEfficiency;
        // Move the biomass eaten but not assimilated by an individual into the organic matter pool
        _MassAccounting[ "organicpool" ][ "carnivory" ] = biomassConsumed * ( 1 - _AssimilationEfficiency ) * _CohortAbundance;
    }
    }
}


void Human::mort(){

    // Variables to hold the mortality rates
    double MortalityRateBackground;
    double MortalityRateSenescence;
    double MortalityRateStarvation;

    // Variable to hold the total abundance lost to all forms of mortality
    double MortalityTotal;

    // Individual body mass including change this time step as a result of other ecological processes
    double BodyMassIncludingChangeThisTimeStep;

    // Individual reproductive mass including change this time step as a result of other ecological processes
    double ReproductiveMassIncludingChangeThisTimeStep;

    // Calculate the body mass of individuals in this cohort including mass gained through eating this time step, up to but not exceeding adult body mass for this cohort. 
    // Should be fine because these deductions are made in the reproduction implementation, but use Math.Min to double check.

    BodyMassIncludingChangeThisTimeStep = 0.0;
    // Loop over all items in the biomass deltas
    for( auto Biomass: _MassAccounting[ "biomass" ] ) {
        // Add the delta biomass to net biomass
        BodyMassIncludingChangeThisTimeStep += Biomass.second;
    }

    BodyMassIncludingChangeThisTimeStep = std::min( _AdultMass, BodyMassIncludingChangeThisTimeStep + _IndividualBodyMass );

    // Temporary variable to hold net reproductive biomass change of individuals in this cohort as a result of other ecological processes
    ReproductiveMassIncludingChangeThisTimeStep = 0.0;

    // Loop over all items in the biomass Cohort::Deltas
    for( auto Biomass: _MassAccounting[ "reproductivebiomass" ] ) {
        // Add the delta biomass to net biomass
        ReproductiveMassIncludingChangeThisTimeStep += Biomass.second;
    }
    std::string TimeUnitImplementation = "Day";
    // Calculate the scalar to convert from the time step units used by this implementation of mortality to the global model time step units
    double DeltaT = Constants::cDay;
    double MortalityRate;
    
    ReproductiveMassIncludingChangeThisTimeStep += _IndividualReproductivePotentialMass;
    // Check to see if the cohort has already been killed by predation etc
    if( BodyMassIncludingChangeThisTimeStep < 1.e-7 ) {
        //MB a small number ! maybe should be larger? (e.g. min cohort body mass))
        //This causes a difference between C# and c++ versions as there is a rounding error issue - changed in C# code to match
        // If individual body mass is not greater than zero, then all individuals become extinct
        MortalityTotal = 0.0;
        // used to be actingCohort->mCohortAbundance, but this leads to a large cancellation in applyEcology
        //so mortality total is changed here  - now we only multiply by mortality total - so values can never be negative
        //BodyMassIncludingChangeThisTimeStep = 0; //MB would be a kludge to exclude negative values below - need mass checking throughout the code
    } else {
        // Calculate background mortality rate
        MortalityRate = 0.001;
        MortalityRateBackground = MortalityRate * DeltaT;

        // If the cohort has matured, then calculate senescence mortality rate, otherwise set rate to zero
        if( _MaturityTimeStep < std::numeric_limits< unsigned >::max( ) ) {
            
          MortalityRate = 0.003;

          double TimeToMaturity = _MaturityTimeStep - _BirthTimeStep;

          // Calculate how many model time steps since the cohort reached maturity
          double AgePostMaturity = _CurrentTimeStep - _MaturityTimeStep;

          // Calculate the time since maturity as a fraction of the time that it took the cohort to reach maturity
          double FractionalAgePostMaturity = AgePostMaturity / ( TimeToMaturity + 1 );

          // Calculate the mortality rate per mortality formulation time step as a function of the exponential of the previous fraction
          double AgeRelatedMortalityRate = MortalityRate * exp( FractionalAgePostMaturity );

          // Convert the mortality rate from formulation time step units to model time step units
          MortalityRateSenescence= AgeRelatedMortalityRate * DeltaT;

        } else {
            MortalityRateSenescence = 0.0;
        }
        // Calculate the starvation mortality rate based on individual body mass and maximum body mass ever
        // achieved by this cohort
        MortalityRate = 0.;
        if( BodyMassIncludingChangeThisTimeStep < _MaximumAchievedBodyMass ) {
          double LogisticInflectionPoint = 0.6;
          double LogisticScalingParameter = 0.05;
          double MaximumStarvationRate = 1;
          // Calculate the first part of the relationship between body mass and mortality rate
          double k = -( BodyMassIncludingChangeThisTimeStep - LogisticInflectionPoint * _MaximumAchievedBodyMass ) / ( LogisticScalingParameter * _MaximumAchievedBodyMass );
          // Calculate mortality rate
          MortalityRate = MaximumStarvationRate / ( 1 + exp( -k ) );
        }
        MortalityRateStarvation = MortalityRate * DeltaT;

        // Calculate the number of individuals that suffer mortality this time step from all sources of mortality
        MortalityTotal = exp( -MortalityRateBackground - MortalityRateSenescence - MortalityRateStarvation );
    }

    // Remove individuals that have died from the delta abundance for this cohort
    _MassAccounting[ "abundance" ][ "mortality" ] = MortalityTotal;

    // Add the biomass of individuals that have died to the delta biomass in the organic pool (including reproductive 
    // potential mass, and mass gained through eating, and excluding mass lost through metabolism)
    _MassAccounting[ "organicpool" ][ "mortality" ] = ( 1 - MortalityTotal ) * _CohortAbundance * ( BodyMassIncludingChangeThisTimeStep + ReproductiveMassIncludingChangeThisTimeStep );
}
void Human::applyEcology(){

    // Variable to calculate net abundance change to check that cohort abundance will not become negative
    double NetAbundanceChange = 0.0;
    // Loop over all abundance deltas
    for( auto& d: _MassAccounting["abundance"] ) {
        // Update net abundance change
        NetAbundanceChange += d.second;
    }
    // Check that cohort abundance will not become negative
    assert( NetAbundanceChange >= 0 && "Cohort abundance < 0" );
    //this is a variation from the original below - 
    //the abundance here only changes in Mortality.h where the change was mortalityTotal=(1-exp(-sum-of-Mortalities))*abundance
    //this is now change to mortalityTotal=exp(-sum-of-mortalities) - avoiding possibilities of negatives arising.
    _CohortAbundance *= NetAbundanceChange;

    // Variable to calculate net biomass change to check that cohort individual body mass will not become negative
    double NetBiomass = 0.0;

    // Loop over all biomass deltas
    for( auto& d: _MassAccounting["biomass"] ) {
        // Update net biomass change
        NetBiomass += d.second;
    }
    double BiomassCheck = 0.0;
    bool NetToBeApplied = true;
    // If cohort abundance is greater than zero, then check that the calculated net biomass will not make individual body mass become negative
    if( _CohortAbundance > 0 ) {

        BiomassCheck = _IndividualBodyMass + NetBiomass;
        if( BiomassCheck < 0 ) {
            std::cout << "Biomass going negative, acting cohort: " << _FunctionalGroupIndex << ", " << getId() << std::endl;
            exit( 1 );
        }
    }

    //Loop over all keys in the deltas sorted list
    for( auto& d: _MassAccounting["biomass"] ) {
        // If cohort abundance is zero, then set cohort individual body mass to zero and reset the biomass delta to zero, 
        // otherwise update cohort individual body mass and reset the biomass delta to zero
        if( _CohortAbundance == 0 ) {
            _IndividualBodyMass = 0.0;
        } else {
            if( NetToBeApplied ) {
                _IndividualBodyMass = _IndividualBodyMass + NetBiomass;
                NetToBeApplied = false;
            }
        }
    }
    // Check that individual body mass is still greater than zero
    assert( _IndividualBodyMass >= 0 && "biomass < 0" );

    // If the current individual body mass is the largest that has been achieved by this cohort, then update the maximum achieved
    // body mass tracking variable for the cohort
    if( _IndividualBodyMass > _MaximumAchievedBodyMass )
        _MaximumAchievedBodyMass = _IndividualBodyMass;

    // Variable to calculate net reproductive biomass change to check that cohort individual body mass will not become negative
    double NetReproductiveBiomass = 0.0;

    // Loop over all reproductive biomass deltas
    for( auto& d: _MassAccounting["reproductivebiomass"] ) {
        // Update net reproductive biomass change
        NetReproductiveBiomass += d.second;
    }

    //Loop over all keys in the abundance deltas sorted list
    for( auto& d: _MassAccounting["reproductivebiomass"] ) {
        // If cohort abundance is zero, then set cohort reproductive body mass to zero and reset the biomass delta to zero, 
        // otherwise update cohort reproductive body mass and reset the biomass delta to zero
        if( _CohortAbundance == 0 ) {
            _IndividualReproductivePotentialMass = 0.0;
        } else {
            _IndividualReproductivePotentialMass += d.second;
        }
    }
    //Note that maturity time step is set in TReproductionBasic
}
/*
void EcologyApply::UpdatePools( GridCell& gcl ) {
    // Loop over all keys in the organic pool deltas sorted list
    for( auto &D: Cohort::mMassAccounting["organicpool"] ) {
        // Check that the delta value is not negative
        if( D.second < 0 ) std::cout << "organic pool " << D.first << " " << D.second << std::endl;

        //assert(D.second >= 0.0 && "A delta value for the organic pool is negative " );
        // Update the organic pool biomass
        Environment::Get( "Organic Pool", gcl ) += D.second;
        //Reset the delta value to zero
    }
    // Loop over all keys in the respiratory pool deltas sorted list
    for( auto &D: Cohort::mMassAccounting["respiratoryCO2pool"] ) {
        // Check that the delta value is not negative
        assert( D.second >= 0.0 && "A delta value for the respiratory CO2 pool is negative" );
        // Update the respiratory CO2 pool
        Environment::Get( "Respiratory CO2 Pool", gcl ) += D.second;
        // Reset the delta value to zero // FIX - Should there be a declaration here?
    }
}
*/
//}




