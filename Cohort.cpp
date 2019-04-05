/*
 *
 * Cohort.cpp
 *
 *  Created on: May 7, 2018
 *      Author: Mike Bithell
 * 
 *  Derived from Original C# code by
 * Emergent Global Patterns of Ecosystem Structure and Function from a Mechanistic General Ecosystem Model , Michael B. J. Harfoot , Tim Newbold, Derek P. Tittensor,  Stephen Emmott, Jon Hutton, Vassily Lyutsarev, Matthew J. Smith, Jörn P. W. Scharlemann, Drew W. Purves PLOS Biology 12(4) 2014 e1001841 https://doi.org/10.1371/journal.pbio.1001841
 */


#include "repast_hpc/RepastProcess.h"
#include "repast_hpc/Random.h"

#include "Groups.h"
#include "Cohort.h"
#include "Stock.h"
#include "Parameters.h"
#include "model.h"

using namespace repast;
    //shared constants - these are defaults - Cohorts::setParameters will get these from model.props
    //herbivores

    double Cohort::_edibleFractionMarine               =1.0;
    double Cohort::_AttackRateExponentMarine           =2.0;
    double Cohort::_HandlingTimeExponentMarine         =0.7;
    double Cohort::_HandlingTimeScalarMarine           =0.7;
    double Cohort::_edibleFractionTerrestrial          =0.1;
    double Cohort::_AttackRateExponentTerrestrial      =2.0;
    double Cohort::_HandlingTimeExponentTerrestrial    =0.7;
    double Cohort::_HandlingTimeScalarTerrestrial      =0.7;
    double Cohort::_HerbivoryRateMassExponent          =1.0;
    double Cohort::_HerbivoryRateConstant              =1.0e-11;
    double Cohort::_ReferenceMass                      =1.0;
 
    //Carnivores

    double Cohort::_HandlingTimeScalar_C               = 0.5;
    double Cohort::_HandlingTimeExponent_C             = 0.7;
    double Cohort::_SearchRateConstant                 = 1E-6;
    double Cohort::_FeedingPreferenceStandardDeviation = 0.7;
    double Cohort::_NumberOfBins                       = 12;
    
    //dispersal constants
    double Cohort::_DispersalSpeedBodyMassScalar       = 0.0278;
    double Cohort::_DispersalSpeedBodyMassExponent     = 0.48;
    //advective dispersal
    double Cohort::_HorizontalDiffusivity                       = 100;
    double Cohort::_AdvectiveModelTimeStepLengthHours           = 18;
    double Cohort::_HorizontalDiffusivityKmSqPerADTimeStep      = _HorizontalDiffusivity / ( 1000 * 1000 ) * 60 * 60 * _AdvectiveModelTimeStepLengthHours;
    // Initialise the advective dispersal temporal scaling to adjust between time steps appropriately
    double Cohort::_AdvectionTimeStepsPerModelTimeStep          = Constants::cDay * 24 / _AdvectiveModelTimeStepLengthHours;
    // Convert velocity from m/s to km/month. Note that if the _TimeUnitImplementation changes, this will also have to change.
    double Cohort::_VelocityUnitConversion                      = 60 * 60 * 24 * Constants::cDay * Constants::cMonth  / 1000;
    //responsive dispersal
    double Cohort::_DensityThresholdScaling                     = 50000;
    double Cohort::_StarvationDispersalBodyMassThreshold        = 0.8;
    
    // Initialise ecological parameters for ectotherms
    // Source: Deutsch et al (2008), Impacts of climate warming on terrestrial ectotherms across latitude, PNAS.
    double Cohort::_TerrestrialWarmingToleranceIntercept = 6.61;
    double Cohort::_TerrestrialWarmingToleranceSlope     = 1.6;
    double Cohort::_TerrestrialTSMIntercept              = 1.51;
    double Cohort::_TerrestrialTSMSlope                  = 1.53;
    
    //Pi!
    double Cohort::_Pi = acos(-1);
    //cell areas are in sq. km.
    double Cohort::_CellAreaToHectares=100;
    
    //metabolism
    //Ectotherms
    // Parameters from fitting to Nagy 1999 Field Metabolic Rates for reptiles - assumes that reptile FMR was measured with animals at their optimal temp of 30degC
    double Cohort::_MetabolismMassExponentEcto = 0.88;
    double Cohort::_NormalizationConstantEcto = 148984000000; //1.4898373851E+11;
    double Cohort::_ActivationEnergyEcto = 0.69; // includes endotherms in hibernation and torpor
    // BMR normalisation constant from Brown et al 2004 - original units of J/s so scale to kJ/d
    double Cohort::_NormalizationConstantBMR = 41918272883; //exp( 20 )*60 * 60 * 24 / 1000;
    double Cohort::_BasalMetabolismMassExponent = 0.69;
    // Currently a very rough estimate based on calorific values of fat, protein and carbohydrate - assumes organism is metabolising mass of 1/4 protein, 1/4 carbohydrate and 1/2 fat 
    double Cohort::_EnergyScalarEcto = 0.036697248; //1 / 27.25;
    
    //Endotherms
    // Parameters from fitting to Nagy 1999 Field Metabolic Rates for mammals and birds, and assuming that these endotherms are metabolising with a body temperature of 310K (37C)
    double Cohort::_MetabolismMassExponentEndo = 0.7;
    double Cohort::_NormalizationConstantEndo = 9.0809083973E+11;
    double Cohort::_ActivationEnergyEndo = 0.69; // includes endotherms in hibernation and torpor
    double Cohort::_BoltzmannConstant = 8.617e-5;
    // Currently a very rough estimate based on calorific values of fat, protein and carbohydrate - assumes organism is metabolising mass of 1/4 protein, 1/4 carbohydrate and 1/2 fat 
    double Cohort::_EnergyScalarEndo = 1 / 27.25;
    // Set the constant to convert temperature in degrees Celsius to Kelvin
    double Cohort::_TemperatureUnitsConvert = 273.0;

    //reproduction
    
    double Cohort::_MassRatioThreshold = 1.5;
    double Cohort::_MassEvolutionProbabilityThreshold = 0.95;
    double Cohort::_MassEvolutionStandardDeviation = 0.05;
    double Cohort::_SemelparityAdultMassAllocation = 0.5;
    
    //mortality
    double Cohort::_MortalityRateBackground  = 0.001;
    double Cohort::_MortalityRateMature      = 0.003;
    double Cohort::_LogisticInflectionPoint  = 0.6;
    double Cohort::_LogisticScalingParameter = 0.05;
    double Cohort::_MaximumStarvationRate    = 1;
    
    //needs to be class variable so that Cohorts can have unique numbers
    unsigned Cohort::_NextID=0;
    //_Accounting is shared by all cohorts - saves *a lot* of memory
    //only works as these are temporaries used only and completely within each call to step() by a cohort
    //NB do not reset within step() (e.g. by having newly reproduced cohorts call ResetAccounts() in setupOffspring() !)
    std::map < std::string, std::map<std::string,double> > Cohort::_Accounting;
//------------------------------------------------------------------------------------------------------------
void Cohort::setParameters(repast::Properties* props){
    //shared constants - static function to read these from parameter file
    //herbivores

     _edibleFractionMarine               = repast::strToDouble(props->getProperty("CohortParameters.edibleFractionMarine"));
     _AttackRateExponentMarine           = repast::strToDouble(props->getProperty("CohortParameters.AttackRateExponentMarine"));
     _HandlingTimeExponentMarine         = repast::strToDouble(props->getProperty("CohortParameters.HandlingTimeExponentMarine"));
     _HandlingTimeScalarMarine           = repast::strToDouble(props->getProperty("CohortParameters.HandlingTimeScalarMarine"));
     _edibleFractionTerrestrial          = repast::strToDouble(props->getProperty("CohortParameters.edibleFractionTerrestrial"));
     _AttackRateExponentTerrestrial      = repast::strToDouble(props->getProperty("CohortParameters.AttackRateExponentTerrestrial"));
     _HandlingTimeExponentTerrestrial    = repast::strToDouble(props->getProperty("CohortParameters.HandlingTimeExponentTerrestrial"));
     _HandlingTimeScalarTerrestrial      = repast::strToDouble(props->getProperty("CohortParameters.HandlingTimeScalarTerrestrial"));
     _HerbivoryRateMassExponent          = repast::strToDouble(props->getProperty("CohortParameters.HerbivoryRateMassExponent"));
     _HerbivoryRateConstant              = repast::strToDouble(props->getProperty("CohortParameters.HerbivoryRateConstant"));
     _ReferenceMass                      = repast::strToDouble(props->getProperty("CohortParameters.ReferenceMass"));
 
    //Carnivores

     _HandlingTimeScalar_C               = repast::strToDouble(props->getProperty("CohortParameters.HandlingTimeScalar_C"));
     _HandlingTimeExponent_C             = repast::strToDouble(props->getProperty("CohortParameters.HandlingTimeExponent_C"));
     _SearchRateConstant                 = repast::strToDouble(props->getProperty("CohortParameters.SearchRateConstant"));
     _FeedingPreferenceStandardDeviation = repast::strToDouble(props->getProperty("CohortParameters.FeedingPreferenceStandardDeviation"));
     _NumberOfBins                       = repast::strToDouble(props->getProperty("CohortParameters.NumberOfBins"));
    
    //dispersal constants
     _DispersalSpeedBodyMassScalar                = repast::strToDouble(props->getProperty("CohortParameters.DispersalSpeedBodyMassScalar"));
     _DispersalSpeedBodyMassExponent              = repast::strToDouble(props->getProperty("CohortParameters.DispersalSpeedBodyMassExponent"));
     _HorizontalDiffusivity                       = repast::strToDouble(props->getProperty("CohortParameters.HorizontalDiffusivity"));
     _AdvectiveModelTimeStepLengthHours           = repast::strToDouble(props->getProperty("CohortParameters.AdvectiveModelTimeStepLengthHours"));
     _HorizontalDiffusivityKmSqPerADTimeStep      = _HorizontalDiffusivity / ( 1000 * 1000 ) * 60 * 60 * _AdvectiveModelTimeStepLengthHours;
     // Initialise the advective dispersal temporal scaling to adjust between time steps appropriately
     _AdvectionTimeStepsPerModelTimeStep          = Constants::cDay * 24 / _AdvectiveModelTimeStepLengthHours;
     // Convert velocity from m/s to km/month. Note that if the _TimeUnitImplementation changes, this will also have to change.
     _VelocityUnitConversion                      = 60 * 60 * 24 * Constants::cDay * Constants::cMonth  / 1000; 
     _DensityThresholdScaling                     = repast::strToDouble(props->getProperty("CohortParameters.DensityThresholdScaling"));
     _StarvationDispersalBodyMassThreshold        = repast::strToDouble(props->getProperty("CohortParameters.StarvationDispersalBodyMassThreshold"));

     // Initialise ecological parameters for ectotherms
     _TerrestrialWarmingToleranceIntercept        = repast::strToDouble(props->getProperty("CohortParameters.TerrestrialWarmingToleranceIntercept"));
     _TerrestrialWarmingToleranceSlope            = repast::strToDouble(props->getProperty("CohortParameters.TerrestrialWarmingToleranceSlope"));
     _TerrestrialTSMIntercept                     = repast::strToDouble(props->getProperty("CohortParameters.TerrestrialTSMIntercept"));
     _TerrestrialTSMSlope                         = repast::strToDouble(props->getProperty("CohortParameters.TerrestrialTSMSlope"));
    
     _Pi                                          = acos(-1.);
     _CellAreaToHectares                          = repast::strToDouble(props->getProperty("CohortParameters.CellAreaToHectares"));

     //metabolim
     //Ectotherms

     _MetabolismMassExponentEcto  = repast::strToDouble(props->getProperty("CohortParameters.MetabolismMassExponentEcto"));
     _NormalizationConstantEcto   = repast::strToDouble(props->getProperty("CohortParameters.NormalizationConstantEcto"));
     _ActivationEnergyEcto        = repast::strToDouble(props->getProperty("CohortParameters.ActivationEnergyEcto")); 
     _NormalizationConstantBMR    = repast::strToDouble(props->getProperty("CohortParameters.NormalizationConstantBMR"));
     _BasalMetabolismMassExponent = repast::strToDouble(props->getProperty("CohortParameters.BasalMetabolismMassExponent"));
     _EnergyScalarEcto            = repast::strToDouble(props->getProperty("CohortParameters.EnergyScalarEcto"));
    
    //Endotherms
     _MetabolismMassExponentEndo  = repast::strToDouble(props->getProperty("CohortParameters.MetabolismMassExponentEndo"));
     _NormalizationConstantEndo   = repast::strToDouble(props->getProperty("CohortParameters.NormalizationConstantEndo"));
     _ActivationEnergyEndo        = repast::strToDouble(props->getProperty("CohortParameters.ActivationEnergyEndo")); 
     _BoltzmannConstant           = repast::strToDouble(props->getProperty("CohortParameters.BoltzmannConstant"));
     _EnergyScalarEndo            = repast::strToDouble(props->getProperty("CohortParameters.EnergyScalarEndo"));
     _TemperatureUnitsConvert     = repast::strToDouble(props->getProperty("CohortParameters.TemperatureUnitsConvert"));

    //reproduction
    
     _MassRatioThreshold                = repast::strToDouble(props->getProperty("CohortParameters.MassRatioThreshold"));
     _MassEvolutionProbabilityThreshold = repast::strToDouble(props->getProperty("CohortParameters.MassEvolutionProbabilityThreshold"));
     _MassEvolutionStandardDeviation    = repast::strToDouble(props->getProperty("CohortParameters.MassEvolutionStandardDeviation"));
     _SemelparityAdultMassAllocation    = repast::strToDouble(props->getProperty("CohortParameters.SemelparityAdultMassAllocation"));
    
    //mortality
     _MortalityRateBackground  = repast::strToDouble(props->getProperty("CohortParameters.MortalityRateBackground"));
     _MortalityRateMature      = repast::strToDouble(props->getProperty("CohortParameters.MortalityRateMature"));
     _LogisticInflectionPoint  = repast::strToDouble(props->getProperty("CohortParameters.LogisticInflectionPoint"));
     _LogisticScalingParameter = repast::strToDouble(props->getProperty("CohortParameters.LogisticScalingParameter"));
     _MaximumStarvationRate    = repast::strToDouble(props->getProperty("CohortParameters.MaximumStarvationRate"));
}
//------------------------------------------------------------------------------------------------------------
void Cohort::ResetAccounts( ) {
    // Initialize delta abundance sorted list with appropriate processes

    _Accounting["abundance"]["mortality"] = 1.0;//NB this applies because of a change from the original model - this value is now a multiplier (reduces possibility of negatives)

    // Initialize delta biomass sorted list with appropriate processes
    _Accounting["biomass"]["metabolism"] = 0.0;
    _Accounting["biomass"]["carnivory"] = 0.0;
    _Accounting["biomass"]["herbivory"] = 0.0;
    _Accounting["biomass"]["reproduction"] = 0.0;

    // Initialize delta reproductive biomass vector with appropriate processes

    _Accounting["reproductivebiomass"]["reproduction"] = 0.0;

    // Initialize organic pool delta vector with appropriate processes
    _Accounting["organicpool"]["herbivory"] = 0.0;
    _Accounting["organicpool"]["carnivory"] = 0.0;
    _Accounting["organicpool"]["mortality"] = 0.0;

    // Initialize respiratory CO2 pool delta vector with appropriate processes
    _Accounting["respiratoryCO2pool"]["metabolism"] = 0.0;
}
//used to create an initial set of cohorts at the start of a run
//------------------------------------------------------------------------------------------------------------
void Cohort::setup(unsigned functionalGroup,unsigned numCohortsThisCell,Environment* e,randomizer* r){
    ResetAccounts( );

    _FunctionalGroupIndex=functionalGroup;
    _Merged                      = false;
    _alive                       = true;

	_Heterotroph=(CohortDefinitions::Get()->Trait(functionalGroup     , "heterotroph/autotroph")  =="heterotroph");   
    _Autotroph  =!_Heterotroph;
    _Endotherm  =(CohortDefinitions::Get()->Trait(functionalGroup     , "endo/ectotherm")         =="endotherm");
    _Ectotherm  =!_Endotherm;
    _Realm      =CohortDefinitions::Get()->Trait(functionalGroup     , "realm");

    _Iteroparous=(CohortDefinitions::Get()->Trait(functionalGroup     , "reproductive strategy")  =="iteroparity");
    _Semelparous=!_Iteroparous;
    _Herbivore=(CohortDefinitions::Get()->Trait(functionalGroup       , "nutrition source")       =="herbivore");
    _Carnivore=(CohortDefinitions::Get()->Trait(functionalGroup       , "nutrition source")       =="carnivore");
    _Omnivore= (CohortDefinitions::Get()->Trait(functionalGroup       , "nutrition source")       =="omnivore");
    _IsPlanktonic= (CohortDefinitions::Get()->Trait(functionalGroup   , "mobility")               =="planktonic");
    _IsFilterFeeder=(CohortDefinitions::Get()->Trait(functionalGroup   , "diet")                   =="allspecial");
    
    _ProportionSuitableTimeActive= CohortDefinitions::Get()->Property(functionalGroup   ,"proportion suitable time active");
    
    _IsMature=false;
    _IndividualReproductivePotentialMass=0;
    
    _AssimilationEfficiency_H=CohortDefinitions::Get()->Property(functionalGroup   ,"herbivory assimilation");
    _AssimilationEfficiency_C=CohortDefinitions::Get()->Property(functionalGroup   ,"carnivory assimilation");
    _BirthTimeStep=0;
    _MaturityTimeStep=std::numeric_limits<unsigned>::max( );
    _MinimumMass=CohortDefinitions::Get()->Property(functionalGroup   ,"minimum mass");
    _MaximumMass=CohortDefinitions::Get()->Property(functionalGroup   ,"maximum mass");

    // A bit faster to use nextDouble()  rather than repast::DoubleUniformGenerator gen = repast::Random::instance()->createUniDoubleGenerator(0, 1);gen.next();
    
    _AdultMass = pow( 10, ( repast::Random::instance()->nextDouble() * ( log10( _MaximumMass ) - log10( 50 * _MinimumMass ) ) + log10( 50 * _MinimumMass ) ) );
    

    NormalGenerator NJ= repast::Random::instance()->createNormalGenerator(0.1,0.02);
    
    _LogOptimalPreyBodySizeRatio = log(std::max( 0.01, NJ.next() ));
    
    double expectedLnAdultMassRatio = 2.24 + 0.13 * log( _AdultMass );
    //in the original code the mean and sd are those of the underlying normal distribution
    //in the boost library they refer to the log distibution - see
    //https://www.boost.org/doc/libs/1_43_0/libs/math/doc/sf_and_dist/html/math_toolkit/dist/dist_ref/dists/lognormal_dist.html
    //NB indepedent calls to create an RNG generator will follow the correct random sequence (i.e. they will not just re-create the initial random value) as this appears to use a singleton (tested -OK!)
    LogNormalGenerator LNJ= repast::Random::instance()->createLogNormalGenerator(exp(expectedLnAdultMassRatio+0.5*0.5/2.), (exp(0.5*0.5)-1)*exp(2*expectedLnAdultMassRatio+0.5*0.5));

    if( _Realm=="terrestrial" ) {
          do {
            _JuvenileMass = _AdultMass  / (1.0 + LNJ.next());
          } while( _AdultMass <= _JuvenileMass || _JuvenileMass < _MinimumMass );
    } else {
          do {
            _JuvenileMass = _AdultMass  / (1.0 + 10 *LNJ.next());
          } while( _AdultMass <= _JuvenileMass || _JuvenileMass < _MinimumMass );
    }

    double NewBiomass = ( 3300. / numCohortsThisCell ) * 100 * 3000 * pow( 0.6, ( log10( _JuvenileMass ) ) ) * ( e->Area() );

    _CohortAbundance = NewBiomass / _JuvenileMass;
    _MaximumAchievedBodyMass=_JuvenileMass;
    _IndividualBodyMass=_JuvenileMass;
    _moved=false;
    _location={0,0};
}
//------------------------------------------------------------------------------------------------------------
//Required by RHPC for cross-core copy - NB "Accounts" do not need to be included as they are instantaneous within a timestep
void Cohort::PullThingsOutofPackage( const AgentPackage& package ) {

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
  
	_Heterotroph=    package._contents._Heterotroph;   
    _Autotroph  =    !_Heterotroph;
    _Endotherm  =    package._contents._Endotherm;
    _Ectotherm  =    !_Endotherm;
    _Realm      =    package._contents._Realm;

    _Iteroparous=    package._contents._Iteroparous;
    _Semelparous=    !_Iteroparous;
    _Herbivore=      package._contents._Herbivore;
    _Carnivore=      package._contents._Carnivore;
    _Omnivore=       package._contents._Omnivore;
    _IsPlanktonic=   package._contents._IsPlanktonic;
    _IsFilterFeeder= package._contents._IsFilterFeeder;
    _MinimumMass=    package._contents._MinimumMass;
    _MaximumMass=    package._contents._MaximumMass;
    
    _ProportionSuitableTimeActive= package._contents._ProportionSuitableTimeActive;
    
    _IsMature=package._contents._IsMature;
    
    _AssimilationEfficiency_H=package._contents._AssimilationEfficiency_H;
    _AssimilationEfficiency_C=package._contents._AssimilationEfficiency_C;
    
    _moved=package._contents._moved;
    _location=package._contents._location;
    _destination=package._contents._destination;
    
}
//------------------------------------------------------------------------------------------------------------
//Required by RHPC for cross-core copy
void Cohort::PushThingsIntoPackage( AgentPackage& package ) {

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
  
	package._contents._Heterotroph=     _Heterotroph;   
    package._contents._Autotroph  =    !_Heterotroph;
    package._contents._Endotherm  =     _Endotherm;
    package._contents._Ectotherm  =    !_Endotherm;
    package._contents._Realm      =     _Realm;

    package._contents._Iteroparous=     _Iteroparous;
    package._contents._Semelparous=    !_Iteroparous;
    package._contents._Herbivore=       _Herbivore;
    package._contents._Carnivore=       _Carnivore;
    package._contents._Omnivore=        _Omnivore;
    package._contents._IsPlanktonic=    _IsPlanktonic;
    package._contents._IsFilterFeeder=  _IsFilterFeeder;
    package._contents._MinimumMass=     _MinimumMass;
    package._contents._MaximumMass=     _MaximumMass;
    
    package._contents._ProportionSuitableTimeActive=  _ProportionSuitableTimeActive;
    
    package._contents._IsMature= _IsMature;
    
    package._contents._AssimilationEfficiency_H= _AssimilationEfficiency_H;
    package._contents._AssimilationEfficiency_C= _AssimilationEfficiency_C;
    
    package._contents._moved=_moved;
    package._contents._location=_location;
    package._contents._destination=_destination;
    
}
//------------------------------------------------------------------------------------------------------------
void Cohort::setupOffspring( Cohort* actingCohort, double juvenileBodyMass, double adultBodyMass, double initialBodyMass, double initialAbundance, unsigned birthTimeStep ) {
     
     _newH=NULL;
    _FunctionalGroupIndex        = actingCohort->_FunctionalGroupIndex;
    _JuvenileMass                = juvenileBodyMass;
    _AdultMass                   = adultBodyMass;
    _IndividualBodyMass          = initialBodyMass;
    _CohortAbundance             = initialAbundance;
    _CurrentTimeStep               = birthTimeStep;
    _BirthTimeStep               = birthTimeStep;
    _MaturityTimeStep            = std::numeric_limits<unsigned>::max( );
    _LogOptimalPreyBodySizeRatio = actingCohort->_LogOptimalPreyBodySizeRatio;
    _MaximumAchievedBodyMass     = juvenileBodyMass;
    _Merged                      = false;
    _alive                       = true;
    _IndividualReproductivePotentialMass = 0;
  
	_Heterotroph=    actingCohort->_Heterotroph;   
    _Autotroph  =    !_Heterotroph;
    _Endotherm  =    actingCohort->_Endotherm;
    _Ectotherm  =    !_Endotherm;
    _Realm      =    actingCohort->_Realm;

    _Iteroparous=    actingCohort->_Iteroparous;
    _Semelparous=    !_Iteroparous;
    _Herbivore=      actingCohort->_Herbivore;
    _Carnivore=      actingCohort->_Carnivore;
    _Omnivore=       actingCohort->_Omnivore;
    _IsPlanktonic=   actingCohort->_IsPlanktonic;
    _IsFilterFeeder= actingCohort->_IsFilterFeeder;
    _MinimumMass=    actingCohort->_MinimumMass;
    _MaximumMass=    actingCohort->_MaximumMass;
    
    _ProportionSuitableTimeActive= actingCohort->_ProportionSuitableTimeActive;
    
    _IsMature=false;
    
    _AssimilationEfficiency_H=CohortDefinitions::Get()->Property(_FunctionalGroupIndex   ,"herbivory assimilation");
    _AssimilationEfficiency_C=CohortDefinitions::Get()->Property(_FunctionalGroupIndex   ,"carnivory assimilation");
    
    _moved=false;
    _location=actingCohort->_location;
    _destination=_location;


}
//------------------------------------------------------------------------------------------------------------
void Cohort::step(Environment* e,vector<Cohort*>& preys,vector<Stock*>& stocks,const unsigned Timestep,MadModel* m) {
    _newH=NULL;//make sure the reproduction pointer has been zeroed out

    if (_CohortAbundance - Parameters::Get( )->GetExtinctionThreshold( ) <= 0)return;
    _CurrentTimeStep=Timestep;

    //note - passing in preys here from above ensures that new cohorts created later do not immediately get eaten.
    //since they *do* get added to the global cohort list straight away.
    //pass in of environment avoids having to find it with an expensive lookup
    ResetAccounts( );
    for (auto A: _Accounting)for (auto B: A.second)if(B.first!="mortality")assert(B.second==0);
    assignTimeActive(e);
    eat(e,preys,stocks,m);
    metabolize(e);
    reproduce(e);
    mort();
    applyEcology(e);

}
//------------------------------------------------------------------------------------------------------------
void Cohort::markForDeath(){
    if (_CohortAbundance - Parameters::Get( )->GetExtinctionThreshold( ) <= 0 || _IndividualBodyMass <= 0){ 

      //mark the cohort but don't kill it yet to avoid any problems with movement code in parallel
      _alive=false;
    }
}
//------------------------------------------------------------------------------------------------------------
void Cohort::moveIt(Environment* e,MadModel* m){

        _moved=false;

        if (!_alive)return;
        _destination=_location;
        vector<int> movement={0,0};
       // Calculate the scalar to convert from the time step units used by this implementation of dispersal to the global model time step units
        double DeltaT = Constants::cMonth;
        double CellArea      = e->Area();

        
        double dispersalSpeed=0.;
        if( _IsPlanktonic || (_Realm=="marine" &&  (_IndividualBodyMass <= Parameters::Get( )->GetPlanktonSizeThreshold( ))))  {
          NormalGenerator NJ= repast::Random::instance()->createNormalGenerator(0,1);
          // Advective dispersal
          //dispersalName = "advective";

          for( int mm = 0; mm < _AdvectionTimeStepsPerModelTimeStep; mm++ ) {
            // Get the u speed and the v speed from the cell data
            double uAdvectiveSpeed = e->uVel();
            assert( uAdvectiveSpeed > -9999 );

            double vAdvectiveSpeed = e->vVel();
            assert( vAdvectiveSpeed > -9999 );

            // Note that this formulation drops the delta t because we set the horizontal diffusivity to be at the same temporal scale as the time step
            
            // Calculate the distance travelled in this dispersal (not global) time step. both advective and diffusive speeds need to have been converted to km / advective model time step
            double uSpeed = uAdvectiveSpeed * _VelocityUnitConversion / _AdvectionTimeStepsPerModelTimeStep + NJ.next() * sqrt( ( 2.0 * _HorizontalDiffusivityKmSqPerADTimeStep ) );
            double vSpeed = vAdvectiveSpeed * _VelocityUnitConversion / _AdvectionTimeStepsPerModelTimeStep + NJ.next() * sqrt( ( 2.0 * _HorizontalDiffusivityKmSqPerADTimeStep ) );
            TryToDisperse( uSpeed,vSpeed,e,m );

          }
        }// Otherwise, if mature do responsive dispersal

        else if( _IsMature ) {

          //dispersalName = "responsive";

          dispersalSpeed=_DispersalSpeedBodyMassScalar * pow( _AdultMass, _DispersalSpeedBodyMassExponent);

          // Check for starvation-driven dispersal
          // A boolean to check whether a cohort has dispersed
          bool cohortHasDispersed = false;

          // Check for starvation driven dispersal. Note that at present we are just tracking starvation for adults
          // Assume a linear relationship between probability of dispersal and body mass loss, up to _StarvationDispersalBodyMassThreshold at which point the cohort will try to disperse every time step
          if( _IndividualBodyMass < _AdultMass ) {

             // If the body mass loss is greater than the starvation dispersal body mass threshold, then the cohort tries to disperse
             if( _IndividualBodyMass / _AdultMass < _StarvationDispersalBodyMassThreshold ) {
                // Cohort tries to disperse
                TryToDisperse( dispersalSpeed,e,m );
                // Note that regardless of whether or not it succeeds,  it is counted as having dispersed for the purposes of not then allowing it to disperse based on its density.
                cohortHasDispersed = true;
                // Otherwise, the cohort has a chance of trying to disperse proportional to its mass lass
             } else {
               // Cohort tries to disperse with a particular probability
               if( ( ( 1.0 - _IndividualBodyMass / _AdultMass ) / ( 1.0 - _StarvationDispersalBodyMassThreshold ) ) > repast::Random::instance()->nextDouble() ) {
                 TryToDisperse( dispersalSpeed,e,m );
                 cohortHasDispersed = true;
               }
             }
          }
          //not starvation - density driven?
          if( !cohortHasDispersed ) {
            // If below the density threshold
            if( ( _CohortAbundance / CellArea ) < _DensityThresholdScaling / _AdultMass ) {
                TryToDisperse( dispersalSpeed,e,m );

           }
          }
        }// If the cohort is immature, run diffusive dispersal
        else {
           dispersalSpeed=_DispersalSpeedBodyMassScalar * pow( _IndividualBodyMass, _DispersalSpeedBodyMassExponent);

                TryToDisperse( dispersalSpeed,e,m );
        }
        
      
}
//------------------------------------------------------------------------------------------------------------
void Cohort::TryToDisperse(double dispersalSpeed, Environment* e,MadModel* m){
    double randomDirection = repast::Random::instance()->nextDouble()* 2 * _Pi;

    // Calculate the u and v components given the dispersal speed
    double uSpeed = dispersalSpeed * cos( randomDirection );
    double vSpeed = dispersalSpeed * sin( randomDirection );
    TryToDisperse(uSpeed, vSpeed,e,m);
 }
  //------------------------------------------------------------------------------------------------------------
void Cohort::TryToDisperse(double uSpeed, double vSpeed,Environment* e, MadModel* m){

      vector<double> uv={0,0};//default to no dispersal
      if (m->_dispersalSelection=="probabilistic")uv=dProb(uSpeed,vSpeed,e);
      else 
      if (m->_dispersalSelection=="direct") uv=dDirect(uSpeed,vSpeed,e);
      
      double signu=uv[0],signv=uv[1]; 
      double x=_destination[0],y=_destination[1];//need to accumulate this over multiple steps or some movements might get missed/not wrap properly

     //only move if realm matches and we haven't exceeded upper and lower latitude bounds (currently no movement across the pole)
     int yw=floor(y+signv);
     if (!m->_noLongitudeWrap){
        if (x + signu < m->_minX){x = x + (m->_maxX - m->_minX);}
        if (x + signu > m->_maxX){x = x - (m->_maxX - m->_minX);}
     }
     int xw=floor(x+signu);
     if (yw >= m->_minY && yw <= m->_maxY){
         if (xw >= m->_minX && xw <= m->_maxX) {
           Environment* E=m->_Env[xw - m->_minX + (m->_maxX - m->_minX + 1)*(yw - m->_minY)];     // get environment of destination cell
           if (E->_Realm==_Realm || _Realm=="all"){// no movement if wrong realm at destination
               _destination[0]=x+signu;_destination[1]=y+signv;
               if (xw!=floor(_location[0]) || yw!=floor(_location[1]))_moved=true;//special treatment needed if we have changed cell
           }          
         }
     } 
   }
//------------------------------------------------------------------------------------------------------------
 vector<double> Cohort::dProb(double uSpeed, double vSpeed,Environment* e){
     //original madingely model dispersal when cohorts are only able to be present at fixed integer cell co-ordinates
    double latCellLength = e->Height();
    double lonCellLength = e->Width();
    double CellArea      = e->Area();

    double AreaOutsideBoth = abs( uSpeed * vSpeed );

    // Calculate the area of the grid cell that is now outside in the u direction (not including the diagonal)
    double AreaOutsideU = abs( uSpeed * latCellLength ) - AreaOutsideBoth;

    // Calculate the proportion of the grid cell that is outside in the v direction (not including the diagonal)
    double AreaOutsideV = abs( vSpeed * lonCellLength ) - AreaOutsideBoth;

    // Convert areas to a probability
    double DispersalProbability = ( AreaOutsideU + AreaOutsideV + AreaOutsideBoth ) / CellArea;
    // Check that we don't have any issues - this is problematic for small cell sizes as the "probability" can often exceed 1.
    if( DispersalProbability > 1 ){
       DispersalProbability=1.;
    }
    // Check to see in which axis the cohort disperses

   // Note that the values in the dispersal array are the proportional area moved outside the grid cell in each direction; we simply compare the random draw to this
   // to determine the direction in which the cohort moves probabilistically
   double RandomValue=repast::Random::instance()->nextDouble();
   double signu=0,signv=0;

   if( DispersalProbability >= RandomValue ) {

      signu = ( uSpeed > 0 ) - ( uSpeed < 0 );
      signv = ( vSpeed > 0 ) - ( vSpeed < 0 );

      // Longitudinally
      if( RandomValue <= AreaOutsideU / CellArea ) {
       signv = 0;
      } else {
       //Latitudinally
       if( RandomValue <= ( AreaOutsideU / CellArea + AreaOutsideV / CellArea ) ) {
           signu = 0;
       }
     }
   }
   return vector<double>({signu,signv});
 }
//------------------------------------------------------------------------------------------------------------
 vector<double> Cohort::dDirect(double uSpeed, double vSpeed,Environment* e){
    //improved dispersla where cohorts can take on fractional cell co-ordinates
    //this is *required* when cross-cell interaction becomes important
    // Calculate the fraction of the grid cell in the u direction 
    double ufrac = ( uSpeed / e->Width() );

    // Calculate the fraction of the grid cell in the v direction
    double vfrac = ( vSpeed / e->Height() );
    return vector<double>({ufrac,vfrac});
 }
//------------------------------------------------------------------------------------------------------------
void Cohort::assignTimeActive(Environment* e){

    double WarmingTolerance = 0;
    double ThermalSafetyMargin = 0;
    double Topt = 0;
    double CTmax = 0;
    double CTmin = 0;
    double AmbientTemp = 0;
    double DTR = 0;


    //Only work on heterotroph cohorts
    if( _Heterotroph ) {
        //Check if this is an endotherm or ectotherm
        if( _Endotherm ) {
            //Assumes the whole timestep is suitable for endotherms to be active - actual time active is therefore the proportion specified for this functional group.
            _ProportionTimeActive = _ProportionSuitableTimeActive;;
        } else {
            //If ectotherm then use realm specific function
            if( _Realm=="marine"){
                 double ProportionTimeSuitableMarine = 1.0;
                _ProportionTimeActive = ProportionTimeSuitableMarine * _ProportionSuitableTimeActive;
            } else {
               AmbientTemp = e->Temperature();
               DTR         = e->DiurnalTemperatureRange();

               //Calculate the Warming tolerance and thermal safety margin given standard deviation of monthly temperature
               WarmingTolerance =    _TerrestrialWarmingToleranceSlope * e->SDTemperature() + _TerrestrialWarmingToleranceIntercept;
               ThermalSafetyMargin = _TerrestrialTSMSlope *              e->SDTemperature() + _TerrestrialTSMIntercept;

               Topt  = ThermalSafetyMargin + e->AnnualTemperature();
               CTmax = WarmingTolerance    + e->AnnualTemperature();

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
               double POver = ( ( _Pi / 2.0 ) - asin( temp ) ) / _Pi;

               if( CTmin - DTmax > 0.0 ) {
                  temp = 1.0;
               } else if( CTmin - DTmin < 0.0 ) {
                  temp = -1.0;
               } else {
                  temp = 2 * ( CTmin - AmbientTemp ) / DTR;
               }
               //Proportion of time for which ambient temperatures are below the critical lower temperature

               double PBelow = 1 - ( ( _Pi / 2.0 ) - asin( temp ) ) / _Pi;

               double  ProportionDaySuitable = 1 - ( POver + PBelow );
               double ProportionTimeSuitableTerrestrial= ProportionDaySuitable;
               _ProportionTimeActive = ProportionTimeSuitableTerrestrial * _ProportionSuitableTimeActive;

            }
        }
    }
}
//------------------------------------------------------------------------------------------------------------
double Cohort::distance(MadAgent* a1, MadAgent* a2,MadModel* m){
    //simply - the number of cell widths, Manhattan style
    double wrappedDistX=abs(a1->_location[0]-a2->_location[0]);
    if (!m->_noLongitudeWrap)wrappedDistX=min(wrappedDistX,(m->_maxX - m->_minX+1)-wrappedDistX);
    return max(wrappedDistX,abs( a1->_location[1] - a2->_location[1]));
}
//------------------------------------------------------------------------------------------------------------
void Cohort::eat(Environment* e,vector<Cohort*>& preys,vector<Stock*>& stocks,MadModel* m){
    
    double  PotentialBiomassEaten=0.,HandlingTime=0.,HandlingTimeScaled=0.,BiomassesEaten=0.,IndividualHerbivoryRate=0.;
    double edibleFraction=0,AttackRateExponent=0, HandlingTimeExponent=0,HandlingTimeScalar=0;
    double _DaysInATimeStep=Constants::cDay;
    double PredatorAbundanceMultipliedByTimeEating=0;
    double CellAreaHectares=e->Area() * _CellAreaToHectares;

    if (_Herbivore || _Omnivore){

       //acumulate handling time - these values will be for herbies - note the issue if these are stored in the class with omnis
       if (_Realm=="marine"){
       edibleFraction       = _edibleFractionMarine;
       AttackRateExponent   = _AttackRateExponentMarine;
       HandlingTimeExponent = _HandlingTimeExponentMarine;
       HandlingTimeScalar   = _HandlingTimeScalarMarine;
       }
       if(_Realm=="terrestrial"){

       edibleFraction       = _edibleFractionTerrestrial;
       AttackRateExponent   = _AttackRateExponentTerrestrial;
       HandlingTimeExponent = _HandlingTimeExponentTerrestrial;
       HandlingTimeScalar   = _HandlingTimeScalarTerrestrial;
       }
       IndividualHerbivoryRate               = _HerbivoryRateConstant    * pow( _IndividualBodyMass, ( _HerbivoryRateMassExponent ) );
       for (auto& stock: stocks){
        if (distance(this,stock,m)<1){
         if( stock->_TotalBiomass > 0.0 ) {
           PotentialBiomassEaten  = IndividualHerbivoryRate   * pow( (stock->_TotalBiomass*edibleFraction) / CellAreaHectares, AttackRateExponent ); // could store this to avoid second calc.?
           HandlingTimeScaled     = HandlingTimeScalar        * pow( ( _ReferenceMass / _IndividualBodyMass ), HandlingTimeExponent );
           HandlingTime          += PotentialBiomassEaten* HandlingTimeScaled; 
         }
        }
       }

    }
    std::map<int, vector<double> > BinnedPreyDensities;
    if (_Carnivore || _Omnivore){
         if( _CohortAbundance > 0 ) {

        // Pre-calculate individual values for this predator to speed things up
        double _MaximumSearchRate = _SearchRateConstant * _IndividualBodyMass; // Ref: Predator-only component of eq. 35 from text S1.?? should be pow(mass, massRateExponent)
        double _DaysEating = _DaysInATimeStep * _ProportionTimeActive;

        double LogOptimalPreySizeRatio = _LogOptimalPreyBodySizeRatio;

        // If a filter feeder, then optimal body size is a value not a ratio: convert it to a ratio to ensure that all calculations work correctly
        if( _IsFilterFeeder ) {
            // Optimal body size is actually a value, not a ratio, so convert it to a ratio based on the present body size
            LogOptimalPreySizeRatio = std::log( std::exp( _LogOptimalPreyBodySizeRatio ) / _IndividualBodyMass ); //LogOptimalPreyBodySizeRatio-log(actingCohort.IndividualBodyMass);
        }
        // Calculate the reference mass scaling ratio
        HandlingTimeScaled = _HandlingTimeScalar_C * std::pow( _ReferenceMass / _IndividualBodyMass, _HandlingTimeExponent_C );
        PredatorAbundanceMultipliedByTimeEating = _CohortAbundance * _DaysEating;

        // Calculate the abundance of prey in each of the prey mass bins

        // Loop through prey functional groups


        for (auto& prey:preys){
            if (distance(this,prey,m)<1){
             // Calculate the difference between the actual body size ratio and the optimal ratio, 
             // and then divide by the standard deviation in log ratio space to determine in 
             // which bin to assign the prey item.
             if (BinnedPreyDensities.count(prey->_FunctionalGroupIndex)==0){BinnedPreyDensities[prey->_FunctionalGroupIndex].resize(_NumberOfBins);
                  for( unsigned binIndex = 0; binIndex < _NumberOfBins; binIndex++ ) BinnedPreyDensities[ prey->_FunctionalGroupIndex][ binIndex ] = 0;}

             if( prey->_IndividualBodyMass > 0 ) {

                int binIndex = ( int )( ( ( std::log( prey->_IndividualBodyMass / _IndividualBodyMass ) - LogOptimalPreySizeRatio ) / ( 0.5 * _FeedingPreferenceStandardDeviation ) ) + ( _NumberOfBins / 2 ) );


                if( ( 0 < binIndex ) && ( binIndex < _NumberOfBins ) ) {
                    BinnedPreyDensities[ prey->_FunctionalGroupIndex ][ binIndex ] += prey->_CohortAbundance / CellAreaHectares;

                }
             }
            }
        }

        // Loop over potential prey functional groups
        for (auto& prey: preys){
            if (distance(this,prey,m)<1){
                prey->_PotentialAbundanceEaten = 0;
                //No Cannibalism
                if( (prey->getId().id() != getId().id()) &&  prey->_IndividualBodyMass > 0) {
                   unsigned binIndex = ( unsigned )( ( ( std::log( prey->_IndividualBodyMass / _IndividualBodyMass ) - LogOptimalPreySizeRatio ) / ( 0.5 * _FeedingPreferenceStandardDeviation ) ) + ( _NumberOfBins / 2 ) );
                    if( (_IsFilterFeeder &&  prey->_IsPlanktonic ) || !_IsFilterFeeder) {

                        if(   ( 0 < binIndex ) && ( binIndex < _NumberOfBins ) ) {
                            // Calculate the potential abundance from this cohort eaten by the acting cohort
                            
                                // Calculate the relative feeding preference from a log-normal distribution with mean equal to the optimal prey to predator ratio and standard deviation as specified
   
                            double preferenceForPrey = std::exp( -( std::pow( ( ( std::log( prey->_IndividualBodyMass / _IndividualBodyMass ) - LogOptimalPreySizeRatio ) / _FeedingPreferenceStandardDeviation ), 2 ) ) ); // Ref: Text S1, eq 36.
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
    }
    if (_Herbivore || _Omnivore){
       //accumulate eaten
       for (auto& stock: stocks){
           if (distance(this,stock,m)<1){
            double InstantFractionEaten = 0.0;
            double EdibleMass=0;
            if( stock->_TotalBiomass > 0.0 ) {
              PotentialBiomassEaten               = IndividualHerbivoryRate   * pow( (stock->_TotalBiomass*edibleFraction) / CellAreaHectares, AttackRateExponent );
              EdibleMass                          = stock->_TotalBiomass*edibleFraction;
              InstantFractionEaten                = _CohortAbundance * ( ( PotentialBiomassEaten / ( 1 + HandlingTime ) ) / EdibleMass );
              BiomassesEaten                      = EdibleMass * ( 1 - exp( -InstantFractionEaten * Constants::cDay*_ProportionTimeActive ) ); // should be min(stock._TotalBiomass,...)
              stock->_TotalBiomass               = stock->_TotalBiomass - BiomassesEaten;
            } 

            // Return the total  biomass of the autotroph stock eaten
       
           if( _CohortAbundance > 0 )_Accounting["biomass"]["herbivory"] += BiomassesEaten * _AssimilationEfficiency_H / _CohortAbundance;

   
           _Accounting["organicpool"]["herbivory"] += BiomassesEaten * ( 1 - _AssimilationEfficiency_H );
          }
       }
    }
    if (_Carnivore || _Omnivore){
    if( _CohortAbundance > 0 ) {
        double biomassConsumed = 0;
        // Loop over potential prey functional groups
        int counter=0;
        for (auto& prey: preys){
           if (distance(this,prey,m)<1){
              counter++;
                // Calculate the actual abundance of prey eaten from this cohort
                double _AbundancesEaten = 0;
                if( prey->_CohortAbundance > 0 && prey->_PotentialAbundanceEaten>0) {
                    // Calculate the actual abundance of prey eaten from this cohort
                    _AbundancesEaten = prey->_CohortAbundance * ( 1 - std::exp( -( PredatorAbundanceMultipliedByTimeEating * ( ( prey->_PotentialAbundanceEaten / ( HandlingTime + 1 ) ) / prey->_CohortAbundance ) ) ) );
                }
                // Remove number of prey eaten from the prey cohort

                prey->_CohortAbundance -= _AbundancesEaten;

                biomassConsumed += ( prey->_IndividualBodyMass + prey->_IndividualReproductivePotentialMass ) * _AbundancesEaten / _CohortAbundance; //per capita

         }
        }

        // Add the biomass eaten and assimilated by an individual to the delta biomass for the acting (predator) cohort
        _Accounting[ "biomass" ][ "carnivory" ] = biomassConsumed * _AssimilationEfficiency_C;
        // Move the biomass eaten but not assimilated by an individual into the organic matter pool
        _Accounting[ "organicpool" ][ "carnivory" ] = biomassConsumed * ( 1 - _AssimilationEfficiency_C ) * _CohortAbundance;
    }
    }

}
//------------------------------------------------------------------------------------------------------------
void Cohort::metabolize(Environment* e){
    if (_Ectotherm){
    

    // Calculate the scalar to convert from the time step units used by this implementation of metabolism to the global  model time step units
    double DeltaT = Constants::cDay;
    double temperature=e->Temperature() + _TemperatureUnitsConvert;
    bool ProportionTimeActiveCalculatedThisTimestep = false;
    double FieldMetabolicLosskJ = _NormalizationConstantEcto    * pow( _IndividualBodyMass, _MetabolismMassExponentEcto )      * exp( -( _ActivationEnergyEcto / ( _BoltzmannConstant * temperature ) ) );
    double BasalMetabolicLosskJ = _NormalizationConstantBMR * pow( _IndividualBodyMass, _BasalMetabolismMassExponent ) * exp( -( _ActivationEnergyEcto / ( _BoltzmannConstant * temperature ) ) );

    // metabolic loss in grams
    double IndividualMetabolicRate = (( _ProportionTimeActive * FieldMetabolicLosskJ ) + ( ( 1 - _ProportionTimeActive ) * ( BasalMetabolicLosskJ ) ) ) * _EnergyScalarEcto;
    _Accounting["biomass"]["metabolism"] = -IndividualMetabolicRate * DeltaT;

    // If metabolic loss is greater than individual body mass after herbivory and predation, then set equal to individual body mass
    _Accounting["biomass"]["metabolism"] = std::max( _Accounting["biomass"]["metabolism"], -( _IndividualBodyMass + _Accounting["biomass"]["carnivory"] + _Accounting["biomass"]["herbivory"] ) );

    // Add total metabolic loss for all individuals in the cohort to delta biomass for metabolism in the respiratory CO2 pool
    _Accounting["respiratoryCO2pool"]["metabolism"] = -_Accounting["biomass"]["metabolism"] * _CohortAbundance;

    }
    if (_Endotherm){


    // Assume all endotherms have a constant body temperature of 37degC
    double EndothermBodyTemperature = 37.0 + _TemperatureUnitsConvert;
    // Calculate the scalar to convert from the time step units used by this implementation of metabolism to the global  model time step units
    double DeltaT = Constants::cDay;
    // Calculate metabolic loss in kJ
    double metabolicLosskJ = _NormalizationConstantEndo * pow( _IndividualBodyMass, _MetabolismMassExponentEndo ) * exp( -( _ActivationEnergyEndo / ( _BoltzmannConstant * EndothermBodyTemperature ) ) );
    // metabolic loss in grams
    double IndividualMetabolicRate= metabolicLosskJ * _EnergyScalarEndo;
    // Calculate metabolic loss for an individual and add the value to the delta biomass for metabolism
    _Accounting[ "biomass" ][ "metabolism" ] = -IndividualMetabolicRate * DeltaT;

    // If metabolic loss is greater than individual body mass after herbivory and predation, then set equal to individual body mass
    _Accounting[ "biomass" ][ "metabolism" ] = std::max( _Accounting[ "biomass" ][ "metabolism" ], -( _IndividualBodyMass + _Accounting[ "biomass" ][ "carnivory" ] + _Accounting[ "biomass" ][ "herbivory" ] ) );

    // Add total metabolic loss for all individuals in the cohort to delta biomass for metabolism in the respiratory CO2 pool
    _Accounting[ "respiratoryCO2pool" ][ "metabolism" ] = -_Accounting[ "biomass" ][ "metabolism" ] * _CohortAbundance;
    }
    //heterotroph - in the original model but not used */
}
//------------------------------------------------------------------------------------------------------------
void Cohort::reproduce(Environment* e){
    _newH=NULL;//memory allocation for new cohorts is dealt with via the model list of agents 
    // Biomass per individual in each cohort to be assigned to reproductive potential
    double BiomassToAssignToReproductivePotential;

    // Net biomass change from other ecological functions this time step
    double NetBiomassFromOtherEcologicalFunctionsThisTimeStep;

    // Reset variable holding net biomass change of individuals in this cohort as a result of other ecological processes
    NetBiomassFromOtherEcologicalFunctionsThisTimeStep = 0.0;

    // Loop over all items in the biomass deltas
    for( auto& Biomass: _Accounting[ "biomass" ] ) {
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
        _Accounting[ "reproductivebiomass" ][ "reproduction" ] += BiomassToAssignToReproductivePotential;
        _Accounting[ "biomass" ][ "reproduction" ] -= BiomassToAssignToReproductivePotential;

    } else {
        // Cohort has not gained sufficient biomass to assign any to reproductive potential, so take no action
    }

    // Run reproductive events. Note that we can't skip juveniles here as they could conceivably grow to adulthood and get enough biomass to reproduce in a single time step
    // due to other ecological processes

    const std:: string _TimeUnitImplementation = "month";
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

    for( auto& Biomass: _Accounting[ "biomass" ] ) {
        // Add the delta biomass to net biomass
        bodyMassIncludingChangeThisTimeStep += Biomass.second;

    }
    bodyMassIncludingChangeThisTimeStep += _IndividualBodyMass;

    // Calculate the reproductive biomass of an individual in this cohort including changes this time step from other ecological processes  
    reproductiveMassIncludingChangeThisTimeStep = 0.0;

    for( auto& ReproBiomass: _Accounting[ "reproductivebiomass" ] ) {
        // Add the delta reproductive biomass to net biomass
        reproductiveMassIncludingChangeThisTimeStep += ReproBiomass.second;
    }
    reproductiveMassIncludingChangeThisTimeStep += _IndividualReproductivePotentialMass;

    if ( _IndividualBodyMass > 1.e-200 ) {
        // Get the current ratio of total individual mass (including reproductive potential) to adult body mass
        currentMassRatio = ( bodyMassIncludingChangeThisTimeStep + reproductiveMassIncludingChangeThisTimeStep ) / _AdultMass;

        // Must have enough mass to hit reproduction threshold criterion, and either (1) be in breeding season, or (2) be a marine cell (no breeding season in marine cells)
        if ( ( currentMassRatio > _MassRatioThreshold ) && (e->Breeding_Season() == 1.0  ||  _Realm=="marine"  ) ) {
            // Iteroparous and semelparous organisms have different strategies

            if ( _Iteroparous ) {
                // Iteroparous organisms do not allocate any of their current non-reproductive biomass to reproduction
                adultMassLost = 0.0;

                // Calculate the number of offspring that could be produced given the reproductive potential mass of individuals
                offspringCohortAbundance = _CohortAbundance * reproductiveMassIncludingChangeThisTimeStep / _JuvenileMass;


            } else {
                // Semelparous organisms allocate a proportion of their current non-reproductive biomass (including the effects of other ecological processes) to reproduction
                adultMassLost = _SemelparityAdultMassAllocation * bodyMassIncludingChangeThisTimeStep;

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
                         
            if( RandomValue > _MassEvolutionProbabilityThreshold ) {
             // Determine the new juvenile body mass
             NormalGenerator NJ= repast::Random::instance()->createNormalGenerator(_JuvenileMass,_MassEvolutionStandardDeviation * _JuvenileMass);
             double RandomValueJ = NJ.next();
             newJuvenileMass = std::max( RandomValueJ, _MinimumMass );

             // Determine the new adult body mass
             NormalGenerator NA= repast::Random::instance()->createNormalGenerator(_AdultMass,_MassEvolutionStandardDeviation * _AdultMass);

             double RandomValueA = NA.next();
             newAdultMass = std::min( RandomValueA, _MaximumMass );

             } 

             // Update cohort abundance in case juvenile mass has been altered through 'evolution'
             offspringCohortAbundance = offspringCohortAbundance * ( _JuvenileMass / newJuvenileMass );
             //take care to get this right using _id - otherwise rank or agent type might not be correct.
             repast::AgentId id(Cohort::_NextID, _id.currentRank(), _id.agentType());

             _newH = new Cohort(id);
             _newH->setupOffspring( this, newJuvenileMass, newAdultMass, newJuvenileMass, offspringCohortAbundance, _CurrentTimeStep);

            // Subtract all of the reproductive potential mass of the parent cohort, which has been used to generate the new
            // cohort, from the delta reproductive potential mass and delta adult body mass
            _Accounting["reproductivebiomass"]["reproduction"] -= reproductiveMassIncludingChangeThisTimeStep;
            _Accounting["biomass"]["reproduction"] -= adultMassLost;
        } else {

            // Organism is not large enough, or it is not the breeding season, so take no action
        }
    }
}

//------------------------------------------------------------------------------------------------------------
void Cohort::mort(){

    // Variables to hold the mortality rates
    double MortalityBackground;
    double MortalitySenescence;
    double MortalityStarvation;

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
    for( auto Biomass: _Accounting[ "biomass" ] ) {
        // Add the delta biomass to net biomass
        BodyMassIncludingChangeThisTimeStep += Biomass.second;
    }

    BodyMassIncludingChangeThisTimeStep = std::min( _AdultMass, BodyMassIncludingChangeThisTimeStep + _IndividualBodyMass );

    // Temporary variable to hold net reproductive biomass change of individuals in this cohort as a result of other ecological processes
    ReproductiveMassIncludingChangeThisTimeStep = 0.0;

    // Loop over all items in the biomass Cohort::Deltas
    for( auto Biomass: _Accounting[ "reproductivebiomass" ] ) {
        // Add the delta biomass to net biomass
        ReproductiveMassIncludingChangeThisTimeStep += Biomass.second;
    }
    std::string TimeUnitImplementation = "Day";
    // Calculate the scalar to convert from the time step units used by this implementation of mortality to the global model time step units
    double DeltaT = Constants::cDay;
    
    
    ReproductiveMassIncludingChangeThisTimeStep += _IndividualReproductivePotentialMass;
    // Check to see if the cohort has already been killed by predation etc
    if( BodyMassIncludingChangeThisTimeStep < 1.e-15 ) {
        //MB a small number ! maybe should be larger? (e.g. min cohort body mass))
        //This causes a difference between C# and c++ versions as there is a rounding error issue - changed in C# code to match
        // If individual body mass is not greater than zero, then all individuals become extinct
        MortalityTotal = 0.0;
        // used to be actingCohort->mCohortAbundance, but this leads to a large cancellation in applyEcology
        //so mortality total is changed here  - now we only multiply by mortality total - so values can never be negative
        //BodyMassIncludingChangeThisTimeStep = 0; //MB would be a kludge to exclude negative values below - need mass checking throughout the code
    } else {
        // Calculate background mortality rate

        MortalityBackground = _MortalityRateBackground * DeltaT;

        // If the cohort has matured, then calculate senescence mortality rate, otherwise set rate to zero
        if( _MaturityTimeStep < std::numeric_limits< unsigned >::max( ) ) {


          double TimeToMaturity = _MaturityTimeStep - _BirthTimeStep;

          // Calculate how many model time steps since the cohort reached maturity
          double AgePostMaturity = _CurrentTimeStep - _MaturityTimeStep;

          // Calculate the time since maturity as a fraction of the time that it took the cohort to reach maturity
          double FractionalAgePostMaturity = AgePostMaturity / ( TimeToMaturity + 1 );

          // Calculate the mortality rate per mortality formulation time step as a function of the exponential of the previous fraction
          double AgeRelatedMortalityRate = _MortalityRateMature * exp( FractionalAgePostMaturity );

          // Convert the mortality rate from formulation time step units to model time step units
          MortalitySenescence= AgeRelatedMortalityRate * DeltaT;

        } else {
            MortalitySenescence = 0.0;
        }
        // Calculate the starvation mortality rate based on individual body mass and maximum body mass ever
        // achieved by this cohort
        double MortalityRateStarvation = 0.;

        if( BodyMassIncludingChangeThisTimeStep - _MaximumAchievedBodyMass < 1.e-15 ) {

          // Calculate the first part of the relationship between body mass and mortality rate
          double k = -( BodyMassIncludingChangeThisTimeStep - _LogisticInflectionPoint * _MaximumAchievedBodyMass ) / ( _LogisticScalingParameter * _MaximumAchievedBodyMass );
          // Calculate mortality rate
          MortalityRateStarvation = _MaximumStarvationRate / ( 1 + exp( -k ) );
        }
        MortalityStarvation = MortalityRateStarvation * DeltaT;

        // Calculate the number of individuals that suffer mortality this time step from all sources of mortality
        MortalityTotal = exp( -MortalityBackground - MortalitySenescence - MortalityStarvation );
    }

    // Remove individuals that have died from the delta abundance for this cohort in apply ecology (by multiplication)
    _Accounting["abundance" ]["mortality"] = MortalityTotal;


    // Add the biomass of individuals that have died to the delta biomass in the organic pool (including reproductive 
    // potential mass, and mass gained through eating, and excluding mass lost through metabolism)
    _Accounting["organicpool"]["mortality"] = ( 1 - MortalityTotal ) * _CohortAbundance * ( BodyMassIncludingChangeThisTimeStep + ReproductiveMassIncludingChangeThisTimeStep );
}
//------------------------------------------------------------------------------------------------------------
void Cohort::applyEcology(Environment* e){
    // Variable to calculate net abundance change to check that cohort abundance will not become negative
    double NetAbundanceChange = 0.0;
    // Loop over all abundance deltas
    for( auto& d: _Accounting["abundance"] ) {
        // Update net abundance change
        NetAbundanceChange += d.second;
    }
    // Check that cohort abundance will not become negative
    assert( NetAbundanceChange >= 0 && "Cohort abundance < 0" );
    //this is a variation from the original below - 
    //the abundance here only changes in mort() where the change was mortalityTotal=(1-exp(-sum-of-Mortalities))*abundance
    //this is now changed to mortalityTotal=exp(-sum-of-mortalities) - avoiding possibilities of negatives arising.
    _CohortAbundance *= NetAbundanceChange;

    // Variable to calculate net biomass change to check that cohort individual body mass will not become negative
    double NetBiomass = 0.0;

    // Loop over all biomass deltas
    for( auto& d: _Accounting["biomass"] ) {
        // Update net biomass change
        NetBiomass += d.second;

    }
    double BiomassCheck = 0.0;
    bool NetToBeApplied = true;
    // If cohort abundance is greater than zero, then check that the calculated net biomass will not make individual body mass become negative
    //if( _CohortAbundance > 0 ) {

        BiomassCheck = _IndividualBodyMass + NetBiomass;

        if( BiomassCheck < 0 ) {
            //std::cout << "Biomass going negative, acting cohort: " << _FunctionalGroupIndex << ", " << getId() << " "<< NetBiomass<< " "<< _IndividualBodyMass<<std::endl;
            //exit( 1 );
            NetBiomass=-_IndividualBodyMass;
        }
    //}

    //Loop over all keys in the deltas sorted list
    for( auto& d: _Accounting["biomass"] ) {
        // If cohort abundance is zero, then set cohort individual body mass to zero and reset the biomass delta to zero, 
        // otherwise update cohort individual body mass and reset the biomass delta to zero
        if( _CohortAbundance == 0 ) {
            _IndividualBodyMass = 0.0;
        } else {
            if( NetToBeApplied ) {
                //careful! NetBiomass only gets applied once, but only because the NetToBeApplied flag is reset to false!
                //needs re-factoring
                _IndividualBodyMass = _IndividualBodyMass + NetBiomass;
                NetToBeApplied = false;
            }
        }
    }
    // Check that individual body mass is still greater than zero
    if (!(_IndividualBodyMass >=0)) cout<<_IndividualBodyMass<<endl;
    //assert( _IndividualBodyMass >= 0 && "biomass < 0" );
    // If the current individual body mass is the largest that has been achieved by this cohort, then update the maximum achieved
    // body mass tracking variable for the cohort
    if( _IndividualBodyMass > _MaximumAchievedBodyMass )
        _MaximumAchievedBodyMass = _IndividualBodyMass;

    // Variable to calculate net reproductive biomass change to check that cohort individual body mass will not become negative
    double NetReproductiveBiomass = 0.0;

    // Loop over all reproductive biomass deltas
    for( auto& d: _Accounting["reproductivebiomass"] ) {
        // Update net reproductive biomass change
        NetReproductiveBiomass += d.second;
    }

    //Loop over all keys in the abundance deltas sorted list
    for( auto& d: _Accounting["reproductivebiomass"] ) {
        // If cohort abundance is zero, then set cohort reproductive body mass to zero and reset the biomass delta to zero, 
        // otherwise update cohort reproductive body mass and reset the biomass delta to zero
        if( _CohortAbundance == 0 ) {
            _IndividualReproductivePotentialMass = 0.0;
        } else {
            _IndividualReproductivePotentialMass += d.second;
        }
    }

    updatePools(e);
    //Note that maturity time step is set in TReproductionBasic

}
//------------------------------------------------------------------------------------------------------------

void Cohort::updatePools( Environment* e) {
    // Loop over all keys in the organic pool deltas sorted list
    for( auto &D: _Accounting["organicpool"] ) {
        // Check that the delta value is not negative
        //if( D.second < 0 ) std::cout << "organic pool " << D.first << " " << D.second << std::endl;

        //assert(D.second >= 0.0 && "A delta value for the organic pool is negative " );
        // Update the organic pool biomass
        e->addToOrganicPool(D.second);
    }
    // Loop over all keys in the respiratory pool deltas sorted list
    for( auto &D: _Accounting["respiratoryCO2pool"] ) {
        // Check that the delta value is not negative
        if( D.second < 0 ) std::cout << "respiratoryCO2pool " << D.first << " " << D.second << std::endl;
        //assert( D.second >= 0.0 && "A delta value for the respiratory CO2 pool is negative" );
        // Update the respiratory CO2 pool
        e->addToRespiratoryCO2Pool(D.second);
    }
}





