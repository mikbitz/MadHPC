#include "CohortMerger.h"
#include "Parameters.h"

double CohortMerger::CalculateDistance( Cohort* cohortA, Cohort* cohortB ) {
    double AdultMassDistance = ( cohortA->_AdultMass - cohortB->_AdultMass ) / cohortA->_AdultMass;
    double JuvenileMassDistance = ( cohortA->_JuvenileMass - cohortB->_JuvenileMass ) / cohortA->_JuvenileMass;
    double CurrentMassDistance = ( cohortA->_IndividualBodyMass - cohortB->_IndividualBodyMass ) / cohortA->_IndividualBodyMass;

    return ( ( AdultMassDistance * AdultMassDistance ) + ( JuvenileMassDistance * JuvenileMassDistance ) + ( CurrentMassDistance * CurrentMassDistance ) );
}

int CohortMerger::MergeToReachThresholdFast( Environment* e ) {
     unsigned maxCohorts=Parameters::Get()->GetMaximumNumberOfCohorts();
    // Set of lists of shortest distances in each functional group
    // set is automatically sorted - multiset allows for elements with the same distance
    std::multiset< CohortPair, CohortPair::Comparator > SortedDistances;
    // How many cohorts to remove to hit the threshold
    unsigned MergeCounter = 0;
    repast::relogo::AgentSet<Cohort> cohorts= e->turtlesHere<Cohort>();
    //std:: vector<Cohort*> cohorts;
    

    
    //break down cell contents into functional groups
    std::map <int, std::vector <Cohort*> > CohortsByFunctionalGroup;
    unsigned count=0;
    for (auto c: cohorts) if (c->_alive){CohortsByFunctionalGroup[c->_FunctionalGroupIndex].push_back(c);count++;}
    //this assumes we are just looking at a total over all cohorts, rather than a limit within each functional group
    unsigned NumberToRemove = count - maxCohorts;
    
    if( NumberToRemove > 0 ) {
        //Loop through functional groups
        for( auto FunctionalGroup: CohortsByFunctionalGroup ) {

                // Loop through cohorts within functional groups
                for( int cc = 0; cc< FunctionalGroup.second.size( ) - 1; cc++ ) {
                    // Loop through comparison cohorts
                    for( int dd = cc + 1; dd < FunctionalGroup.second.size( ); dd++ ) {
                        CohortPair PairwiseDistance( FunctionalGroup.second[cc], FunctionalGroup.second[dd], repast::Random::instance()->nextDouble() );
                        SortedDistances.insert( PairwiseDistance );
                    }
                }
            
        }

        auto I = SortedDistances.begin( );
        while( MergeCounter < NumberToRemove && I != SortedDistances.end( ) ) {
            Cohort* CohortToMergeFrom = ( I->mCohortA );
            Cohort* CohortToMergeTo =   ( I->mCohortB );

            if( CohortToMergeFrom->_CohortAbundance > 0 && CohortToMergeTo->_CohortAbundance > 0 ) {
                // Add the abundance of the second cohort to that of the first

                CohortToMergeTo->_CohortAbundance += CohortToMergeFrom->_CohortAbundance * CohortToMergeFrom->_IndividualBodyMass / CohortToMergeTo->_IndividualBodyMass;
                // Add the reproductive potential mass of the second cohort to that of the first
                CohortToMergeTo->_IndividualReproductivePotentialMass += CohortToMergeFrom->_IndividualReproductivePotentialMass * CohortToMergeFrom->_CohortAbundance / CohortToMergeTo->_CohortAbundance;
                // Set the abundance of the second cohort to zero
                CohortToMergeFrom->_CohortAbundance = 0.0;
                // Designate both cohorts as having merged
                CohortToMergeTo->_Merged = true;
                CohortToMergeFrom->_Merged = true;
                MergeCounter++;
            }
            ++I;
        }
    }
    
    return MergeCounter;
}
