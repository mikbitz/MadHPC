#include "CohortMerger.h"
#include "Parameters.h"

double CohortMerger::CalculateDistance( Cohort* cohortA, Cohort* cohortB ) {
    double AdultMassDistance = ( cohortA->_AdultMass - cohortB->_AdultMass ) / cohortA->_AdultMass;
    double JuvenileMassDistance = ( cohortA->_JuvenileMass - cohortB->_JuvenileMass ) / cohortA->_JuvenileMass;
    double CurrentMassDistance = ( cohortA->_IndividualBodyMass - cohortB->_IndividualBodyMass ) / cohortA->_IndividualBodyMass;

    return ( ( AdultMassDistance * AdultMassDistance ) + ( JuvenileMassDistance * JuvenileMassDistance ) + ( CurrentMassDistance * CurrentMassDistance ) );
}

int CohortMerger::MergeToReachThresholdFast( Environment* e) {
     unsigned maxCohorts=Parameters::Get()->GetMaximumNumberOfCohorts();
    // Set of lists of shortest distances in each functional group
    // set is automatically sorted - multiset allows for elements with the same distance
    std::multiset< CohortPair, CohortPair::Comparator > SortedDistances;
    // How many cohorts to remove to hit the threshold
    vector<Cohort*> cohorts;
  
   
    //break down cell contents into functional groups
    std::map <int, std::vector <Cohort*> > CohortsByFunctionalGroup;
    unsigned count=0;
    for (auto c: cohorts) {if (c->_alive) {CohortsByFunctionalGroup[c->_FunctionalGroupIndex].push_back(c);count++;}   }

    //this assumes we are just looking at a total over all cohorts, rather than a limit within each functional group

    int NumberToRemove = count - maxCohorts;

    int MergeCounter = 0;

    if( NumberToRemove > 0 ) {

        //Loop through functional groups
        for( auto& FunctionalGroup: CohortsByFunctionalGroup ) {

                // Loop through cohorts within functional groups

                for( int cc = 0; cc< FunctionalGroup.second.size( ) - 1; cc++ ) {
                    // Loop through comparison cohorts
                    //FunctionalGroup.second[cc]->_alive=false;
                    for( int dd = cc + 1; dd < FunctionalGroup.second.size( ); dd++ ) {
                        //unctionalGroup.second[dd]->_alive=false;
                        CohortPair PairwiseDistance( FunctionalGroup.second[cc], FunctionalGroup.second[dd], repast::Random::instance()->nextDouble() );
                        SortedDistances.insert( PairwiseDistance );
                    }
                }
            
        }

        auto I = SortedDistances.begin( );

        while( MergeCounter < NumberToRemove && I != SortedDistances.end( ) ) {
            Cohort* CohortToMergeFrom = ( I->mCohortA );
            Cohort* CohortToMergeTo =   ( I->mCohortB );
            //be careful not to merge cohorts that may already have been (randomly) merged away
            if( CohortToMergeTo->_alive && CohortToMergeFrom->_alive){
                CohortToMergeFrom->_alive=false;

                // Add the abundance of the second cohort to that of the first
                CohortToMergeTo->_CohortAbundance += CohortToMergeFrom->_CohortAbundance * CohortToMergeFrom->_IndividualBodyMass / CohortToMergeTo->_IndividualBodyMass;
                // Add the reproductive potential mass of the second cohort to that of the first
                CohortToMergeTo->_IndividualReproductivePotentialMass += CohortToMergeFrom->_IndividualReproductivePotentialMass * CohortToMergeFrom->_CohortAbundance / CohortToMergeTo->_CohortAbundance;
                // Set the abundance of the second cohort to zero
                CohortToMergeFrom->_CohortAbundance = 0.0;

                // Designate both cohorts as having merged
                I->mCohortA->_Merged = false;
                I->mCohortB->_Merged = true;
                MergeCounter++;
            }
            ++I;
        }
    }

    return MergeCounter;
}
