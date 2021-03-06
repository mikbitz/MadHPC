/*
 *
 * CohortMerger.h
 *
 *  Created on: May 7, 2018
 *      Author: Mike Bithell
 * 
 *  * Emergent Global Patterns of Ecosystem Structure and Function from a Mechanistic General Ecosystem Model , Michael B. J. Harfoot , Tim Newbold, Derek P. Tittensor,  Stephen Emmott, Jon Hutton, Vassily Lyutsarev, Matthew J. Smith, Jörn P. W. Scharlemann, Drew W. Purves PLOS Biology 12(4) 2014 e1001841 https://doi.org/10.1371/journal.pbio.1001841
 */
#ifndef COHORTMERGE_H
#define COHORTMERGE_H

#include "Cohort.h"
#include "CohortPair.h"
#include "Environment.h"


/** \brief Merges cohorts with similar properties */
class CohortMerger {
public:
    CohortMerger( ){;}

    /** \brief Calculate the distance between two cohorts in multi-dimensional trait space (body mass, adult mass, juvenile mass)
    @param Cohort1 The first cohort to calculate distance to 
    @param Cohort2 The cohort to compare to 
    @returns The relative distance in trait space */
    static double CalculateDistance( Cohort*, Cohort* );

    /** \brief Merge cohorts until below a specified threshold number of cohorts in each grid cell
    @param gridCellCohorts The cohorts within this grid cell 
    @param TotalNumberOfCohorts The total number of cohorts in this grid cell 
    @param TargetCohortThreshold The target threshold to reduce the number of cohorts to 
    @return The number of cohorts that have been merged */
    static int MergeToReachThresholdFast(vector<Cohort*> ) ;

};

#endif
