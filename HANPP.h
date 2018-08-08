/*
 *
 * HANPP.h
 *
 *  Created on: May 7, 2018
 *      Author: Mike Bithell
 * 
 *  Derived from Origianl C# code by
 * Emergent Global Patterns of Ecosystem Structure and Function from a Mechanistic General Ecosystem Model , Michael B. J. Harfoot , Tim Newbold, Derek P. Tittensor,  Stephen Emmott, Jon Hutton, Vassily Lyutsarev, Matthew J. Smith, Jörn P. W. Scharlemann, Drew W. Purves PLOS Biology 12(4) 2014 e1001841 https://doi.org/10.1371/journal.pbio.1001841
 */
#ifndef HANPP_H
#define HANPP_H

#include "Constants.h"
#include "Environment.h"
/** \file HANPP.h
 * \brief the HANPP header file
 */

/** \brief Removes autotroph matter appropriated by humans from a grid cell's autotroph stocks
\remarks   Assumes that autotroph matter is appropriated evenly from different stocks in proportion to their biomass */
using namespace std;

class HumanAutotrophMatterAppropriation {
public:
    //----------------------------------------------------------------------------------------------
    //Methods
    //----------------------------------------------------------------------------------------------

    /** \brief Constructor for human appropriation of autotroph matter */
    HumanAutotrophMatterAppropriation( );
    //----------------------------------------------------------------------------------------------

    /** \brief Remove human appropriated matter from the grid cell autotroph stocks
    @param gcl The current grid cell 
    @param humanNPPExtraction The type of NPP extraction to apply: 'no' = no removal; 'hanpp' = appropriated NPP estimate from input map; or proportion of total NPP 
    @param gridCellStocks The stocks in the current grid cell 
    @param actingStock The position of the acting stock in the jagged array of grid cell stocks 
    @param currentTimestep The current model time step */
    double RemoveHumanAppropriatedMatter(Environment* , double,double, unsigned);

};
#endif
