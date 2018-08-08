/*
 *
 * AutotrophProcessor.h
 *
 *  Created on: May 7, 2018
 *      Author: Mike Bithell
 * 
 *  Derived from Origianl C# code by
 * Emergent Global Patterns of Ecosystem Structure and Function from a Mechanistic General Ecosystem Model , Michael B. J. Harfoot , Tim Newbold, Derek P. Tittensor,  Stephen Emmott, Jon Hutton, Vassily Lyutsarev, Matthew J. Smith, Jörn P. W. Scharlemann, Drew W. Purves PLOS Biology 12(4) 2014 e1001841 https://doi.org/10.1371/journal.pbio.1001841
 */
#include "AutotrophProcessor.h"

AutotrophProcessor::AutotrophProcessor( ) {
    mPhytoplanktonConversionRatio = 10;
    mMsqToKmSqConversion = 1000000;
}

double AutotrophProcessor::ConvertNPPToAutotroph( Environment* LocalEnvironment) {
    // Get NPP from the cell environment
    double NPP = LocalEnvironment->NPP();
    double MissingValue=-9999;
    // If NPP is a missing value then set to zero
    if( NPP == MissingValue ) NPP = 0.0; else{

        // Check that the units of oceanic NPP are gC per m2 per day
        //assert(params.Units["OceanNPP"] == "gC/m2/day" && "Oceanic NPP data are not in the correct units for this formulation of the model");

        //Convert to g/cell/month
        NPP *= mMsqToKmSqConversion;

        //Multiply by cell area to get g/cell/day
        NPP *= LocalEnvironment->Area();

        //Convert to g wet matter, assuming carbon content of phytoplankton is 10% of wet matter
        NPP *= mPhytoplanktonConversionRatio;

    }
    return NPP;

}
