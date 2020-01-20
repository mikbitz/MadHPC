/*
 *
 * Environment.cpp
 *
 *  Created on: May 7, 2018
 *      Author: Mike Bithell
 * 
 *  Derived from Original C# code by
 * Emergent Global Patterns of Ecosystem Structure and Function from a Mechanistic General Ecosystem Model , Michael B. J. Harfoot , Tim Newbold, Derek P. Tittensor,  Stephen Emmott, Jon Hutton, Vassily Lyutsarev, Matthew J. Smith, Jörn P. W. Scharlemann, Drew W. Purves PLOS Biology 12(4) 2014 e1001841 https://doi.org/10.1371/journal.pbio.1001841
 */
#include "Environment.h"
#include "CohortMerger.h"

#include "ClimateVariablesCalculator.h"

#include "TimeStep.h"
#include "Parameters.h"
#include "UtilityFunctions.h"
#include "Stock.h"
#include "DataLayerSet.h"


Environment::Environment() {}
Environment::Environment(int minX,int maxX,int minY,int maxY) {
    for(int x=minX;x<=maxX;x++) {
        for(int y=minY;y<=maxY;y++) {
            EnvironmentCell* E=new EnvironmentCell(x,y);
            _Cells[x][y]=E;
        }
    }
}

void Environment::update(unsigned CurrentTimeStep){
  TimeStep::instance()->SetTime(CurrentTimeStep);
}


