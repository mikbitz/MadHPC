/*
 *
 * Environment.h
 *
 *  Created on: May 7, 2018
 *      Author: Mike Bithell
 * 
 *  Derived from Origianl C# code by
 * Emergent Global Patterns of Ecosystem Structure and Function from a Mechanistic General Ecosystem Model , Michael B. J. Harfoot , Tim Newbold, Derek P. Tittensor,  Stephen Emmott, Jon Hutton, Vassily Lyutsarev, Matthew J. Smith, Jörn P. W. Scharlemann, Drew W. Purves PLOS Biology 12(4) 2014 e1001841 https://doi.org/10.1371/journal.pbio.1001841
 */
#ifndef ENVIRONMENT
#define ENVIRONMENT
#include "EnvironmentCell.h"
#include "Constants.h"
#include "Types.h"
#include "TimeStep.h"
#include <vector>
class Environment {
public:
	Environment();
    Environment(int,int,int,int);
	~Environment() {}
	std::map<int, EnvironmentCell*> operator[](unsigned i){return _Cells[i];}
    void update(unsigned );
private:

    
 std::map< int, std::map<int, EnvironmentCell*> >  _Cells;
};
#endif


