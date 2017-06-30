#ifndef GROUPS_H
#define GROUPS_H

#include "FunctionalGroupDefinitions.h"
#include "Constants.h"
#include "assert.h"

class StockDefinitions : public FunctionalGroupDefinitions {
public:
//--------------------------------------------------------------------------
static void Initialise( std::string fileName) {
    mThis = new StockDefinitions(fileName );
}
//--------------------------------------------------------------------------
static StockDefinitions* Get( ) {
    assert(mThis!=NULL && "Stock group definitions must be initialised before use!");
    return mThis;
}
private:
static StockDefinitions* mThis;
StockDefinitions( std::string fileName ): FunctionalGroupDefinitions(fileName) {}
//-------------------------------------------------------------------------
~StockDefinitions( ) {
    if( mThis != NULL ) delete mThis;
}
};


//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
class CohortDefinitions : public FunctionalGroupDefinitions {
public:
//--------------------------------------------------------------------------
static void Initialise( std::string fileName) {
    mThis = new CohortDefinitions(fileName );
}
//--------------------------------------------------------------------------
static CohortDefinitions* Get( ) {
    assert(mThis!=NULL && "Cohort group definitions must be initialised before use!");
    return mThis;
}
private:
static CohortDefinitions* mThis ;
CohortDefinitions( std::string fileName ): FunctionalGroupDefinitions(fileName) {}
//--------------------------------------------------------------------------
~CohortDefinitions( ) {
    if( mThis != NULL ) delete mThis;
}
};


#endif
