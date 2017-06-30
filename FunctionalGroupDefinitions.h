#ifndef FUNCTIONALGROUPDEFINITIONS
#define FUNCTIONALGROUPDEFINITIONS
#include<string>
#include<map>

/** \brief Reads in and performs look-ups on functional group definitions */
class FunctionalGroupDefinitions {

public:
std::string Trait(unsigned,std::string);
double Property(unsigned,std::string);
~FunctionalGroupDefinitions( );

unsigned size();
//This is a base class for a singleton so has a protected constructor
protected:
    /** \brief Setup the functional group definitions: reads in the specified functional group definition file, 
      Constructs lookup tables, mass ranges and initial cohort numbers in each functional group
      @param fileName The name of the functional group definition file to be read in */
    FunctionalGroupDefinitions(){}
    FunctionalGroupDefinitions(std::string );

    std::map<unsigned, std::map< std::string,std::string> > _Traits;
    std::map<unsigned, std::map< std::string,double> >      _Properties;
    unsigned count;
   
};
#endif
