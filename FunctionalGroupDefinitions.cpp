#include "FunctionalGroupDefinitions.h"
#include "assert.h"
#include "Constants.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include <algorithm>


//--------------------------------------------------------------------------
FunctionalGroupDefinitions::FunctionalGroupDefinitions( std::string fileName ) {
    //std::cout << "Reading \"" << fileName << "\" functional group definitions." << std::endl;
    
    std::ifstream infile( fileName.c_str( ) );

    if( infile.is_open( ) ) {
        std::string l;
        std::vector<std::string> header, category;
        getline( infile, l );
        //trim off newline character
        l.pop_back( );
        std::istringstream s( l );
        //split out the comma-separated and underscore separated header
        while( s.good( ) ) {
            std::string tmp;
            getline( s, tmp, ',' );
            std::transform( tmp.begin( ), tmp.end( ), tmp.begin( ), ::tolower );
            //split out the header category (definition,property or note)
            std::istringstream splt( tmp );
            std::string dp, op;
            getline( splt, dp, '_' );
            category.push_back( dp );
            getline( splt, op, '_' );
            header.push_back( op );
        }
        count = 0;
        //retrieve the lines defining each functional group
        while( infile.good( ) ) {
            std::string l, data;
            getline( infile, l );
            if( infile.good( ) ) {
                l.pop_back( );
                if( l.length( ) > 1 ) {
                    std::istringstream s( l );
                    //step through the columns for this functional group
                    for( unsigned i = 0; i < header.size( ); i++ ) {
                        getline( s, data, ',' );
                        std::transform( data.begin( ), data.end( ), data.begin( ), ::tolower );

                        if( category[i] == "definition" ) {
                            //for each trait, store the value for a given functional group
                            //indexed by functional group number
                            _Traits[ count][header[ i ] ]= data ;
                        }
                        //Otherwise get the value for the given property
                        //for this functional group
                        if( category[i] == "property" ) {
                            _Properties[count][ header[ i ] ]=atof( data.c_str( ) ) ;
                        }
                    }
                }
                count++;
            }
        }
    } else {
        std::cout << "Something wrong with functional group definitions file " << fileName << std::endl;
    }
    infile.close( );
}
//--------------------------------------------------------------------------
std::string FunctionalGroupDefinitions::Trait(unsigned i,std::string s){std::transform(s.begin( ), s.end( ), s.begin( ), ::tolower );return _Traits[i][s];}
//--------------------------------------------------------------------------
double FunctionalGroupDefinitions::Property(unsigned i,std::string s){std::transform( s.begin( ), s.end( ), s.begin( ), ::tolower );return _Properties[i][s];}
//--------------------------------------------------------------------------
unsigned FunctionalGroupDefinitions::size() {
    return count;
}
