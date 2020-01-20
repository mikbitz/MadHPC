#include "FileReader.h"

#include "Convertor.h"
#include "Constants.h"

#include "DataLayerSet.h"
#include "Parameters.h"

#include <netcdf>

/**
\file Filereader.cpp
\brief This is the Netcdf-file reading implementation file
*/
//------------------------------------------------------------------------------------------------------------

FileReader::FileReader(bool v ):_verbose(v) {

}
//------------------------------------------------------------------------------------------------------------
FileReader::~FileReader( ) {

}
//------------------------------------------------------------------------------------------------------------
bool FileReader::ReadFiles() {

    //this file is expected to list the internal variable name for the model, the file path below the root (including filename extension) and the default variable name in the file to be read in. These are specified in the EnvironmentVariablesFile header.
    bool success = ReadFileNames(Parameters::instance()->GetDataDecriptorFileName() );
    assert(success);
    success = ReadInputDataFiles( );
    
    return success;
}
//------------------------------------------------------------------------------------------------------------
bool FileReader::ReadInputDataFiles() {
    
    bool success=false;
    
    for(unsigned i=0;i< _FileDescriptor.size();i++) {
        
        std::string filePath = Parameters::instance()->GetRootDataDirectory( );
        
        filePath.append( _FileDescriptor[i]["FilePath"]);
        if (_verbose)std::cout << "Reading NetCDF file \"" << filePath << "\"..." << std::endl;
        
        try {
            netCDF::NcFile inputNcFile( filePath.c_str(), netCDF::NcFile::read ); // Open the file for read access
            
            std::multimap< std::string, netCDF::NcVar > multiMap = inputNcFile.getVars( );
            
            bool hasTime=(multiMap.find("time")!=multiMap.end());
            
            auto name=_FileDescriptor[i]["InternalName"];
            
            auto data = GetDataField       (inputNcFile, _FileDescriptor[i]["DefaultVariableName"]); 
            auto lon  = GetCoordinateVector(inputNcFile, "lon" );
            auto lat  = GetCoordinateVector(inputNcFile, "lat" );
            
            if(hasTime){
                auto time = GetCoordinateVector(inputNcFile, "time");
                DataLayerSet::Data( )->addLayer(name,data,lon,lat,time);
            }else{
                DataLayerSet::Data( )->addLayer(name,data,lon,lat);
            }
            inputNcFile.close();
            
        } catch( netCDF::exceptions::NcException& e ) {
            std::cout<<e.what( )<<std::endl;
            std::cout << "ERROR> File path \"" << filePath << "\" is invalid." << std::endl;
        }
        success = true;
    }
    

    return success;
}
//------------------------------------------------------------------------------------------------------------
std::vector<float> FileReader::GetCoordinateVector(netCDF::NcFile& file, std::string name) {
    auto varNc=file.getVar(name);
    assert(!varNc.isNull());
    std::vector<float>var(varNc.getDims()[0].getSize());
    varNc.getVar(var.data());
    return var;
}
//------------------------------------------------------------------------------------------------------------
float* FileReader::GetDataField(netCDF::NcFile& file, std::string name) {
    auto dataNc=file.getVar(name);
    assert (!dataNc.isNull());
    auto dims=dataNc.getDims();
    unsigned dataSize=1;
    for( unsigned dimIndex = 0; dimIndex < dims.size( ); ++dimIndex ) {
        dataSize *= dims[ dimIndex ].getSize( );
    }
    float* data=new float[dataSize];
    dataNc.getVar( data );
    return data;
}
//------------------------------------------------------------------------------------------------------------
bool FileReader::ReadFileNames( const std::string& filePath ) {

    bool success = false;

    if(_verbose)std::cout << "Reading text file \"" << filePath << "\"..." << std::endl;
    std::ifstream fileStream( filePath.c_str( ), std::ios::in );

    if( fileStream.is_open( ) ) {
        std::string readLine;
        std::getline( fileStream, readLine );
        //StringToWords returns a string split into a vector using the Delimiter value
        _Headings = Convertor::Get( )->StringToWords( readLine, Constants::cDataDelimiterValue );
        assert(_Headings[0]=="InternalName" && _Headings[1]=="FilePath" && _Headings[2]=="DefaultVariableName");
        unsigned lineCount = 0;

        while( std::getline( fileStream, readLine ) ) {
            if( readLine[ 0 ] != Constants::cTextFileCommentCharacter ) {
                auto descriptor = ( Convertor::Get( )->StringToWords( readLine, Constants::cDataDelimiterValue ) );
                _DataDescriptor["InternalName"]       = descriptor[0];
                _DataDescriptor["FilePath"]           = descriptor[1];
                _DataDescriptor["DefaultVariableName"]= descriptor[2];
                _FileDescriptor.push_back(_DataDescriptor);
            }
        }
        success = true;
        fileStream.close( );

    } else {
        std::cout << "File path \"" << filePath << "\" is invalid." << std::endl;
    }
    return success;

}
//------------------------------------------------------------------------------------------------------------

