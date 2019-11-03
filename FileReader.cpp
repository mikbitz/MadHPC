#include "FileReader.h"

#include "Convertor.h"
#include "Constants.h"
#include "InputData.h"
#include "DataLayerSet.h"
#include "Parameters.h"

#include <netcdf>

FileReader::FileReader( ) {

}

FileReader::~FileReader( ) {

}
bool FileReader::ReadFiles(repast::Properties& props ) {

    bool success = false;

    success = Parameters::Get()->Initialise( props );
    //extra output if requested
    _verbose=(props.getProperty("verbose")=="true");

    if( success == true )
        success = ReadInputDataFiles( props );
    
    return success;
}

bool FileReader::ReadTextFile( const std::string& filePath ) {

    bool success = false;

    ClearMetadata( );

    if(_verbose)std::cout << "Reading text file \"" << filePath << "\"..." << std::endl;
    std::ifstream fileStream( filePath.c_str( ), std::ios::in );

    if( fileStream.is_open( ) ) {
        std::string readLine;
        unsigned lineCount = 0;

        while( std::getline( fileStream, readLine ) ) {
            if( lineCount > 0 && readLine[ 0 ] != Constants::cTextFileCommentCharacter ) {
                //StringToWords returns a string split into a vector using the Delimiter value
                mMetadata.push_back( Convertor::Get( )->StringToWords( readLine, Constants::cDataDelimiterValue ) );
            } else if( lineCount == 0 ) {
                mMetadataHeadings = Convertor::Get( )->StringToWords( readLine, Constants::cDataDelimiterValue );
            }
            ++lineCount;
        }
        success = true;
        fileStream.close( );

    } else {
        std::cout << "File path \"" << filePath << "\" is invalid." << std::endl;
    }

    return success;
}

bool FileReader::ReadInputDataFiles(repast::Properties& props ) {
    //this file is expected to list the internal vaiable name for the model, the file path below the root (including filename extension) and the default variable name in the file to be read in.
    bool success = ReadTextFile(props.getProperty("input.DataDirectory") + "/" + props.getProperty("input.EnvironmentVariablesFile") );

    if( success == true ) {

        if( mMetadata.size( ) > 0 ) {

            for( unsigned environmentalDataFileIndex = 0; environmentalDataFileIndex < mMetadata.size( ); ++environmentalDataFileIndex ) {

                std::string filePath = Parameters::Get( )->GetRootDataDirectory( );
                filePath.append( Convertor::Get( )->ToString( Parameters::Get( )->GetDataGridCellSize( ) ) );
                filePath.append( "deg/" );
                filePath.append( mMetadata[ environmentalDataFileIndex ][ Constants::eFilePath ] );
                if (_verbose)std::cout << "Reading NetCDF file \"" << filePath << "\"..." << std::endl;

                try {
                    netCDF::NcFile inputNcFile( filePath.c_str(), netCDF::NcFile::read ); // Open the file for read access

                    std::multimap< std::string, netCDF::NcVar > multiMap = inputNcFile.getVars( );
                    
                    bool hasTime=(multiMap.find("time")!=multiMap.end());
                    
                    auto name=mMetadata[ environmentalDataFileIndex ][ Constants::eInternalName ];
                    
                    auto dataNc=inputNcFile.getVar(mMetadata[ environmentalDataFileIndex ][ Constants::eDefaultVariableName ]);
                    assert (!dataNc.isNull());
                    auto dims=dataNc.getDims();
                    unsigned dataSize=1;
                    for( unsigned dimIndex = 0; dimIndex < dims.size( ); ++dimIndex ) {
                        dataSize *= dims[ dimIndex ].getSize( );
                    }
                    float* data=new float[dataSize];
                    dataNc.getVar( data );

                    auto lon = Get1DNcFloatVector(inputNcFile, "lon" );
                    auto lat = Get1DNcFloatVector(inputNcFile, "lat" );

                    if(hasTime){
                        auto time= Get1DNcFloatVector(inputNcFile, "time");
                        DataLayerSet::Data( )->addLayer(name,data,lon,lat,time);
                    }else{
                        DataLayerSet::Data( )->addLayer(name,data,lon,lat);
                    }


                } catch( netCDF::exceptions::NcException& e ) {
                    std::cout<<e.what( )<<std::endl;
                    std::cout << "ERROR> File path \"" << filePath << "\" is invalid." << std::endl;
                }
            }

            success = true;

        } else {
            success = false;
        }
    }

    return success;
}
std::vector<float> FileReader::Get1DNcFloatVector(netCDF::NcFile& file, std::string name) {
    auto varNc=file.getVar(name);
    assert(!varNc.isNull());
    std::vector<float>var(varNc.getDims()[0].getSize());
    varNc.getVar(var.data());
    return var;
}

void FileReader::ClearMetadata( ) {
    for( unsigned rowIndex = 0; rowIndex < mMetadata.size( ); ++rowIndex ) {
        mMetadata[ rowIndex ].clear( );
    }
    mMetadata.clear( );
    mMetadataHeadings.clear( );
}
