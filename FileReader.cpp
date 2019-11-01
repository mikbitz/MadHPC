#include "FileReader.h"

#include "Convertor.h"
#include "Constants.h"
#include "InputData.h"
#include "DataLayerSet.h"
#include "Parameters.h"
#include "DataCoords.h"

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

            Types::InputDataPointer initialInputData = new InputData( );

            for( unsigned environmentalDataFileIndex = 0; environmentalDataFileIndex < mMetadata.size( ); ++environmentalDataFileIndex ) {

                std::string filePath = Parameters::Get( )->GetRootDataDirectory( );
                filePath.append( Convertor::Get( )->ToString( Parameters::Get( )->GetDataGridCellSize( ) ) );
                filePath.append( "deg/" );
                filePath.append( mMetadata[ environmentalDataFileIndex ][ Constants::eFilePath ] );
                if (_verbose)std::cout << "Reading NetCDF file \"" << filePath << "\"..." << std::endl;

                try {
                    netCDF::NcFile inputNcFile( filePath.c_str(), netCDF::NcFile::read ); // Open the file for read access

                    std::multimap< std::string, netCDF::NcVar > multiMap = inputNcFile.getVars( );

                    // Outer variable loop
                    for( std::multimap<std::string, netCDF::NcVar>::iterator it = multiMap.begin( ); it != multiMap.end( ); ++it ) {
                        std::string variableName = ( *it ).first;
                        netCDF::NcVar variableNcVar = ( *it ).second;
                        std::vector< netCDF::NcDim > varDims = variableNcVar.getDims( );

                        Types::UnsignedVector variableDimensions;
                        unsigned variableSize = 1;

                        // Inner variable dimension loop
                        for( unsigned dimIndex = 0; dimIndex < varDims.size( ); ++dimIndex ) {
                            variableDimensions.push_back( varDims[ dimIndex ].getSize( ) );
                            variableSize *= varDims[ dimIndex ].getSize( );
                        }

                        float* variableData = new float[ variableSize ];
                        variableNcVar.getVar( variableData );
                        bool isDefault = variableName == Convertor::Get( )->ToLowercase( mMetadata[ environmentalDataFileIndex ][ Constants::eDefaultVariableName ] );
                        initialInputData->AddVariableToDatum( mMetadata[ environmentalDataFileIndex ][ Constants::eInternalName ], variableName, variableDimensions, variableSize, variableData, isDefault );
                    }

                } catch( netCDF::exceptions::NcException& e ) {
                    std::cout<<e.what( )<<std::endl;
                    std::cout << "ERROR> File path \"" << filePath << "\" is invalid." << std::endl;
                }
            }

            DataLayerSet::Get( )->SetDataLayers( initialInputData );
            success = true;
        } else {
            success = false;
        }
    }

    return success;
}

void FileReader::ClearMetadata( ) {
    for( unsigned rowIndex = 0; rowIndex < mMetadata.size( ); ++rowIndex ) {
        mMetadata[ rowIndex ].clear( );
    }
    mMetadata.clear( );
    mMetadataHeadings.clear( );
}
