#ifndef FILEREADER
#define	FILEREADER

#include "Types.h"
#include "Constants.h"

/**
\file Filereader.h
\brief This is the Netcdf-file reading class header file
*/
/**
\brief Read in the driving datafiles from netcdf 
 */

class FileReader {
public:
/**
* \brief Constructor - call from main.
* 
* Simply allows for the setting of verbose output
* @param v boolean variable, true for verbose output.
*/

    FileReader(bool);
/**
\brief Destructor - the class allocates no permanent internal storage, so little to do
*/
//------------------------------------------------------------------------------------------------------------
    ~FileReader( );
/**
* \brief Call this routine to read the entire input file list. Called from main.
* 
* Calls ReadFileNames in order to gather the list of files to be read, then calls ReadInputDatafiles
* to open the list one at a time, read in data and store it in the DataLayerSet singleton.
 */
//------------------------------------------------------------------------------------------------------------
    bool ReadFiles();
//------------------------------------------------------------------------------------------------------------
private:
//------------------------------------------------------------------------------------------------------------
/**
*\brief Read a (comma) delimited text file that describes the files to be read.
*
*There is assumed to be a single header line, follwed by one line to describe each file.
*The resulting data is store in _Headers and the _FileDescriptor vector
*@param filePath The full path to the file
*/
    bool ReadFileNames( const std::string& );    
//------------------------------------------------------------------------------------------------------------
/**
* \brief  Reads data files. Creates storage for data and hands this off data to a DataLayerSet
*
* 
*/
    bool ReadInputDataFiles();  
//------------------------------------------------------------------------------------------------------------
/**
\brief  Get a one-dimensional vector from the netcdf file corresponding to one of the data co-ordinates
*/    
    std::vector<float> GetCoordinateVector(netCDF::NcFile& , 
                                           std::string );
//------------------------------------------------------------------------------------------------------------
/**
\brief  Get the actual data from the netcdf file, assumed to be 2D, possibly changing with time
*/
    float*             GetDataField       (netCDF::NcFile& , 
                                           std::string );
//------------------------------------------------------------------------------------------------------------
/** 
\brief store headings from the description file 
 */
    std::vector<std::string> _Headings;
/** 
\brief store the descpription for a single file that needs to be read
*/
    std::map<std::string,std::string>_DataDescriptor;
/** 
\brief store the set of DataDescrptors, one for each file to be read
*/
    std::vector<std::map<std::string,std::string>>_FileDescriptor;
/**
\brief store whether output will be verbose - e.g. write out full file paths to stdout
 */    
    bool _verbose;
};

#endif

