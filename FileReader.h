#ifndef FILEREADER
#define	FILEREADER

#include "Types.h"
#include "Constants.h"

class FileReader {
public:
    FileReader( );
    ~FileReader( );

    bool ReadFiles( );
    bool ReadInputParameters( );
private:
    bool ReadTextFile( const std::string& );
    

    bool SetUpOutputVariables( );
    bool ReadInputDataFiles( );
    
    void ClearMetadata( );
    
    std::string mFilePath;
    Types::StringVector mMetadataHeadings;
    Types::StringMatrix mMetadata;
    
    Types::IntegerVector mAllocatedCellIndices;
};

#endif

