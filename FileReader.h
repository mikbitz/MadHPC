#ifndef FILEREADER
#define	FILEREADER

#include "Types.h"
#include "Constants.h"
#include "repast_hpc/Properties.h"
class FileReader {
public:
    FileReader( );
    ~FileReader( );

    bool ReadFiles(repast::Properties&);

private:
    bool ReadTextFile( const std::string& );
    
    bool ReadInputDataFiles(repast::Properties& );
    
    void ClearMetadata( );
    
    std::string mFilePath;
    Types::StringVector mMetadataHeadings;
    Types::StringMatrix mMetadata;
    
    Types::IntegerVector mAllocatedCellIndices;

    bool _verbose;
};

#endif

