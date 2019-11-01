#ifndef DATALAYERSET
#define	DATALAYERSET

#include "Types.h"


class DataLayerSet {
public:
    ~DataLayerSet( );
    static Types::DataLayerSetPointer Get( );
    
    Types::DataLayerPointer GetDataLayerWithName( const std::string& );
    Types::VariablePointer GetDefaultVariableFor( const std::string& );
    void SetDataLayers( const Types::InputDataPointer );
    
    float GetDataAtCellIndexFor( const std::string, const unsigned );
    float GetDataAtCellLonLatFor( const std::string, const double, const double );
private:
    DataLayerSet( );
    
    static Types::DataLayerSetPointer mThis;
    
    Types::DataLayerMap mDataLayerMap;
};

#endif

