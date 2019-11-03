#include "DataLayerSet.h"

DataLayerSet* DataLayerSet::_DataSet = NULL;
//------------------------------------------------------------------------------------------------------------
DataLayerSet* DataLayerSet::Data( ) {
    if( _DataSet == NULL ) {
        _DataSet = new DataLayerSet( );
    }
    return _DataSet;
}
//------------------------------------------------------------------------------------------------------------
DataLayerSet::~DataLayerSet( ) {

    for( auto iter :Layers ) {
        delete iter.second;
    }
}
//------------------------------------------------------------------------------------------------------------
DataLayerSet::DataLayerSet( ) {

}
//------------------------------------------------------------------------------------------------------------
void DataLayerSet::addLayer(std::string name,float *data,std::vector<float>lon,std::vector<float>lat){
    Layer* d=new Layer2D(data,lon,lat);
    Layers[name]=d;
}
//------------------------------------------------------------------------------------------------------------
void DataLayerSet::addLayer(std::string name,float *data,std::vector<float>lon,std::vector<float>lat,std::vector<float>time){
    Layer* d=new Layer2DWithTime(data,lon,lat,time);
    Layers[name]=d;
}
//------------------------------------------------------------------------------------------------------------
float DataLayerSet::GetDataAtLonLatFor( const std::string name, const double Longitude, const double Latitude ) {
    
    assert(Layers.find( name ) != Layers.end( ));
    
    return Layers[name]->GetDataAtLonLat(Longitude,Latitude);
}
