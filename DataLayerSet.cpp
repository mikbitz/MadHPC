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
    Layer* d=new Layer2D(data,lon,lat,Parameters::instance()->GetSpatialInterpolation());
    Layers[name]=d;
}
//------------------------------------------------------------------------------------------------------------
void DataLayerSet::addLayer(std::string name,float *data,std::vector<float>lon,std::vector<float>lat,std::vector<float>time){
    Layer* d=new Layer2DWithTime(data,lon,lat,time,Parameters::instance()->GetSpatialInterpolation());
    Layers[name]=d;
}
//------------------------------------------------------------------------------------------------------------
void DataLayerSet::addLayer(std::string name, Layer* d){
    Layers[name]=d;
}
//------------------------------------------------------------------------------------------------------------
Layer* DataLayerSet::GetLayer( const std::string name) {
    
    assert(Layers.find( name ) != Layers.end( ));
    
    return Layers[name];
}
//------------------------------------------------------------------------------------------------------------
float DataLayerSet::GetDataAtLonLatFor( const std::string name, const double Longitude, const double Latitude ) {
    
    assert(Layers.find( name ) != Layers.end( ));
    
    return Layers[name]->GetDataAtLonLat(Longitude,Latitude);
}
