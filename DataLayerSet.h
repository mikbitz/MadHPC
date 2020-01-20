#ifndef DATALAYERSET_H
#define	DATALAYERSET_H

#include "Layers.h"


class DataLayerSet {

public:
//------------------------------------------------------------------------------------------------------------
//Remove all layers - thsi destroys all associated storage
    ~DataLayerSet( );
//------------------------------------------------------------------------------------------------------------
//Function to return the class pointer to itself 
    static DataLayerSet* Data( );
//------------------------------------------------------------------------------------------------------------
//return just one float value from a named layer at a given Longitude and Latitude  
    float GetDataAtLonLatFor( const std::string, 
                              const double, 
                              const double );
//------------------------------------------------------------------------------------------------------------
//add a time indepedent longitude-latitude layer
    void addLayer(std::string ,
                  float *,
                  std::vector<float>,
                  std::vector<float>);
//------------------------------------------------------------------------------------------------------------
//add a longitude-latitude layer with time dependence
    void addLayer(std::string ,
                  float *,
                  std::vector<float>,
                  std::vector<float>,
                  std::vector<float>);
//------------------------------------------------------------------------------------------------------------
// add a pre-existing Layer
    void addLayer(std::string ,
                  Layer *);
//------------------------------------------------------------------------------------------------------------
// return a pointer to a pre-existing layer
    Layer* GetLayer( const std::string name);
private:
    
//private constructor - clients get access to this class by calling Data() above
    DataLayerSet( );
//------------------------------------------------------------------------------------------------------------
//The class is a singleton so needs a pointer to itself
    static DataLayerSet* _DataSet;
//------------------------------------------------------------------------------------------------------------
//named layers with strings as keys
    std::map<std::string,Layer*> Layers;
};

#endif

