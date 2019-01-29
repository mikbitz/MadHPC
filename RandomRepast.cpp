#include "RandomRepast.h"

RandomRepast::RandomRepast(){}
RandomRepast::~RandomRepast(){}
    void RandomRepast::SetSeed( const unsigned ){  
        //initializeRandom(props, &world);
    }

    void RandomRepast::SetSeedFromSystemTime( ){
        //corresponds to global.random.seed=AUTO in model.props
    };


    double RandomRepast::GetUniform( ){
        repast::DoubleUniformGenerator gen = repast::Random::instance()->createUniDoubleGenerator(0, 1);
        return gen.next( );
    }

    unsigned RandomRepast::GetUint( ){
        repast::IntUniformGenerator  UI=repast::Random::instance()->createUniIntGenerator(0,U_MAX);
        return UI.next();
    }


    double RandomRepast::GetNormal( const double& mean , const double& standardDeviation  ){
        repast::NormalGenerator NJ=repast::Random::instance()->createNormalGenerator(mean,standardDeviation);
        return NJ.next();
    }


    double RandomRepast::GetLogNormal( const double& mean, const double& standardDeviation){
        repast::LogNormalGenerator LNJ= repast::Random::instance()->createLogNormalGenerator(exp(mean), exp(standardDeviation));
        return LNJ.next();
    }
