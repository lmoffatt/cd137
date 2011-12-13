#include "LevenbergMarquardt_Test.h"
#include "MatrixInverse.h"
#include <cmath>
#include <iostream>
std::vector<double> sampleFunction::yfit(const std::vector<double>& parameters)
{
    const std::size_t numSamples=100;
    double dt=0.02;
    std::vector<double> result(numSamples);

    for (std::size_t i=0; i<result.size();++i)
    {
        double t=dt*i;
        result[i]=parameters[0]+t*parameters[1]+t*t*parameters[2]+
                t*t*t*parameters[3];
    }

    return result;

}

void LevenbergMarquardt_Test(){
    std::cerr<<"aqui\n";

    sampleFunction f;
    std::cerr<<"aqui\n";

    std::size_t numPar=4;

    std::vector<double> originalParam(numPar);
    originalParam[0]=4;
    originalParam[1]=2.5;
    originalParam[2]=6.3;
    originalParam[3]=2;
    std::vector<double> data=f.yfit(originalParam);

    std::vector<double> initialParam(numPar);

    initialParam[0]=47.01;
    initialParam[1]=51;
    initialParam[2]=1;
    initialParam[3]=0.3;
    LevenbergMarquardt LM(&f,data,initialParam);

    LM.optimize();

    std::cout<<originalParam;
    std::cout<<LM.OptimParameters();

}

