#ifndef LEVENBERGMARQUARDT_TEST_H
#define LEVENBERGMARQUARDT_TEST_H


#include "LevenbergMarquardt.h"

class sampleFunction: public ABC_function
{
public:
    virtual std::vector<double> yfit(const std::vector<double>& parameters);
    virtual std::vector<double> yfit(const Parameters& parameters)const{}


};


void LevenbergMarquardt_Test();

#endif // LEVENBERGMARQUARDT_TEST_H
