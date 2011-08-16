#ifndef OPTIMIZATIONRESULTS_H
#define OPTIMIZATIONRESULTS_H

#include "Includes/LevenbergMarquardt.h"
#include "Includes/Cell_simulator.h"

class OptimizationResults
{
public:
    SimParameters OptimalParameters()const;
    SimParameters InitialParameters()const;

    Experiment Data()const;
    Experiment FittedData()const;

    std::size_t numEval()const;
    std::size_t numIter()const;

    double SS()const;

    friend class Cell_simulator;

private:
    SimParameters OptimalParameters_;
    SimParameters InitialParameters_;

    Experiment Data_;
    Experiment FittedData_;

    std::size_t numEval_;
    std::size_t numIter_;

    double SS_;

};
std::ostream& operator<<(std::ostream&, const OptimizationResults& o);

#endif // OPTIMIZATIONRESULTS_H
