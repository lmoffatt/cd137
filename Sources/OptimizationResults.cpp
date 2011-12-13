#include "Includes/OptimizationResults.h"

SimParameters OptimizationResults::OptimalParameters()const
{
    return OptimalParameters_;
}
SimParameters OptimizationResults::InitialParameters()const
{
    return InitialParameters_;
}

Experiment OptimizationResults::Data()const
{
    return Data_;
}
Experiment OptimizationResults::FittedData()const
{
    return FittedData_;
}

std::size_t OptimizationResults::numEval()const
{
    return numEval_;
}
std::size_t OptimizationResults::numIter()const
{
    return numIter_;
}

double OptimizationResults::SS()const
{
    return SS_;
}

std::ostream& operator<<(std::ostream& s, const OptimizationResults& o)
{
    s<<"\n Original Data\n";
    s<<o.Data();
    s<<"\n Initial Parameters\n";
    s<<o.InitialParameters();
    s<<"\n Fitted Data \n";
    s<<o.FittedData();
    s<<"\n Fitted Parameters\n";
    s<<o.OptimalParameters();
    s<<"\n SS\n";
    s<<o.SS();
  return s;
}

