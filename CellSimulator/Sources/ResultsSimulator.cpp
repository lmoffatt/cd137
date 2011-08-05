#include "ResultsSimulator.h"

Results ResultsSimulator::Simulate(const SimParameters& simPar,
                                   const Results& results)
{
    Cell_simulator Cell(simPar);
    return Cell.Simulate(simPar,results);


}


