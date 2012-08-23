#ifndef RESULTSSIMULATOR_H
#define RESULTSSIMULATOR_H
#include "Results.h"
#include "SimParameters.h"
#include "Cell_simulator.h"
#include "Includes/Experiment.h"

class ResultsSimulator {
public:
    Experiment Simulate(const SimParameters& simPar,
                     const Experiment& results,
                        bool isDirectInteraction);

    double SumSquare(const SimParameters& SimPar,
                     const Results& results);



};


#endif // RESULTSSIMULATOR_H
