#include "ResultsSimulator.h"
#include <vector>
Experiment ResultsSimulator::Simulate(const SimParameters& simPar,
                 const Experiment& E,
                                      bool isDirectInteraction)
{
    Cell_simulator sim(simPar,E, isDirectInteraction);
    Experiment out;
    for (std::size_t i=0; i<E.size();i++)

    {
        Results simRes=sim.Simulate(simPar,E.Treatment_i(i),E.Result_i(i));
        out.push_back(E.Treatment_i(i),simRes);
    };
    return out;


}


