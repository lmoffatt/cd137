#ifndef CELL_SIMULATOR_H_INCLUDED
#define CELL_SIMULATOR_H_INCLUDED
#include <iostream>
#include <string>
#include "Includes/SimParameters.h"
#include "Includes/Media.h"
#include "Includes/APC.h"
#include "Includes/NK.h"
#include "Includes/LT.h"
#include "Results.h"

class Cell_simulator
{
    public:

        void ask_parameters();
        Cell_simulator(const SimParameters& sp);
        void run();

        Results Simulate(const SimParameters& simPar,
                         const Results& results);

        void update(double time_step);


        Cell_simulator(){}

    private:
        Media   m;
        APC_cells APC;
        NK_cells NK;
        LT_cells LT;

    double time_step_d;
    double sim_duration_d;
    double trun_d;
    std::string filename;


};




#endif // CELL_SIMULATOR_H_INCLUDE
