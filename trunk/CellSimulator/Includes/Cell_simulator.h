#ifndef CELL_SIMULATOR_H_INCLUDED
#define CELL_SIMULATOR_H_INCLUDED
#include <iostream>
#include <string>
#include "Includes/SimParameters.h"
#include "Includes/Media.h"
#include "Includes/APC.h"
#include "Includes/NK.h"
#include "Includes/LT.h"
#include "Includes/OptimizationResults.h"
#include "Includes/LevenbergMarquardt.h"

#include "Experiment.h"
class Cell_simulator: public ABC_function
{
public:

    void ask_parameters();
    void run();

    Cell_simulator& applyParameters(const SimParameters& sp,
				    const Treatment& tr);
    Cell_simulator(const SimParameters& sp,
                   const Experiment& E);

    Results Simulate(const SimParameters& simPar,
                     const Treatment& protocol,
                     const Results& results);

    Experiment Simulate(const SimParameters& simPar,
                        const Experiment& exp);


    OptimizationResults Optimize(const SimParameters& simPar,
				 const Experiment& exp);


    void update(double time_step);



   virtual std::vector<double> yfit (const std::vector<double>& param);






    Cell_simulator(const Cell_simulator& other);

    friend void swap(Cell_simulator& one, Cell_simulator& other);

    Cell_simulator& operator=(const Cell_simulator& other);

    Cell_simulator();

    void reset(const SimParameters& sp,const Treatment& tr);

private:
    Media   m;
    APC_cells APC;
    NK_cells NK;
    LT_cells LT;

    double time_step_d;
    double sim_duration_d;
    double trun_d;
    std::string filename;


    Experiment experiment_;
    Experiment fitExperiment_;
    SimParameters initialPar_;
    SimParameters fitPar_;




    void formatInputsforLM();
    void runLM();
    void formatOutputsfromLM();
    OptimizationResults buildOptimizationResults();

};




#endif // CELL_SIMULATOR_H_INCLUDE
