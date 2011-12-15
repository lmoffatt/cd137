#ifndef CELL_SIMULATOR_H_INCLUDED
#define CELL_SIMULATOR_H_INCLUDED
#include <iostream>
#include <string>
#include "Includes/SimParameters.h"
#include "Includes/Parameters.h"

#include "Includes/Media.h"
#include "Includes/APC.h"
#include "Includes/NK.h"
#include "Includes/LT.h"
#include "Includes/LevenbergMarquardt.h"
#include "Includes/BayesIteration.h"


#include "Experiment.h"

class OptimizationResults;
class Cell_simulator: public ABC_model
{
public:
    ~Cell_simulator(){}
//    void ask_parameters();
    void run();



    Cell_simulator& applyParameters(const SimParameters& sp,
				    const Treatment& tr);


    virtual Cell_simulator& setData(const ABC_data& conditions){};



    Cell_simulator(const SimParameters& sp,
                   const Experiment& E);

    Results Simulate(const SimParameters& simPar,
                     const Treatment& protocol,
                     const Results& results);

    Experiment Simulate(const SimParameters& simPar,
                        const Experiment& exp);


    OptimizationResults Optimize(const SimParameters& priorPar,
                                 const SimParameters& simPar,
				 const Experiment& exp,
				 double range,
				 std::size_t numStarts);



    // versions with Parameters instead of SimParameters

    Cell_simulator& applyParameters(const Parameters& sp,
                                    const Treatment& tr);





    Cell_simulator(const Parameters& sp,
                   const Experiment& E);

    Results Simulate(const Parameters& simPar,
                     const Treatment& protocol,
                     const Results& results);

    Experiment Simulate(const Parameters& simPar,
                        const Experiment& exp);


    OptimizationResults Optimize(const Parameters& priorPar,
                                 const Parameters& simPar,
                                 const Experiment& exp,
                                 double range,
                                 std::size_t numStarts);
















    void update(double time_step);



   virtual std::vector<double> yfit (const std::vector<double>& param);

    virtual std::vector<double> yfit (const Parameters& param);


    std::vector<double> difParam(const std::vector<double>& param);

    std::vector<double> getData(const std::vector<double>& param);



    Cell_simulator(const Cell_simulator& other);

    friend void swap(Cell_simulator& one, Cell_simulator& other);

    Cell_simulator& operator=(const Cell_simulator& other);

    Cell_simulator();

  //  void reset(const SimParameters& sp,const Treatment& tr);

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

    Parameters prior_;
    Parameters current_;




    void formatInputsforLM();
    void runLM();
    void formatOutputsfromLM();

};




#endif // CELL_SIMULATOR_H_INCLUDE
