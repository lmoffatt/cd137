#ifndef CELL_SIMULATOR_H_INCLUDED
#define CELL_SIMULATOR_H_INCLUDED
#include <iostream>
#include <string>
#include <vector>

#include "Includes/SimParameters.h"
#include "Includes/Parameters.h"

#include "Includes/Media.h"
#include "Includes/APC.h"
#include "Includes/NK.h"
#include "Includes/LT.h"
#include "Includes/LevenbergMarquardt.h"
#include "Includes/BayesIteration.h"
#include "RungeKutta4.h"

#include "Experiment.h"

class OptimizationResults;
class Cell_simulator: public ABC_model, public ABC_ODE
{
public:
    ~Cell_simulator(){}
//    void ask_parameters();

    void run();


    std::ostream& run(std::ostream& f);
    std::ostream& run(std::ostream& s, const Parameters& par)const;


    Cell_simulator& applyParameters(const SimParameters& sp,
                    const Treatment& tr);


    virtual Cell_simulator& setData(const ABC_data& conditions){}



    Cell_simulator(const SimParameters& sp,
                   const Experiment& E, bool isDirectInteraction);

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

    Cell_simulator& applyParameters(const Parameters& current,
                                    const Treatment& tr);



    static Parameters getStandardParameters();
    static Parameters getMinimalParameters();



    Cell_simulator(const Parameters& prior,
                   const Parameters& current,
                   const Experiment& E,
                   bool isDirectInteraction);

    Results Simulate(const Parameters& current,
                     const Treatment& protocol,
                     const Results& results);

    Experiment Simulate(const Parameters& current,
                        const Experiment& exp);


    void Optimize(const Parameters& current,
                  const Experiment& exp,
                  const std::string& filename);



    virtual std::ostream& put(std::ostream &s,const Parameters& parameters) const;


    virtual std::vector<double> Derivative(double t, std::vector<double> y);


    std::vector<double> getState() const;

    void setState(const std::vector<double>& y);






    void update(double time_step);



   virtual std::vector<double> yfit (const std::vector<double>& param);

    virtual std::vector<double> yfit (const Parameters& param)const;


    std::vector<double> difParam(const std::vector<double>& param);

    std::vector<double> getData(const std::vector<double>& param);

    bool hasDirectInteraction()const;
    void setDirectInteraction(bool);

    Cell_simulator(const Cell_simulator& other);

    friend void swap(Cell_simulator& one, Cell_simulator& other);

    Cell_simulator& operator=(const Cell_simulator& other);

    Cell_simulator();

  //  void reset(const SimParameters& sp,const Treatment& tr);

private:
    bool directInteraction;
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
