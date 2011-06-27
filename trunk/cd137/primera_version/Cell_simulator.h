#ifndef CELL_SIMULATOR_H_INCLUDED
#define CELL_SIMULATOR_H_INCLUDED
#include <iostream>
#include <string>
#include "SimParameters.h"
#include "Media.h"
#include "APC.h"
#include "NK.h"
#include "LT.h"

class Cell_simulator
{
    public:

        void ask_parameters();
        Cell_simulator(const SimParameters& sp);
        void run();
        void update(double time_step);

        Cell_simulator(){};

    private:
        LT_cells LT;
        APC_cells APC;
        Media   m;
        NK_cells NK;

    double time_step_d;
    double trun_d;
    double sim_duration_d;
    double max_num_cells_;
    double init_num_APC_cells;
    double init_num_NK_cells;
    double init_num_LT_cells;
    double LT_num_specific;
    double Ag;
    double Ag_internalization_rate;
    double APC_max_proliferation_rate_;
    double NK_max_proliferation_rate_;
    double LT_max_no_receptor_prol_rate_;
    double LT_max_free_prol_rate_;
    double LT_max_bound_prol_rate_;
    double APC_no_to_free_rate_per_Ag_;
    double APC_free_to_bound_rate_per_LT_;
    double APC_exh_rate;
    double NK_no_to_free_rate_per_Ag_;
    double NK_free_to_bound_rate_per_LT_;
    double NK_exh_rate;
    double LT_no_to_free_rate_per_APC_;
    double LT_free_to_bound_rate_per_APC_;
    double LT_IFN_no_rec_prod_rate_;
    double LT_IFN_free_prod_rate_;
    double LT_IFN_bound_prod_rate_;
    double LT_TNF_no_rec_prod_rate_;
    double LT_TNF_free_prod_rate_;
    double LT_TNF_bound_prod_rate_;
    double APC_IFN_free_prod_rate_;
    double APC_IFN_Ag_prod_rate_;
    double APC_IFN_bound_prod_rate_;
    double APC_TNF_free_prod_rate_;
    double APC_TNF_Ag_prod_rate_;
    double APC_TNF_bound_prod_rate_;
    double NK_IFN_free_prod_rate_;
    double NK_IFN_Ag_prod_rate_;
    double NK_IFN_bound_prod_rate_;
    double NK_TNF_free_prod_rate_;
    double NK_TNF_Ag_prod_rate_;
    double NK_TNF_bound_prod_rate_;
    double TNF_deg;

   std::string filename;

};




#endif // CELL_SIMULATOR_H_INCLUDE
