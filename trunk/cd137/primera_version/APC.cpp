#include "APC.h"
/// Main step for APC


double APC_cells::num() const
    {
        return num_Ag_d+num_free_d+num_LT_bound_d+num_exhausted_d;
    }

double APC_cells::IFNgamma_production_rate() const
    {
        double sum=IFN_free_prod_rate_d*num_free_d+
                   IFN_Ag_prod_rate_d*num_Ag_d+
                   IFN_bound_prod_rate_d*num_LT_bound_d;

        return sum;
    }

double APC_cells::TNF_production_rate() const
    {
        double sum=TNF_free_prod_rate_d*num_free_d+
                   TNF_Ag_prod_rate_d*num_Ag_d+
                   TNF_bound_prod_rate_d*num_LT_bound_d;

        return sum;
    }


double& APC_cells::num_Ag()
    {
        return num_Ag_d;
    }

const double& APC_cells::num_Ag()const
    {
        return num_Ag_d;
    }

double APC_cells::num_bound() const
    {
        return num_LT_bound_d;
    }

double& APC_cells::num_exhausted ()
    {
        return num_exhausted_d;
    }

const double& APC_cells::num_exhausted()const
    {
        return num_exhausted_d;
    }

const double& APC_cells::no_to_free_rate_per_Ag()const
    {
        return APC_no_to_free_rate_per_Ag_d;
    }

double& APC_cells::no_to_free_rate_per_Ag()
    {
        return APC_no_to_free_rate_per_Ag_d;
    }
