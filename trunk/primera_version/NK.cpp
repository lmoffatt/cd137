#include "NK.h"

double NK_cells::num() const
    {
        return NK_num_Ag_d+NK_num_free_d+NK_num_LT_bound_d+NK_num_exhausted_d;
    }

double NK_cells::IFNgamma_production_rate() const
    {
        double sum=NK_IFN_free_prod_rate_d*NK_num_free_d+
                   NK_IFN_Ag_prod_rate_d*NK_num_Ag_d+
                   NK_IFN_bound_prod_rate_d*NK_num_LT_bound_d+
                   NK_num_exhausted_d*0;
        return sum;
    }

double NK_cells::TNF_production_rate() const
    {
        double sum=NK_TNF_free_prod_rate_d*NK_num_free_d+
                   NK_TNF_Ag_prod_rate_d*NK_num_Ag_d+
                   NK_TNF_bound_prod_rate_d*NK_num_LT_bound_d+
                   NK_num_exhausted_d*0;
        return sum;
    }


double& NK_cells::NK_num_Ag()
    {
        return NK_num_Ag_d;
    }

const double& NK_cells::NK_num_Ag()const
    {
        return NK_num_Ag_d;
    }

double& NK_cells::NK_num_bound()
    {
        return NK_num_LT_bound_d;
    }

const double& NK_cells::NK_num_bound() const
    {
        return NK_num_LT_bound_d;
    }

double& NK_cells::NK_exhausted ()
    {
        return NK_num_exhausted_d;
    }

const double& NK_cells::NK_exhausted()const
    {
        return NK_num_exhausted_d;
    }

const double& NK_cells::NK_no_to_free_rate_per_Ag()const
    {
        return NK_no_to_free_rate_per_Ag_d;
    }

double& NK_cells::NK_no_to_free_rate_per_Ag()
    {
        return NK_no_to_free_rate_per_Ag_d;
    }
