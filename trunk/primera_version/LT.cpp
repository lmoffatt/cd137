#include "LT.h"


double LT_cells::num() const
    {
        return num_Agsp_bound_receptor_d+num_Agsp_free_receptor_d+num_Agsp_no_receptor_d+num_non_Agsp_d;
    }

double LT_cells::IFNgamma_production_rate() const
    {
        return (num_non_Agsp_d+num_Agsp_no_receptor_d)*IFN_no_rec_prod_rate_d+
                num_Agsp_free_receptor_d*IFN_free_prod_rate_d+
                num_Agsp_bound_receptor_d*IFN_bound_prod_rate_d;
    };

double LT_cells::TNF_production_rate() const
    {
        return (num_non_Agsp_d+num_Agsp_no_receptor_d)*TNF_no_rec_prod_rate_d+
                num_Agsp_free_receptor_d*TNF_free_prod_rate_d+
                num_Agsp_bound_receptor_d*TNF_bound_prod_rate_d;
    };

double LT_cells::num_cells_not_Ag_specific()const
    {
        return num_non_Agsp_d;
    };

double LT_cells::num_cells_not_expressing_receptor()const
    {
        return num_Agsp_no_receptor_d;
    };

double LT_cells::num_cells_expressing_receptor_and_free()const
    {
        return num_Agsp_free_receptor_d;
    };

double LT_cells::num_cells_expressing_receptor()const
    {
        return num_Agsp_free_receptor_d+num_Agsp_bound_receptor_d;
    };

double LT_cells::num_cells_expressing_receptor_and_bound()const
    {
        return num_Agsp_bound_receptor_d;
    };
