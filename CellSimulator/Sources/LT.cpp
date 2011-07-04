#include "Includes/Media.h"
#include "Includes/APC.h"
#include "Includes/LT.h"
#include "Includes/NK.h"
/// Main step for LT
void LT_cells::update(double time_step,const Media& m,const APC_cells& a,const NK_cells& NK)

{
    double proliferation_ratio=(m.Max_num_cells()-m.num_cells())/m.Max_num_cells();


    /// cells not sensitive to the Ag proliferate passively
    num_non_Agsp_d+=time_step*num_non_Agsp_d*proliferation_ratio*LT_max_no_receptor_prol_rate_d;


    /// Ag specific cells proliferate and som of them start to express the receptor
    num_Agsp_no_receptor_d+=time_step*num_Agsp_no_receptor_d*
                            (proliferation_ratio*LT_max_no_receptor_prol_rate_d-LT_no_to_free_rate_per_APC_d*a.num_Ag());


    num_Agsp_free_receptor_d+=time_step*num_Agsp_no_receptor_d*LT_no_to_free_rate_per_APC_d*a.num_Ag()+
                              time_step*num_Agsp_free_receptor_d*
                              (proliferation_ratio*LT_max_free_prol_rate_d-LT_free_to_bound_rate_per_APC_d*a.num_Ag());

    /// (Monocytes, NK or LT can interact only with one cell)
    num_Agsp_bound_receptor_d+=time_step*num_Agsp_free_receptor_d*LT_free_to_bound_rate_per_APC_d*a.num_Ag()*NK.NK_num_Ag()+
                               time_step*num_Agsp_bound_receptor_d*proliferation_ratio*LT_max_bound_prol_rate_d;

}



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