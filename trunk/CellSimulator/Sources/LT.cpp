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


    /// Ag specific cells proliferate and some of them interact with APC and get activated and express the receptor
    num_Agsp_no_receptor_d+=time_step*num_Agsp_no_receptor_d*
                            (proliferation_ratio*LT_max_no_receptor_prol_rate_d-LT_no_to_free_rate_per_APC_d*a.num_Ag());

    /// Acà hay algo que clarificar (La cèlula T no interactúa solo una vez con la APC???)
    num_Agsp_free_receptor_d+=time_step*num_Agsp_no_receptor_d*LT_no_to_free_rate_per_APC_d*a.num_Ag()+
                              time_step*num_Agsp_free_receptor_d*
                              (proliferation_ratio*LT_max_free_prol_rate_d-LT_free_to_bound_rate_per_APC_d*a.num_Ag()
                               -LT_mAb_binding_rate_d*m.Ab());

    /// (Monocytes, NK or LT can interact only with one cell) (There are not LT exhausted at the times of experiment)
    num_Agsp_bound_receptor_d+=time_step*num_Agsp_free_receptor_d*LT_free_to_bound_rate_per_APC_d*a.num_Ag()+
                               time_step*num_Agsp_free_receptor_d*LT_free_to_bound_rate_per_APC_d*NK.NK_num_Ag()+
                               time_step*num_Agsp_bound_receptor_d*proliferation_ratio*LT_max_bound_prol_rate_d;


    /// LT interact with blocking mAb and grow as LT free rates (There are not LT exhausted at the times of experiment)
    num_blocked_d+= num_Agsp_free_receptor_d*time_step*LT_mAb_binding_rate_d*m.Ab() +
                    time_step*num_blocked_d*proliferation_ratio*LT_max_blocked_prol_rate_d;
};

void LT_cells::reset(const SimParameters& sp,
                      const Treatment& tr)
    {
        num_non_Agsp_d=sp.init_ratio_LT_cells_*tr.init_cells;
        num_Agsp_no_receptor_d=sp.LT_ratio_specific_;
        num_Agsp_free_receptor_d=0;
        num_Agsp_bound_receptor_d=0;
        num_blocked_d=0;

    }


double LT_cells::num() const
    {
        return num_Agsp_bound_receptor_d+num_Agsp_free_receptor_d+num_Agsp_no_receptor_d+num_non_Agsp_d+num_blocked_d;
    }

double LT_cells::IFNgamma_production_rate() const
    {
        return (num_non_Agsp_d+num_Agsp_no_receptor_d)*IFN_no_rec_prod_rate_d+
                num_Agsp_free_receptor_d*IFN_free_prod_rate_d+
                num_Agsp_bound_receptor_d*IFN_bound_prod_rate_d+
                num_blocked_d*IFN_blocked_prod_rate_d;
    };

double LT_cells::TNF_production_rate() const
    {
        return (num_non_Agsp_d+num_Agsp_no_receptor_d)*TNF_no_rec_prod_rate_d+
                num_Agsp_free_receptor_d*TNF_free_prod_rate_d+
                num_Agsp_bound_receptor_d*TNF_bound_prod_rate_d+
                num_blocked_d*TNF_blocked_prod_rate_d;
    };

double LT_cells::num_cells_not_Ag_specific()const
    {
        return num_non_Agsp_d;
    };

double LT_cells::num_cells_not_expressing_receptor()const
    {
        return num_Agsp_no_receptor_d;
    };

double LT_cells::num_blocked()const
    {
        return num_blocked_d;
    };

double LT_cells::LT_percentage_cell_expressing_receptor()const
    {
    return num_cells_expressing_receptor()/num()*100;
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
