#include "Includes/Media.h"
#include "Includes/APC.h"
#include "Includes/NK.h"
#include "Includes/LT.h"

/// main step for the NK cells
void NK_cells::update(double time_step,const Media& m,const APC_cells& APC, const LT_cells& LT)
{
    /// the proliferation is one for cero cells,
    /// 0 when the number of cells equals the maximum
    /// and negative (apoptosis) when there are more cells than the maximum
    double proliferation_ratio=(m.Max_num_cells()-m.num_cells())/m.Max_num_cells();

    /// we update each subpopulation of cells independently and we take into account the transition from one state to the other

    /// the number of free cells (no Ag) they proliferate according to the cell concentration (factor proliferation ratio)
    /// and some of them are "lost" since they internalize the Ag
    NK_num_free_d+=NK_num_free_d*time_step*
                (proliferation_ratio*NK_max_proliferation_rate_d-NK_no_to_free_rate_per_Ag_d*m.Ag());

    /** the cells that have internalize the Ag proliferate in the same way than the free
    they grow also by the free cells that internalize the Ag
    they shrink by the cells that interact with the LT cells*/
    NK_num_Ag_d+=NK_num_free_d*time_step*NK_no_to_free_rate_per_Ag_d*m.Ag()+
                 NK_num_Ag_d*time_step*
                 (proliferation_ratio*NK_max_proliferation_rate_d-
                 NK_free_to_bound_rate_per_LT_d*LT.num_cells_expressing_receptor_and_free()-
                 NK_free_to_bound_rate_per_LT_d*APC.num_Ag()-
                 NK_free_to_bound_rate_per_LT_d*NK_num_Ag_d-
                 NK_Ab_binding_rate_d*m.Ab());


    /// the cells that have interacted with LT grow accordingly with the number of cells that have internalized the Ag and the
    /// number of LT cells, monocytes and NK expressing the receptor with the same affinity (Monocytes, NK or LT can
    /// interact only with one cell). We are supposing that all activated cells express receptor and ligand.
    NK_num_LT_bound_d+= NK_num_LT_bound_d*time_step*
                        (proliferation_ratio*NK_max_proliferation_rate_d -
                        NK_exh_rate_d)+
                        NK_num_Ag_d*time_step*
                       (NK_free_to_bound_rate_per_LT_d*LT.num_cells_expressing_receptor_and_free()+
                        NK_free_to_bound_rate_per_LT_d*APC.num_Ag()+
                        NK_free_to_bound_rate_per_LT_d*NK_num_Ag_d);


    /// NK cells can interact with other cells or with the blocking mAb
    NK_blocked_d+=NK_num_Ag_d*time_step*NK_Ab_binding_rate_d*m.Ab()+
                  NK_blocked_d*time_step*
                  (proliferation_ratio*NK_max_proliferation_rate_d-
                   NK_exh_rate_d);

    /// the cells that have interacted with LT get exhausted
    NK_num_exhausted_d+=NK_num_exhausted_d*time_step*proliferation_ratio*NK_max_proliferation_rate_d+
                        NK_num_LT_bound_d*NK_exh_rate_d*time_step+
                        NK_blocked_d*time_step*NK_exh_rate_d;


};


double NK_cells::num() const
    {
        return NK_num_Ag_d+NK_num_free_d+NK_num_LT_bound_d+NK_num_exhausted_d;
    }

double NK_cells::IFNgamma_production_rate() const
    {
        double sum=NK_IFN_free_prod_rate_d*NK_num_free_d+
                   NK_IFN_Ag_prod_rate_d*NK_num_Ag_d+
                   NK_IFN_bound_prod_rate_d*NK_num_LT_bound_d+
                   NK_IFN_blocked_prod_rate_d*NK_blocked_d+
                   NK_num_exhausted_d*0;
        return sum;
    }

double NK_cells::TNF_production_rate() const
    {
        double sum=NK_TNF_free_prod_rate_d*NK_num_free_d+
                   NK_TNF_Ag_prod_rate_d*NK_num_Ag_d+
                   NK_TNF_bound_prod_rate_d*NK_num_LT_bound_d+
                   NK_TNF_blocked_prod_rate_d*NK_blocked_d+
                   NK_num_exhausted_d*0;
        return sum;
    }


double& NK_cells::NK_num_free()
    {
        return NK_num_free_d;
    }

const double& NK_cells::NK_num_free()const
    {
        return NK_num_free_d;
    }


double& NK_cells::NK_num_Ag()
    {
        return NK_num_Ag_d;
    }

const double& NK_cells::NK_num_Ag()const
    {
        return NK_num_Ag_d;
    }

double& NK_cells::NK_blocked()
    {
        return NK_blocked_d;
    }

const double& NK_cells::NK_blocked()const
    {
        return NK_blocked_d;
    }

double NK_cells::percentage_cell_expressing_receptor()const
    {
        return NK_num_Ag_d/num()*100;
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
