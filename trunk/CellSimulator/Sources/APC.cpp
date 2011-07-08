#include "Includes/Media.h"
#include "Includes/APC.h"
#include "Includes/LT.h"
#include "Includes/NK.h"
/// Main step for APC
void APC_cells::update(double time_step,const Media& m, const NK_cells& NK, const LT_cells& LT)
{
    /// the proliferation is one for cero cells,
    /// 0 when the number of cells equals the maximum
    /// and negative (apoptosis) when there are more cells than the maximum
    double proliferation_ratio=(m.Max_num_cells()-m.num_cells())/m.Max_num_cells();

    /// we update each subpopulation of cells independently and we take into account the transition from one state to the other

    /// the number of free cells (no Ag) they proliferate according to the cell concentration (factor proliferation ratio)
    /// and some of them are "lost" since they internalize the Ag
    num_free_d+=num_free_d*time_step*
                (proliferation_ratio*APC_max_proliferation_rate_d - APC_no_to_free_rate_per_Ag_d*m.Ag());

    /** the cells that have internalize the Ag proliferate in the same way than the free
    they grow also by the free cells that internalize the Ag
    they shrink by the cells that interact with the LT cells*/
    num_Ag_d+=num_free_d*time_step*APC_no_to_free_rate_per_Ag_d*m.Ag()+
              num_Ag_d*time_step*
              (proliferation_ratio*APC_max_proliferation_rate_d
               -APC_free_to_bound_rate_per_LT_d*LT.num_cells_expressing_receptor_and_free()*NK.NK_num_Ag()-num_Ag_d*time_step*APC_Ab_binding_rate_d*m.Ab());

    /// the cells that have interacted with LT grow accordingly with the number of cells that have internalized the Ag and the
    /// number of APC cells expressing the ligand and receptor and bound with monocytes, NK or LT cells with the same rate (Monocytes, NK or LT can
    /// interact only with one cell). We are supposing that all activated cells express receptor and ligand.
    num_LT_bound_d+=num_Ag_d*time_step*APC_free_to_bound_rate_per_LT_d*LT.num_cells_expressing_receptor_and_free()+
                    num_Ag_d*time_step*APC_free_to_bound_rate_per_LT_d*NK.NK_num_Ag()+
                    num_Ag_d*time_step*APC_free_to_bound_rate_per_LT_d*num_Ag_d+
                    num_LT_bound_d*time_step*(proliferation_ratio*APC_max_proliferation_rate_d -
                            APC_exh_rate_d);

    /// the cells that binds the Ab
    num_blocked_d+=num_Ag_d*time_step*APC_Ab_binding_rate_d*m.Ab()+num_blocked_d*time_step*(proliferation_ratio*APC_max_proliferation_rate_d - APC_exh_rate_d);


    /// the cells that have interacted with LT get exhausted
    num_exhausted_d+=num_LT_bound_d*APC_exh_rate_d*time_step+num_exhausted_d*time_step*proliferation_ratio*APC_max_proliferation_rate_d

                     ;
};

double APC_cells::num() const
    {
        return num_Ag_d+num_free_d+num_LT_bound_d+num_exhausted_d+num_blocked_d;
    }


double APC_cells::IFNgamma_production_rate() const
    {
        double sum=IFN_free_prod_rate_d*num_free_d+
                   IFN_Ag_prod_rate_d*num_Ag_d+
                   IFN_bound_prod_rate_d*num_LT_bound_d+
                   IFN_blocked_prod_rate_d*num_blocked_d;

        return sum;
    }

double APC_cells::TNF_production_rate() const
    {
        double sum=TNF_free_prod_rate_d*num_free_d+
                   TNF_Ag_prod_rate_d*num_Ag_d+
                   TNF_bound_prod_rate_d*num_LT_bound_d+
                   TNF_blocked_prod_rate_d*num_blocked_d;

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

double APC_cells::percentage_cell_expressing_receptor()const
{
return num_Ag_d/num()*100;
};

const double& APC_cells::num_blocked() const
    {
        return num_blocked_d;
    }

double& APC_cells::num_blocked()
    {
      return num_blocked_d;
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
