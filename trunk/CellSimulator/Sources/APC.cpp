#include "Includes/Media.h"
#include "Includes/APC.h"
#include "Includes/LT.h"
#include "Includes/NK.h"




APC_cells::APC_cells(double APC_init,
          double APC_max_proliferation_rate_,
          double APC_no_to_free_rate_per_Ag_ ,
          double free_to_bound_rate_per_LT,
          double APC_Ab_binding_rate,
          double APC_exh_rate,
          double IFN_free_prod_rate_,
          double IFN_Ag_prod_rate_,
          double IFN_bound_prod_rate_,
          double IFN_blocked_prod_rate_,
          double TNF_free_prod_rate_,
          double TNF_Ag_prod_rate_,
          double TNF_bound_prod_rate_,
          double TNF_blocked_prod_rate_
          ):

          num_free_d(APC_init),
          num_Ag_d(0),
          num_LT_bound_d(0),
          num_blocked_d (0),
          num_exhausted_d(0),
          IFN_free_prod_rate_d(IFN_free_prod_rate_),
          IFN_Ag_prod_rate_d(IFN_Ag_prod_rate_),
          IFN_bound_prod_rate_d(IFN_bound_prod_rate_),
          IFN_blocked_prod_rate_d(IFN_blocked_prod_rate_),
          TNF_free_prod_rate_d(TNF_free_prod_rate_),
          TNF_Ag_prod_rate_d(TNF_Ag_prod_rate_),
          TNF_bound_prod_rate_d(TNF_bound_prod_rate_),
          TNF_blocked_prod_rate_d (TNF_blocked_prod_rate_),
          APC_max_proliferation_rate_d(APC_max_proliferation_rate_),
          APC_no_to_free_rate_per_Ag_d(APC_no_to_free_rate_per_Ag_),
          APC_free_to_bound_rate_per_LT_d (free_to_bound_rate_per_LT),
          APC_Ab_binding_rate_d (APC_Ab_binding_rate),
          APC_exh_rate_d (APC_exh_rate)
          {}


APC_cells::APC_cells(const SimParameters& sp,
          const Treatment& tr):

          num_free_d(sp.init_ratio_APC_cells_*tr.init_cells),
          num_Ag_d(0),
          num_LT_bound_d(0),
          num_blocked_d (0),
          num_exhausted_d(0),
          IFN_free_prod_rate_d(sp.APC_IFN_free_prod_rate_),
          IFN_Ag_prod_rate_d(sp.APC_IFN_Ag_prod_rate_),
          IFN_bound_prod_rate_d(sp.APC_IFN_bound_prod_rate_),
          IFN_blocked_prod_rate_d(sp.APC_IFN_blocked_prod_rate_),
          TNF_free_prod_rate_d(sp.APC_TNF_free_prod_rate_),
          TNF_Ag_prod_rate_d(sp.APC_TNF_Ag_prod_rate_),
          TNF_bound_prod_rate_d(sp.APC_TNF_bound_prod_rate_),
          TNF_blocked_prod_rate_d (sp.APC_TNF_blocked_prod_rate_),
          APC_max_proliferation_rate_d(sp.APC_max_proliferation_rate_),
          APC_no_to_free_rate_per_Ag_d(sp.APC_no_to_free_rate_per_Ag_),
          APC_free_to_bound_rate_per_LT_d (sp.APC_free_to_bound_rate_per_LT_),
          APC_Ab_binding_rate_d (sp.APC_Ab_binding_rate_),
          APC_exh_rate_d (sp.APC_exh_rate)
          {}

APC_cells::APC_cells(){}




APC_cells::APC_cells(const APC_cells& other):
    num_free_d(other.num_free_d),
    num_Ag_d(other.num_Ag_d),
    num_LT_bound_d(other.num_LT_bound_d),
    num_blocked_d (other.num_blocked_d),
    num_exhausted_d(other.num_exhausted_d),
    IFN_free_prod_rate_d(other.IFN_free_prod_rate_d),
    IFN_Ag_prod_rate_d(other.IFN_Ag_prod_rate_d),
    IFN_bound_prod_rate_d(other.IFN_bound_prod_rate_d),
    IFN_blocked_prod_rate_d(other.IFN_blocked_prod_rate_d),
    TNF_free_prod_rate_d(other.TNF_free_prod_rate_d),
    TNF_Ag_prod_rate_d(other.TNF_Ag_prod_rate_d),
    TNF_bound_prod_rate_d(other.TNF_bound_prod_rate_d),
    TNF_blocked_prod_rate_d (other.TNF_blocked_prod_rate_d),
    APC_max_proliferation_rate_d(other.APC_max_proliferation_rate_d),
    APC_no_to_free_rate_per_Ag_d(other.APC_no_to_free_rate_per_Ag_d),
    APC_free_to_bound_rate_per_LT_d (other.APC_free_to_bound_rate_per_LT_d),
    APC_Ab_binding_rate_d (other.APC_Ab_binding_rate_d),
    APC_exh_rate_d (other.APC_exh_rate_d)
    {}

APC_cells& operator=(const APC_cells& other)
{
    if (this!=&other)
    {
        APC_cells tmp(other);
        swap(*this,tmp);
    }
    return *this;
}

friend void swap(APC_cells& one, APC_cells& other)
{
    std::swap(one.num_free_d,other.num_free_d);
    std::swap(one.num_Ag_d,other.num_Ag_d);
    std::swap(one.num_LT_bound_d,other.num_LT_bound_d);
    std::swap(one.num_blocked_d ,other.num_blocked_d);
    std::swap(one.num_exhausted_d,other.num_exhausted_d);
    std::swap(one.IFN_free_prod_rate_d,other.IFN_free_prod_rate_d);
    std::swap(one.IFN_Ag_prod_rate_d,other.IFN_Ag_prod_rate_d);
    std::swap(one.IFN_bound_prod_rate_d,other.IFN_bound_prod_rate_d);
    std::swap(one.IFN_blocked_prod_rate_d,other.IFN_blocked_prod_rate_d);
    std::swap(one.TNF_free_prod_rate_d,other.TNF_free_prod_rate_d);
    std::swap(one.TNF_Ag_prod_rate_d,other.TNF_Ag_prod_rate_d);
    std::swap(one.TNF_bound_prod_rate_d,other.TNF_bound_prod_rate_d);
    std::swap(one.TNF_blocked_prod_rate_d ,other.TNF_blocked_prod_rate_d);
    std::swap(one.APC_max_proliferation_rate_d,other.APC_max_proliferation_rate_d);
    std::swap(one.APC_no_to_free_rate_per_Ag_d,other.APC_no_to_free_rate_per_Ag_d);
    std::swap(one.APC_free_to_bound_rate_per_LT_d ,other.APC_free_to_bound_rate_per_LT_d);
    std::swap(one.APC_Ab_binding_rate_d ,other.APC_Ab_binding_rate_d);
    std::swap(one.APC_exh_rate_d ,other.APC_exh_rate_d);

}




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
               (proliferation_ratio*APC_max_proliferation_rate_d-
               APC_free_to_bound_rate_per_LT_d*LT.num_cells_expressing_receptor_and_free()-
               APC_free_to_bound_rate_per_LT_d*NK.NK_num_Ag()-
               APC_free_to_bound_rate_per_LT_d*num_Ag_d-
               APC_Ab_binding_rate_d*m.Ab());

    /// the cells that have interacted with LT grow accordingly with the number of cells that have internalized the Ag and the
    /// number of APC cells expressing the ligand and receptor and bound with monocytes, NK or LT cells with the same rate (Monocytes, NK or LT can
    /// interact only with one cell). We are supposing that all activated cells express receptor and ligand. We are supposing that probabiliities of
    /// interactionts between cells are similar.
    num_LT_bound_d+=num_LT_bound_d*time_step*(proliferation_ratio*APC_max_proliferation_rate_d -
                                              APC_exh_rate_d)+
                    num_Ag_d*time_step*
                   (APC_free_to_bound_rate_per_LT_d*LT.num_cells_expressing_receptor_and_free()+
                    APC_free_to_bound_rate_per_LT_d*NK.NK_num_Ag()+
                    APC_free_to_bound_rate_per_LT_d*num_Ag_d);


    /// the cells that binds the Ab
    num_blocked_d+=num_Ag_d*time_step*APC_Ab_binding_rate_d*m.Ab()+
                   num_blocked_d*time_step*
                   (proliferation_ratio*APC_max_proliferation_rate_d -
                    APC_exh_rate_d);



    /// the cells that have interacted with LT get exhausted
    num_exhausted_d+=num_LT_bound_d*APC_exh_rate_d*time_step+
                     num_blocked_d*time_step*APC_exh_rate_d+
                     num_exhausted_d*time_step*proliferation_ratio*APC_max_proliferation_rate_d;


}
void APC_cells::reset(const SimParameters& sp,
                      const Treatment& tr)
{
    num_free_d=sp.init_ratio_APC_cells_*tr.init_cells;
    num_Ag_d=0;
    num_LT_bound_d=0;
    num_blocked_d=0;
    num_exhausted_d=0;
 }

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

double& APC_cells::num_free()
    {
        return num_free_d;
    }

const double& APC_cells::num_free()const
    {
        return num_free_d;
    }


double APC_cells::num_bound() const
    {
        return num_LT_bound_d;
    }

double APC_cells::percentage_cell_expressing_receptor()const
{
return (num_Ag_d+num_LT_bound_d+num_blocked_d)/num()*100;
}


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
