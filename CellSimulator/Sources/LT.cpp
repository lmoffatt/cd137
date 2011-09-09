#include "Includes/Media.h"
#include "Includes/APC.h"
#include "Includes/LT.h"
#include "Includes/NK.h"


LT_cells::LT_cells(double num_LT_init_,
                   double LT_num_specific_,
                   double LT_max_no_receptor_prol_rate_,
                   double LT_max_free_prol_rate_,
                   double LT_max_bound_prol_rate_,
                   double LT_max_blocked_prol_rate_,
                   double IFN_no_rec_prod_rate_,
                   double IFN_free_prod_rate_,
                   double IFN_bound_prod_rate_,
                   double IFN_blocked_prod_rate_,
                   double TNF_no_rec_prod_rate_,
                   double TNF_free_prod_rate_,
                   double TNF_bound_prod_rate_,
                   double TNF_blocked_prod_rate_,
                   double LT_no_to_free_rate_per_APC_,
                   double LT_free_to_bound_rate_per_APC_,
                   double LT_mAb_binding_rate_,
                   double LT_exh_rate_):

    num_non_Agsp_d(num_LT_init_),
    num_Agsp_no_receptor_d(LT_num_specific_),
    num_Agsp_free_receptor_d(0),
    num_Agsp_bound_receptor_d(0),
    num_blocked_d(0),
    IFN_no_rec_prod_rate_d(IFN_no_rec_prod_rate_),
    IFN_free_prod_rate_d(IFN_free_prod_rate_),
    IFN_bound_prod_rate_d(IFN_bound_prod_rate_),
    IFN_blocked_prod_rate_d(IFN_blocked_prod_rate_),
    TNF_no_rec_prod_rate_d(TNF_no_rec_prod_rate_),
    TNF_free_prod_rate_d(TNF_free_prod_rate_),
    TNF_bound_prod_rate_d(TNF_bound_prod_rate_),
    TNF_blocked_prod_rate_d (TNF_blocked_prod_rate_),
    LT_max_no_receptor_prol_rate_d(LT_max_no_receptor_prol_rate_),
    LT_max_free_prol_rate_d(LT_max_free_prol_rate_),
    LT_max_bound_prol_rate_d(LT_max_bound_prol_rate_),
    LT_max_blocked_prol_rate_d(LT_max_blocked_prol_rate_),
    LT_no_to_free_rate_per_APC_d(LT_no_to_free_rate_per_APC_),
    LT_free_to_bound_rate_per_APC_d (LT_free_to_bound_rate_per_APC_),
    LT_mAb_binding_rate_d (LT_mAb_binding_rate_)
    {}




LT_cells::LT_cells(const LT_cells& other):
    num_non_Agsp_d(other.num_non_Agsp_d),
    num_Agsp_no_receptor_d(other.num_Agsp_no_receptor_d),
    num_Agsp_free_receptor_d(other.num_Agsp_free_receptor_d),
    num_Agsp_bound_receptor_d(other.num_Agsp_bound_receptor_d),
    num_blocked_d(other.num_blocked_d),
    IFN_no_rec_prod_rate_d(other.IFN_no_rec_prod_rate_d),
    IFN_free_prod_rate_d(other.IFN_free_prod_rate_d),
    IFN_bound_prod_rate_d(other.IFN_bound_prod_rate_d),
    IFN_blocked_prod_rate_d(other.IFN_blocked_prod_rate_d),
    TNF_no_rec_prod_rate_d(other.TNF_no_rec_prod_rate_d),
    TNF_free_prod_rate_d(other.TNF_free_prod_rate_d),
    TNF_bound_prod_rate_d(other.TNF_bound_prod_rate_d),
    TNF_blocked_prod_rate_d (other.TNF_blocked_prod_rate_d),
    LT_max_no_receptor_prol_rate_d(other.LT_max_no_receptor_prol_rate_d),
    LT_max_free_prol_rate_d(other.LT_max_free_prol_rate_d),
    LT_max_bound_prol_rate_d(other.LT_max_bound_prol_rate_d),
    LT_max_blocked_prol_rate_d(other.LT_max_blocked_prol_rate_d),
    LT_no_to_free_rate_per_APC_d(other.LT_no_to_free_rate_per_APC_d),
    LT_free_to_bound_rate_per_APC_d (other.LT_free_to_bound_rate_per_APC_d),
    LT_mAb_binding_rate_d (other.LT_mAb_binding_rate_d)
    {}


LT_cells&
LT_cells::operator=(const LT_cells& other)
{
    if (this!=&other)
    {
        LT_cells tmp(other);
        swap(*this,tmp);
    }
    return *this;
}

void swap(LT_cells& one, LT_cells& other)
{
    std::swap(one.num_non_Agsp_d,other.num_non_Agsp_d);
    std::swap(one.num_Agsp_no_receptor_d,other.num_Agsp_no_receptor_d);
    std::swap(one.num_Agsp_free_receptor_d,other.num_Agsp_free_receptor_d);
    std::swap(one.num_Agsp_bound_receptor_d,other.num_Agsp_bound_receptor_d);
    std::swap(one.num_blocked_d,other.num_blocked_d);
    std::swap(one.IFN_no_rec_prod_rate_d,other.IFN_no_rec_prod_rate_d);
    std::swap(one.IFN_free_prod_rate_d,other.IFN_free_prod_rate_d);
    std::swap(one.IFN_bound_prod_rate_d,other.IFN_bound_prod_rate_d);
    std::swap(one.IFN_blocked_prod_rate_d,other.IFN_blocked_prod_rate_d);
    std::swap(one.TNF_no_rec_prod_rate_d,other.TNF_no_rec_prod_rate_d);
    std::swap(one.TNF_free_prod_rate_d,other.TNF_free_prod_rate_d);
    std::swap(one.TNF_bound_prod_rate_d,other.TNF_bound_prod_rate_d);
    std::swap(one.TNF_blocked_prod_rate_d,other.TNF_blocked_prod_rate_d);
    std::swap(one.LT_max_no_receptor_prol_rate_d,other.LT_max_no_receptor_prol_rate_d);
    std::swap(one.LT_max_free_prol_rate_d,other.LT_max_free_prol_rate_d);
    std::swap(one.LT_max_bound_prol_rate_d,other.LT_max_bound_prol_rate_d);
    std::swap(one.LT_max_blocked_prol_rate_d,other.LT_max_blocked_prol_rate_d);
    std::swap(one.LT_no_to_free_rate_per_APC_d,other.LT_no_to_free_rate_per_APC_d);
    std::swap(one.LT_free_to_bound_rate_per_APC_d,other.LT_free_to_bound_rate_per_APC_d);
    std::swap(one.LT_mAb_binding_rate_d,other.LT_mAb_binding_rate_d);


}



/// Main step for LT
void LT_cells::update(double time_step,const Media& m,const APC_cells& a,const NK_cells& NK)

{
    double proliferation_ratio=(m.Max_num_cells()-m.num_cells())/m.Max_num_cells();


    /// cells not sensitive to the Ag proliferate passively
    num_non_Agsp_d+=time_step*num_non_Agsp_d*proliferation_ratio*LT_max_no_receptor_prol_rate_d;


    /// Ag specific cells proliferate and some of them interact with APC and get activated and express the receptor
    num_Agsp_no_receptor_d+=time_step*num_Agsp_no_receptor_d*
                            (proliferation_ratio*LT_max_no_receptor_prol_rate_d-LT_no_to_free_rate_per_APC_d*(a.num_Ag()+a.num_bound()));

    /// Acà hay algo que clarificar (La cèlula T no interactúa solo una vez con la APC???)
    num_Agsp_free_receptor_d+=time_step*num_Agsp_no_receptor_d*LT_no_to_free_rate_per_APC_d*(a.num_Ag()+ a.num_bound())+
                              time_step*num_Agsp_free_receptor_d*
                              (proliferation_ratio*LT_max_free_prol_rate_d-LT_free_to_bound_rate_per_APC_d*
                               -LT_mAb_binding_rate_d*m.Ab());

    /// (Monocytes, NK or LT can interact only with one cell) (There are not LT exhausted at the times of experiment)
    num_Agsp_bound_receptor_d+=time_step*num_Agsp_free_receptor_d*LT_free_to_bound_rate_per_APC_d*+
                               time_step*num_Agsp_bound_receptor_d*proliferation_ratio*LT_max_bound_prol_rate_d-
                               num_Agsp_bound_receptor_d*exh_rate_d*time_step;


    /// LT interact with blocking mAb and grow as LT free rates (There are not LT exhausted at the times of experiment)
    num_blocked_d+= num_Agsp_free_receptor_d*time_step*LT_mAb_binding_rate_d*m.Ab() +
                    time_step*num_blocked_d*proliferation_ratio*LT_max_blocked_prol_rate_d-
                    num_blocked_d*exh_rate_d*time_step;
    /// LT exhausted
    num_exhausted_d+= num_Agsp_bound_receptor_d*exh_rate_d*time_step+
                      num_blocked_d*exh_rate_d*time_step;

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

double LT_cells::num_exhausted () const
{
    return num_exhausted_d;
}
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
    }





std::ostream& operator<<(std::ostream& s, const LT_cells& c)
{

   s<<"\n num_non_Agsp_d \t"<<c.num_non_Agsp_d;

   s<<"\n num_Agsp_no_receptor_d \t"<<c.num_Agsp_no_receptor_d;
   s<<"\n num_Agsp_free_receptor_d \t"<<c.num_Agsp_free_receptor_d;
   s<<"\n num_Agsp_bound_receptor_d \t"<<c.num_Agsp_bound_receptor_d;
   s<<"\n num_blocked_d \t"<<c.num_blocked_d;
   if (0)
   {
   s<<"///\n----------------------------------\n";
   s<<"those are parameters that do not vary\n\n";

   s<<"\n IFN_no_rec_prod_rate_d \t"<<c.IFN_no_rec_prod_rate_d;
   s<<"\n IFN_free_prod_rate_d \t"<<c.IFN_free_prod_rate_d;
   s<<"\n IFN_bound_prod_rate_d \t"<<c.IFN_bound_prod_rate_d;
   s<<"\n IFN_blocked_prod_rate_d \t"<<c.IFN_blocked_prod_rate_d;
   s<<"\n TNF_no_rec_prod_rate_d \t"<<c.TNF_no_rec_prod_rate_d;
   s<<"\n TNF_free_prod_rate_d \t"<<c.TNF_free_prod_rate_d;
   s<<"\n TNF_bound_prod_rate_d \t"<<c.TNF_bound_prod_rate_d;
   s<<"\n TNF_blocked_prod_rate_d \t"<<c.TNF_blocked_prod_rate_d;
   s<<"\n LT_max_no_receptor_prol_rate_d \t"<<c.LT_max_no_receptor_prol_rate_d;
   s<<"\n LT_max_free_prol_rate_d \t"<<c.LT_max_free_prol_rate_d;
   s<<"\n LT_max_bound_prol_rate_d \t"<<c.LT_max_bound_prol_rate_d;
   s<<"\n LT_max_blocked_prol_rate_d \t"<<c.LT_max_blocked_prol_rate_d;
   s<<"\n LT_no_to_free_rate_per_APC_d \t"<<c.LT_no_to_free_rate_per_APC_d;
   s<<"\n LT_free_to_bound_rate_per_APC_d \t"<<c.LT_free_to_bound_rate_per_APC_d;
   s<<"\n LT_mAb_binding_rate_d \t"<<c.LT_mAb_binding_rate_d;
}
   return s;
}


