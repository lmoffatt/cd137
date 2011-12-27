#include "Includes/Media.h"
#include "Includes/APC.h"
#include "Includes/LT.h"
#include "Includes/NK.h"
#include <cmath>

Media::Media(double IFNgamma_init,
             double TNF_init,
             double Tymidine_incorprated_init,
             double TNF_deg_init,
             double IFN_deg_init,
             double init_num_cells,
             double init_Ag,
             double init_Ab,
             double Prol_TymTr_init
                ):

    IFNgamma_d(IFNgamma_init),
    TNF_d(TNF_init),
    Tymidine_incorporated_d (Tymidine_incorprated_init),
    TNF_deg_d(TNF_deg_init),
    IFN_deg_d(IFN_deg_init),
    num_cells_d(init_num_cells),
    Ag_d(init_Ag),
    Ab_d(init_Ab),
    TymidineTriteate_d(0),
    Prol_TymTr_d(Prol_TymTr_init)

   {}

/*

Media::Media(const SimParameters& sp,
             const Treatment& tr):

    IFNgamma_d(0),
    TNF_d(0),
    Tymidine_incorporated_d(0),
    TNF_deg_d(sp.TNF_deg_),
    IFN_deg_d(sp.IFN_deg_),
    num_cells_d(tr.init_cells),
    Ag_d(tr.Ag),
    Ab_d(tr.Ab),
    TymidineTriteate_d(sp.TymidineTriteate_),
    Prol_TymTr_d(sp.Prol_TymTr_)

{}
*/
Media::Media(const Parameters& p,
             const Treatment& tr):

    IFNgamma_d(0),
    TNF_d(0),
    Tymidine_incorporated_d(0),
    TNF_deg_d(p.mean("TNF_deg")),
    IFN_deg_d(p.mean("IFN_deg")),
    num_cells_d(tr.init_cells),
    Ag_d(tr.Ag),
    Ab_d(tr.Ab),
    TymidineTriteate_d(0),
    Prol_TymTr_d(p.mean("Prol_TymTr")){}

Media::Media(const Media& other):
    IFNgamma_d(other.IFNgamma_d),
    TNF_d(other.TNF_d),
    Tymidine_incorporated_d(other.Tymidine_incorporated_d),
    TNF_deg_d(other.TNF_deg_d),
    IFN_deg_d(other.IFN_deg_d),
    num_cells_d(other.num_cells_d),
    Ag_d(other.Ag_d),
    Ab_d(other.Ab_d),
    TymidineTriteate_d(other.TymidineTriteate_d),
    Prol_TymTr_d(other.Prol_TymTr_d)

  {}


Media&
Media::operator=(const Media& other)
{
    if (this!=&other)
    {
        Media tmp(other);
        swap(*this,tmp);
    }
    return *this;
}


void swap(Media& one, Media& other)
{   std::swap(one.IFNgamma_d,other.IFNgamma_d);
    std::swap(one.TNF_d,other.TNF_d);
    std::swap(one.Tymidine_incorporated_d,other.Tymidine_incorporated_d),
    std::swap(one.TNF_deg_d,other.TNF_deg_d);
    std::swap(one.IFN_deg_d, other.IFN_deg_d);
    std::swap(one.num_cells_d,other.num_cells_d);
    std::swap(one.Ag_d,other.Ag_d);
    std::swap(one.Ab_d,other.Ab_d);
    std::swap(one.TymidineTriteate_d,other.TymidineTriteate_d);
    std::swap(one.Prol_TymTr_d,other.Prol_TymTr_d);
 }



double& Media::IFNgamma()
    {
        return IFNgamma_d;
    }

const double& Media::IFNgamma()const
    {
        return IFNgamma_d;
    }

double& Media::TNF()
    {
        return TNF_d;
    }

const double& Media::TNF()const
    {
        return TNF_d;
    }

double& Media::TNF_degradation ()
   {
   return TNF_deg_d;
   }

const double& Media::TNF_degradation () const
   {
   return TNF_deg_d;
   }

double& Media::IFN_degradation ()
   {
   return IFN_deg_d;
   }

const double& Media::IFN_degradation () const
   {
   return IFN_deg_d;
   }



double& Media::num_cells()
    {
        return num_cells_d;
    }

const double& Media::num_cells() const
    {
        return num_cells_d;
    }

double& Media::Ag()
    {
        return Ag_d;
    }

const double& Media::Ag()const
    {
        return Ag_d;
    }

double& Media::Ab()
    {
        return Ab_d;
    }

const double& Media::Ab()const
    {
        return Ab_d;
    }

double& Media::TymidineTriteate()
    {
    return TymidineTriteate_d;
    }
const double& Media::TymidineTriteate() const
    {
    return TymidineTriteate_d;
    }

double& Media::Prol_TymTr()
   {
    return Prol_TymTr_d;
   }

const double& Media::Prol_TymTr() const
    {
     return Prol_TymTr_d;
    }

/// Tymidine incorporated by NK cells
double& Media::Tymidine_incorporated()
{ return Tymidine_incorporated_d;}

const double& Media::Tymidine_incorporated()const
 { return Tymidine_incorporated_d;}

/// Media Main step of the simulation
/// it updates the state of the media and the cells populations for a given time period
void Media::update(double& time_step,double t_run,const APC_cells& APC ,const NK_cells& NK,const LT_cells& LT)
{

    /// IFN is increased by the production rate of each population;

    double IFNgamma_delta=(APC.APC_IFNgamma_production_rate()+
                NK.NK_IFNgamma_production_rate()+
                LT.LT_IFNgamma_production_rate()-
                IFNgamma_d*IFN_degradation())*time_step;
    IFNgamma_d+=IFNgamma_delta;
    /// TNF is increased by the production rate of each population;

    double TNF_delta=(APC.APC_TNF_production_rate()+
           NK.NK_TNF_production_rate()+
           LT.TNF_production_rate()-
           TNF_d*TNF_degradation())*time_step;
    TNF_d+=TNF_delta;
    /// The total number of cells is the adittion of APC + NK + LT
    num_cells_d=APC.num_APC()+NK.num_NK()+LT.num_LT();
    /// Tymidine Pulse at 114

    double Tymidine_incorporated_delta;

    if (t_run<104)
     {TymidineTriteate_d=0;}
      else
    {
        TymidineTriteate_d=1;
    }

     Tymidine_incorporated_delta=APC.APC_TymTr_incorporated()
             +NK.NK_TymTr_incorporated()
             +LT.LT_TymTr_incorporated();
    Tymidine_incorporated_d+=Tymidine_incorporated_delta;
}


 std::ostream& operator<<(std::ostream& s, const Media& c)
{
    s<<"\n IFNgamma_d \t"<<c.IFNgamma_d;
    s<<"\n TNF_d \t"<<c.TNF_d;
    s<<"\n Tymidine Incorprated\t"<<c.Tymidine_incorporated_d;
    s<<"\n num_cells_d \t"<<c.num_cells_d;
    if(0)
    s<<"\n Ag_d \t"<<c.Ag_d;
    s<<"\n Ab_d \t"<<c.Ab_d;
    s<<"\n TNF deg \t"<<c.TNF_deg_d;
    s<<"\n IFN deg \t"<<c.IFN_deg_d;
    s<<"\n TymidineTriteate_d \t"<<c.TymidineTriteate_d;
    s<<"\n Prol_TymTr_d \t"<<c.Prol_TymTr_d;


  return s;
}



