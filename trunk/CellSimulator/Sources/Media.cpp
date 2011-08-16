#include "Includes/Media.h"
#include "Includes/APC.h"
#include "Includes/LT.h"
#include "Includes/NK.h"

Media::Media( double max_num_cells,
              double init_num_cells,
              double Ag_=0,
              double Ab_=0,
              double IFNgamma_init=0,
              double TNF_init=0,
              double TNF_deg_init=0,
              double IFN_deg_init=0):
    //   double internalization_init=0):

    IFNgamma_d(IFNgamma_init),
    TNF_d(TNF_init),
    max_num_cells_d(max_num_cells),
    num_cells_d(init_num_cells),
    Ag_d(Ag_),
    Ab_d(Ab_),
    TNF_deg(TNF_deg_init),
    IFN_deg(IFN_deg_init)
  /*  Ag_internalization_rate (internalization_init)*/ {}

Media::Media(const SimParameters& sp,
             const Treatment& tr):
    //   double internalization_init=0):

    IFNgamma_d(0),
    TNF_d(0),
    max_num_cells_d(sp.max_num_cells_),
    num_cells_d(tr.init_cells),
    Ag_d(tr.Ag),
    Ab_d(tr.Ab),
    TNF_deg(sp.TNF_deg),
    IFN_deg(sp.IFN_deg)
{}


Media::Media(const Media& other):
    IFNgamma_d(other.IFNgamma_d),
    TNF_d(other.TNF_d),
    max_num_cells_d(other.max_num_cells_d),
    num_cells_d(other.num_cells_d),
    Ag_d(other.Ag_d),
    Ab_d(other.Ab_d),
    TNF_deg(other.TNF_deg),
    IFN_deg(other.IFN_deg)
  /*  Ag_internalization_rate (internalization_init)*/ {}


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
    std::swap(one.max_num_cells_d,other.max_num_cells_d);
    std::swap(one.num_cells_d,other.num_cells_d);
    std::swap(one.Ag_d,other.Ag_d);
    std::swap(one.Ab_d,other.Ab_d);
    std::swap(one.TNF_deg,other.TNF_deg);
    std::swap(one.IFN_deg,other.IFN_deg);
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
   return TNF_deg;
   }

const double& Media::IFN_degradation () const
   {
   return IFN_deg;
   }

double& Media::IFN_degradation ()
   {
   return IFN_deg;
   }

const double& Media::TNF_degradation () const
   {
   return TNF_deg;
   }

double& Media::Max_num_cells()
    {
        return max_num_cells_d;
    }

const double& Media::Max_num_cells() const
    {
        return max_num_cells_d;
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
/// Media Main step of the simulation
/// it updates the state of the media and the cells populations for a given time period
void Media::update(double time_step,const APC_cells& APC, const NK_cells& NK, const LT_cells& LT)
{
    /// IFN is increased by the production rate of each population;
    IFNgamma_d+=APC.IFNgamma_production_rate()*time_step+
		NK.IFNgamma_production_rate()*time_step +
		LT.IFNgamma_production_rate()*time_step -
                IFNgamma_d*time_step*IFN_degradation();
    /// TNF is increased by the production rate of each population;
    TNF_d+=APC.TNF_production_rate()*time_step+
	   NK.TNF_production_rate()*time_step +
	   LT.TNF_production_rate()*time_step -
           TNF_d*time_step*TNF_degradation();
    /// The total number of cells is the adittion of APC + NK + LT
    num_cells_d=APC.num()+NK.num()+LT.num();

}


 std::ostream& operator<<(std::ostream& s, const Media& c)
{
    s<<"\n IFNgamma_d \t"<<c.IFNgamma_d;
    s<<"\n TNF_d \t"<<c.TNF_d;
    s<<"\n max_num_cells_d \t"<<c.max_num_cells_d;
    s<<"\n num_cells_d \t"<<c.num_cells_d;
    s<<"\n Ag_d \t"<<c.Ag_d;
    s<<"\n Ab_d \t"<<c.Ab_d;
    s<<"\n TNF_deg \t"<<c.TNF_deg;
    s<<"\n IFN_deg \t"<<c.IFN_deg;
  return s;
}



