#include "Includes/Media.h"
#include "Includes/APC.h"
#include "Includes/LT.h"
#include "Includes/NK.h"


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
void Media::update(double time_step,const APC_cells& APC_, const NK_cells& NK_, const LT_cells& LT_)
{
    /// IFN is increased by the production rate of each population;
    IFNgamma_d+=APC_.IFNgamma_production_rate()*time_step+
                NK_.IFNgamma_production_rate()*time_step +
                LT_.IFNgamma_production_rate()*time_step ;
    /// TNF is increased by the production rate of each population;
    TNF_d+=APC_.TNF_production_rate()*time_step+
           NK_.TNF_production_rate()*time_step +
           LT_.TNF_production_rate()*time_step -
           TNF_d*time_step*TNF_degradation();
    /// The total number of cells is the adittion of APC + NK + LT
    num_cells_d=APC_.num()+NK_.num()+LT_.num();

};

Media::Media( double max_num_cells_,
              double init_num_cells,
              double Ag_=0,
              double Ab_=0,
              double IFNgamma_init=0,
              double TNF_init=0,
              double TNF_deg_init=0):
           //   double internalization_init=0):

              IFNgamma_d(IFNgamma_init),
              TNF_d(TNF_init),
              max_num_cells_d(max_num_cells_),
              num_cells_d(init_num_cells),
              Ag_d(Ag_),
              Ab_d(Ab_),
              TNF_deg(TNF_deg_init)
            /*  Ag_internalization_rate (internalization_init)*/ {}