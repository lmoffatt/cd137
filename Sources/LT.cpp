#include "Includes/Media.h"
#include "Includes/APC.h"
#include "Includes/LT.h"
#include "Includes/NK.h"


LT_cells::LT_cells(/// 1) Init number of LT
                   /*1*/ double ratio_init_LTns_,
                   /*2*/ double ratio_initLTspecific_,

                   /// 2) IFN Poductions rates of each type of LT
                   /*3*/ double IFN_LTns_prod_rate_,
                   /*4*/ double IFN_LTbo_prod_rate_,
                   /*5*/ double IFN_LTbl_prod_rate_,

                   /// 3) TNF Poductions rates of each type of LT
                   /*6*/ double TNF_LTns_prod_rate_,
                   /*7*/ double TNF_LTbo_prod_rate_,
                   /*8*/ double TNF_LTbl_prod_rate_,

                   /// 4) Percentages of IFN productions of each type of LT
                   /*9*/ double percentage_IFN_LTns_prod_rate_,
                   /*10*/ double percentage_IFN_LTbo_prod_rate_,
                   /*11*/ double percentage_IFN_LTbl_prod_rate_,


                   /// 5)Percentages of TNF productions of each type of LT
                   /*12*/ double percentage_TNF_LTns_prod_rate_,
                   /*13*/ double percentage_TNF_LTbo_prod_rate_,
                   /*14*/ double percentage_TNF_LTbl_prod_rate_,

                   /// 6) Proliferation rates
                   /*15*/ double LTns_proliferation_rate_,
                   /*16*/ double LTbo_proliferation_rate_,
                   /*17*/ double LTbl_proliferation_rate_,

                   /// 7) Apoptosis rates
                   /*18*/ double LTns_apop_rate_,
                   /*19*/ double LTbo_apop_rate_,
                   /*20*/ double LTbl_apop_rate_,
                   /*21*/ double LTexh_apop_rate_,

                   /// 8) constant saturation of TNF for apoptosis
                   /*22*/ double Ks_LT_m_TNF_,

                   /// 9) Percentages of cell expressing receptor
                   /*23*/ double LTns_expressing_receptor_,

                   /// 10) Apoptosis rate for TNF
                   /*24*/ double u_LT_TNF_,

                   /// 11) LT exh rate
                   /*25*/ double LT_exh_rate_,

                   /// 12) apoptosis related parameters
                   /*26*/ double t_apop_meas_,
                   /*27*/ double t_duration_apoptosis_

                   ):



    LTns_d(ratio_init_LTns_),
    LT0_d(ratio_initLTspecific_),
    LTbo_d(0),
    LTbl_d(0),
    LTexh_d(0),
    LT_TymTr_incorporated_d(0),
    Total_cells_in_apoptosis_d(0),

    IFN_LTns_prod_rate_d(IFN_LTns_prod_rate_),
    IFN_LTbo_prod_rate_d(IFN_LTbo_prod_rate_),
    IFN_LTbl_prod_rate_d(IFN_LTbl_prod_rate_),

    TNF_LTns_prod_rate_d(TNF_LTns_prod_rate_),
    TNF_LTbo_prod_rate_d(TNF_LTbo_prod_rate_),
    TNF_LTbl_prod_rate_d(TNF_LTbl_prod_rate_),

    percentage_IFN_LTns_prod_rate_d(percentage_IFN_LTns_prod_rate_),
    percentage_IFN_LTbo_prod_rate_d( percentage_IFN_LTbo_prod_rate_),
    percentage_IFN_LTbl_prod_rate_d(percentage_IFN_LTbl_prod_rate_),

    percentage_TNF_LTns_prod_rate_d(percentage_TNF_LTns_prod_rate_),
    percentage_TNF_LTbo_prod_rate_d(percentage_TNF_LTbo_prod_rate_),
    percentage_TNF_LTbl_prod_rate_d(percentage_TNF_LTbl_prod_rate_),

    LTns_proliferation_rate_d(LTns_proliferation_rate_),
    LTbo_proliferation_rate_d(LTbo_proliferation_rate_),
    LTbl_proliferation_rate_d(LTbl_proliferation_rate_),

    LTns_apop_rate_d(LTns_apop_rate_),
    LTbo_apop_rate_d(LTbo_apop_rate_),
    LTbl_apop_rate_d(LTbl_apop_rate_),
    LTexh_apop_rate_d(LTexh_apop_rate_),

    Ks_LT_m_TNF_d(Ks_LT_m_TNF_),

    LTns_expressing_receptor_d(LTns_expressing_receptor_),

    u_LT_TNF_d(u_LT_TNF_),

    LT_exh_rate_d(LT_exh_rate_),

    t_apop_meas_d (t_apop_meas_),
    t_duration_apoptosis_d(t_duration_apoptosis_)

    {}




LT_cells::LT_cells(const LT_cells& other):
    LTns_d(other.LTns_d),
    LT0_d(other.LT0_d),
    LTbo_d(other.LTbo_d),
    LTbl_d(other.LTbl_d),
    LTexh_d(other.LTexh_d),
    LT_TymTr_incorporated_d(other.LT_TymTr_incorporated_d),
    Total_cells_in_apoptosis_d(other.Total_cells_in_apoptosis_d),

    IFN_LTns_prod_rate_d(other.IFN_LTns_prod_rate_d),
    IFN_LTbo_prod_rate_d(other.IFN_LTbo_prod_rate_d),
    IFN_LTbl_prod_rate_d(other.IFN_LTbl_prod_rate_d),

    TNF_LTns_prod_rate_d(other.TNF_LTns_prod_rate_d),
    TNF_LTbo_prod_rate_d(other.TNF_LTbo_prod_rate_d),
    TNF_LTbl_prod_rate_d(other.TNF_LTbl_prod_rate_d),

    percentage_IFN_LTns_prod_rate_d(other.percentage_IFN_LTns_prod_rate_d),
    percentage_IFN_LTbo_prod_rate_d(other. percentage_IFN_LTbo_prod_rate_d),
    percentage_IFN_LTbl_prod_rate_d(other.percentage_IFN_LTbl_prod_rate_d),

    percentage_TNF_LTns_prod_rate_d(other.percentage_TNF_LTns_prod_rate_d),
    percentage_TNF_LTbo_prod_rate_d(other.percentage_TNF_LTbo_prod_rate_d),
    percentage_TNF_LTbl_prod_rate_d(other.percentage_TNF_LTbl_prod_rate_d),

    LTns_proliferation_rate_d(other.LTns_proliferation_rate_d),
    LTbo_proliferation_rate_d(other.LTbo_proliferation_rate_d),
    LTbl_proliferation_rate_d(other.LTbl_proliferation_rate_d),

    LTns_apop_rate_d(other.LTns_apop_rate_d),
    LTbo_apop_rate_d(other.LTbo_apop_rate_d),
    LTbl_apop_rate_d(other.LTbl_apop_rate_d),
    LTexh_apop_rate_d(other.LTexh_apop_rate_d),

    Ks_LT_m_TNF_d(other.Ks_LT_m_TNF_d),

    LTns_expressing_receptor_d(other.LTns_expressing_receptor_d),

    u_LT_TNF_d(other.u_LT_TNF_d),

    LT_exh_rate_d(other.LT_exh_rate_d),
    t_apop_meas_d (other.t_apop_meas_d),
    t_duration_apoptosis_d(other.t_duration_apoptosis_d)

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
    std::swap(one.LTns_d,other.LTns_d);
    std::swap(one.LT0_d,other.LT0_d);
    std::swap(one.LTbo_d,other.LTbo_d);
    std::swap(one.LTbl_d,other.LTbl_d);
    std::swap(one.LTexh_d,other.LTexh_d);
    std::swap(one.LT_TymTr_incorporated_d,other.LT_TymTr_incorporated_d);
    std::swap(one.Total_cells_in_apoptosis_d,other.Total_cells_in_apoptosis_d);

    std::swap(one.IFN_LTns_prod_rate_d,other.IFN_LTns_prod_rate_d);
    std::swap(one.IFN_LTbo_prod_rate_d,other.IFN_LTbo_prod_rate_d);
    std::swap(one.IFN_LTbl_prod_rate_d,other.IFN_LTbl_prod_rate_d);

    std::swap(one.TNF_LTns_prod_rate_d,other.TNF_LTns_prod_rate_d);
    std::swap(one.TNF_LTbo_prod_rate_d,other.TNF_LTbo_prod_rate_d);
    std::swap(one.TNF_LTbl_prod_rate_d,other.TNF_LTbl_prod_rate_d);

    std::swap(one.percentage_IFN_LTns_prod_rate_d,other.percentage_IFN_LTns_prod_rate_d);
    std::swap(one.percentage_IFN_LTbo_prod_rate_d,other. percentage_IFN_LTbo_prod_rate_d);
    std::swap(one.percentage_IFN_LTbl_prod_rate_d,other.percentage_IFN_LTbl_prod_rate_d);

    std::swap(one.percentage_TNF_LTns_prod_rate_d,other.percentage_TNF_LTns_prod_rate_d);
    std::swap(one.percentage_TNF_LTbo_prod_rate_d,other.percentage_TNF_LTbo_prod_rate_d);
    std::swap(one.percentage_TNF_LTbl_prod_rate_d,other.percentage_TNF_LTbl_prod_rate_d);

    std::swap(one.LTns_proliferation_rate_d,other.LTns_proliferation_rate_d);
    std::swap(one.LTbo_proliferation_rate_d,other.LTbo_proliferation_rate_d);
    std::swap(one.LTbl_proliferation_rate_d,other.LTbl_proliferation_rate_d);

    std::swap(one.LTns_apop_rate_d,other.LTns_apop_rate_d);
    std::swap(one.LTbo_apop_rate_d,other.LTbo_apop_rate_d);
    std::swap(one.LTbl_apop_rate_d,other.LTbl_apop_rate_d);
    std::swap(one.LTexh_apop_rate_d,other.LTexh_apop_rate_d);

    std::swap(one.Ks_LT_m_TNF_d,other.Ks_LT_m_TNF_d);

    std::swap(one.LTns_expressing_receptor_d,other.LTns_expressing_receptor_d);

    std::swap(one.u_LT_TNF_d,other.u_LT_TNF_d);

    std::swap(one.LT_exh_rate_d,other.LT_exh_rate_d);

    std::swap(one.t_apop_meas_d ,other.t_apop_meas_d);
    std::swap(one.t_duration_apoptosis_d,other.t_duration_apoptosis_d);




}



/// Main step for LT
void LT_cells::update(double time_step, double t_run, const Media& m,const APC_cells& APC,const NK_cells& NK)

{    
    /// cells not sensitive to the Ag proliferate passively
    LTns_d+=(-LTns_apop_rate_d*LTns_d+LTns_proliferation_rate_d*LTns_d)*time_step;


    /// Ag specific cells proliferate and some of them interact with APC and get activated and express the receptor
    LT0_d+=(-LTns_apop_rate_d*LT0_d+LTns_proliferation_rate_d*LT0_d-APC.APC_LT_1()*LT0_d*(APC.APCa()/(APC.APCa()+ APC.KsAPC_LT()))-
            APC.APC_LT_2()*LT0_d*(APC.APCbo()/(APC.APCbo()+ APC.KsAPC_LT()))-
            APC.APC_LT_2()*LT0_d*(APC.APCbo_Ab()/(APC.APCbo_Ab()+ APC.KsAPC_LT()))-
            APC.APC_LT_1()*LT0_d*(APC.APCbl()/(APC.APCbl()+ APC.KsAPC_LT())))*
            time_step;

    /// Cells interact only once with APC and can recieve signaling by CD137 or not.
    LTbo_d+=(APC.APC_LT_1()*LT0_d*(APC.APCa()/(APC.APCa()+ APC.KsAPC_LT()))+
            APC.APC_LT_2()*LT0_d*(APC.APCbo()/(APC.APCbo()+ APC.KsAPC_LT()))+
            LTbo_proliferation_rate_d*LTbo_d-LTbo_apop_rate_d*LTbo_d-
            LTbo_d*u_LT_TNF_d*(m.TNF()/(m.TNF()+Ks_LT_m_TNF_d))-LTbo_d*LT_exh_rate_d)*time_step;


    LTbl_d+=(APC.APC_LT_2()*LT0_d*(APC.APCbo_Ab()/(APC.APCbo_Ab()+ APC.KsAPC_LT()))+
             APC.APC_LT_1()*LT0_d*(APC.APCbl()/(APC.APCbl()+ APC.KsAPC_LT()))+
             LTbl_proliferation_rate_d*LTbl_d-LTbl_apop_rate_d*LTbl_d-
             LTbl_d*u_LT_TNF_d*(m.TNF()/(m.TNF()+Ks_LT_m_TNF_d))-LTbl_d*LT_exh_rate_d)*time_step;

    /// LT get exhausted after a period of time
    LTexh_d+=LTbo_d*LT_exh_rate_d+LTbl_d*LT_exh_rate_d-LTexh_d*LTexh_apop_rate_d-
             LTexh_d*u_LT_TNF_d*(m.TNF()/(m.TNF()+Ks_LT_m_TNF_d))*time_step;

    if (m.TymidineTriteate()>0)
    LT_TymTr_incorporated_d+=(LTns_proliferation_rate_d*LTns_d+LTns_proliferation_rate_d*LT0_d+LTbo_proliferation_rate_d*LTbo_d+
                              +LTbl_proliferation_rate_d*LTbl_d)*m.Prol_TymTr();

    if ((t_run>t_apop_meas_d-t_duration_apoptosis_d)&&(t_run<=t_apop_meas_d))
        Total_cells_in_apoptosis_d+=(LTns_apop_rate_d*LTns_d+LTns_apop_rate_d*LT0_d+
                                     LTbo_apop_rate_d*LTbo_d+LTbo_d*u_LT_TNF_d*(m.TNF()/(m.TNF()+Ks_LT_m_TNF_d))+
                                     LTbl_apop_rate_d*LTbl_d+LTbl_d*u_LT_TNF_d*(m.TNF()/(m.TNF()+Ks_LT_m_TNF_d))+
                                     LTexh_d*u_LT_TNF_d*(m.TNF()/(m.TNF()+Ks_LT_m_TNF_d)))*time_step;
}

/*void LT_cells::reset(const SimParameters& sp,
                      const Treatment& tr)
    {
        LTns_d=sp.ratio_init_LTns_*tr.init_cells;
        LT0_d=sp.ratio_initLTspecific_*tr.init_cells;
        LTbo_d=0;
        LTbl_d=0;
        LTexh_d=0;

    }
    */

/// 1) Total number of LT
double& LT_cells::num_LT()
    {
      double sum=LTns_d+LT0_d+LTbo_d+LTbl_d+LTexh_d;
      return sum;
    }

const double& LT_cells::num_LT() const
    {
      double sum=LTns_d+LT0_d+LTbo_d+LTbl_d+LTexh_d;
      return sum;
    }

/// 2) números de células (5)
double& LT_cells::LTns()
    {
        return LTns_d;
    }

const double& LT_cells::LTns() const
    {
        return LTns_d;
    }

double& LT_cells::LT0()
    {
        return LT0_d;
    }

const double& LT_cells::LT0() const
    {
        return LT0_d;
    }

double& LT_cells::LTbo()
    {
        return LTbo_d;
    }

const double& LT_cells::LTbo() const
    {
        return LTbo_d;
    }

double& LT_cells::LTbl()
    {
        return LTbl_d;
    }

const double& LT_cells::LTbl() const
    {
        return LTbl_d;
    }

double& LT_cells::LTexh()
    {
        return LTexh_d;
    }

const double& LT_cells::LTexh() const
    {
        return LTexh_d;
    }

/// 3) Percentage of cells expressing
double& LT_cells::LT_percentage_cell_expressing_receptor()
{
    double sum=100*(percentage_IFN_LTns_prod_rate_d*(LTns_d+LT0_d)+
               percentage_IFN_LTbo_prod_rate_d*LTbo_d+
                    percentage_IFN_LTbl_prod_rate_d*LTbl_d)/num_LT();
    return sum;
}

const double& LT_cells::LT_percentage_cell_expressing_receptor() const
{
    double sum=100*(percentage_IFN_LTns_prod_rate_d*(LTns_d+LT0_d)+
               percentage_IFN_LTbo_prod_rate_d*LTbo_d+
                    percentage_IFN_LTbl_prod_rate_d*LTbl_d)/num_LT();
    return sum;
}
/// 4) Total cytokines production
double& LT_cells::LT_IFNgamma_production_rate()
   {
      double sum=LTns_d*percentage_IFN_LTns_prod_rate_d*IFN_LTns_prod_rate_d+
               LT0_d*percentage_IFN_LTns_prod_rate_d*IFN_LTns_prod_rate_d+
               LTbo_d*percentage_IFN_LTbo_prod_rate_d*IFN_LTbo_prod_rate_d+
               LTbl_d*percentage_IFN_LTbl_prod_rate_d*IFN_LTbl_prod_rate_d;
      return sum;
  }


const double& LT_cells::LT_IFNgamma_production_rate() const
   {
      double sum=LTns_d*percentage_IFN_LTns_prod_rate_d*IFN_LTns_prod_rate_d+
               LT0_d*percentage_IFN_LTns_prod_rate_d*IFN_LTns_prod_rate_d+
               LTbo_d*percentage_IFN_LTbo_prod_rate_d*IFN_LTbo_prod_rate_d+
               LTbl_d*percentage_IFN_LTbl_prod_rate_d*IFN_LTbl_prod_rate_d;
      return sum;
  }

double& LT_cells::percentage_LT_IFN_production()
   {
      double sum=LTns_d*percentage_IFN_LTns_prod_rate_d+
               LT0_d*percentage_IFN_LTns_prod_rate_d+
               LTbo_d*percentage_IFN_LTbo_prod_rate_d+
               LTbl_d*percentage_IFN_LTbl_prod_rate_d;
      return sum;
  }

const double& LT_cells::percentage_LT_IFN_production() const
   {
      double sum=LTns_d*percentage_IFN_LTns_prod_rate_d+
               LT0_d*percentage_IFN_LTns_prod_rate_d+
               LTbo_d*percentage_IFN_LTbo_prod_rate_d+
               LTbl_d*percentage_IFN_LTbl_prod_rate_d;
      return sum;
  }

double& LT_cells::TNF_production_rate()
   {
      double sum=LTns_d*percentage_TNF_LTns_prod_rate_d*TNF_LTns_prod_rate_d+
                 LT0_d*percentage_TNF_LTns_prod_rate_d*TNF_LTns_prod_rate_d+
                 LTbo_d*percentage_TNF_LTbo_prod_rate_d*TNF_LTbo_prod_rate_d+
                 LTbl_d*percentage_TNF_LTbl_prod_rate_d*TNF_LTbl_prod_rate_d;
      return sum;
  }

const double& LT_cells::TNF_production_rate() const
   {
      double sum=LTns_d*percentage_TNF_LTns_prod_rate_d*TNF_LTns_prod_rate_d+
               LT0_d*percentage_TNF_LTns_prod_rate_d*TNF_LTns_prod_rate_d+
               LTbo_d*percentage_TNF_LTbo_prod_rate_d*TNF_LTbo_prod_rate_d+
               LTbl_d*percentage_TNF_LTbl_prod_rate_d*TNF_LTbl_prod_rate_d;
      return sum;
  }

double& LT_cells::percentage_LT_TNF_production()
   {
      double sum=LTns_d*percentage_TNF_LTns_prod_rate_d+
               LT0_d*percentage_TNF_LTns_prod_rate_d+
               LTbo_d*percentage_TNF_LTbo_prod_rate_d+
               LTbl_d*percentage_TNF_LTbl_prod_rate_d;
      return sum;
  }

const double& LT_cells::percentage_LT_TNF_production() const
   {
      double sum=LTns_d*percentage_TNF_LTns_prod_rate_d+
               LT0_d*percentage_TNF_LTns_prod_rate_d+
               LTbo_d*percentage_TNF_LTbo_prod_rate_d+
               LTbl_d*percentage_TNF_LTbl_prod_rate_d;
      return sum;
  }

/// 5) Tymidine incorporated by LT cells
double& LT_cells::LT_TymTr_incorporated()
{
    return LT_TymTr_incorporated_d;
}
    const double&  LT_cells::LT_TymTr_incorporated()const

{
     return LT_TymTr_incorporated_d;
}


/// 6) Percentage of LT cells undergoing apoptosis
       double& LT_cells::percentage_apoptotic_LT_cells()
       { double sum=100.0*Total_cells_in_apoptosis_d/num_LT();
         return sum;
       }

       const double& LT_cells::percentage_apoptotic_LT_cells() const
       { double sum=100.0*Total_cells_in_apoptosis_d/num_LT();
         return sum;
       }




std::ostream& operator<<(std::ostream& s, const LT_cells& c)
{   
   s<<"\n LT ns \t"<<c.LTns_d;
   s<<"\n LT 0 \t"<<c.LT0_d;
   s<<"\n LT bound \t"<<c.LTbo_d;
   s<<"\n LT blocked \t"<<c.LTbl_d;
   s<<"\n LT exhausted \t"<<c.LTexh_d;
   s<<"\n Tym incorporated by LT \t"<<c.LT_TymTr_incorporated_d;
   s<<"\n Total cells in apoptosis \t"<<c.Total_cells_in_apoptosis_d;
   if (0)
   {
   s<<"///\n----------------------------------\n";
   s<<"those are parameters that do not vary\n\n";

   s<<"\n IFN_LTns_prod_rate_d \t"<<c.IFN_LTns_prod_rate_d;
   s<<"\n IFN_LTbo_prod_rate_d \t"<<c.IFN_LTbo_prod_rate_d;
   s<<"\n IFN_LTbl_prod_rate_d \t"<<c.IFN_LTbl_prod_rate_d;

   s<<"\n TNF_LTns_prod_rate_d \t"<<c.TNF_LTns_prod_rate_d;
   s<<"\n TNF_LTbo_prod_rate_d \t"<<c.TNF_LTbo_prod_rate_d;
   s<<"\n TNF_LTbl_prod_rate_d \t"<<c.TNF_LTbl_prod_rate_d;

   s<<"\n percentage_IFN_LTns_prod_rate_d \t"<<c.percentage_IFN_LTns_prod_rate_d;
   s<<"\n percentage_IFN_LTbo_prod_rate_d \t"<<c. percentage_IFN_LTbo_prod_rate_d;
   s<<"\n percentage_IFN_LTbl_prod_rate_d \t"<<c.percentage_IFN_LTbl_prod_rate_d;

   s<<"\n percentage_TNF_LTns_prod_rate_d \t"<<c.percentage_TNF_LTns_prod_rate_d;
   s<<"\n percentage_TNF_LTbo_prod_rate_d \t"<<c.percentage_TNF_LTbo_prod_rate_d;
   s<<"\n percentage_TNF_LTbl_prod_rate_d \t"<<c.percentage_TNF_LTbl_prod_rate_d;

   s<<"\n LTns_proliferation_rate_d \t"<<c.LTns_proliferation_rate_d;
   s<<"\n LTbo_proliferation_rate_d \t"<<c.LTbo_proliferation_rate_d;
   s<<"\n LTbl_proliferation_rate_d \t"<<c.LTbl_proliferation_rate_d;

   s<<"\n LTns_apop_rate_d \t"<<c.LTns_apop_rate_d;
   s<<"\n LTbo_apop_rate_d \t"<<c.LTbo_apop_rate_d;
   s<<"\n LTbl_apop_rate_d \t"<<c.LTbl_apop_rate_d;
   s<<"\n LTexh_apop_rate_d \t"<<c.LTexh_apop_rate_d;

   s<<"\n Ks_LT_m_TNF_d \t"<<c.Ks_LT_m_TNF_d;

   s<<"\n LTns_expressing_receptor_d \t"<<c.LTns_expressing_receptor_d;

   s<<"\n u_LT_TNF_d \t"<<c.u_LT_TNF_d;

   s<<"\n LT_exh_rate_d \t"<<c.LT_exh_rate_d;

   s<<"\n LT apoptosis measure \t"<<c.t_apop_meas_d;
   s<<"\n LT apoptosis duration \t"<<c.t_duration_apoptosis_d;



}
   return s;
}


