#include "Includes/Media.h"
#include "Includes/APC.h"
#include "Includes/LT.h"
#include "Includes/NK.h"
#include <cmath>


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
                   //*11*/ double percentage_IFN_LTbl_prod_rate_,
                   /// 5)Percentages of TNF productions of each type of LT
                   /*12*/ double percentage_TNF_LTns_prod_rate_,
                   /*13*/ double percentage_TNF_LTbo_prod_rate_,
                   //*14*/ double percentage_TNF_LTbl_prod_rate_,
                   /// 6) Proliferation rates
                   /*15*/ double LTns_proliferation_rate_,
                   /*16*/ double LTbo_proliferation_rate_,
                   /*17*/ double LTbl_proliferation_rate_,
                   /// 7) Apoptosis rates
                   /*18*/ double LTns_apop_rate_,
                   /*19*/ double LTbo_apop_rate_,
                   /*20*/ double LTbl_apop_rate_,
                   /// 8) constant saturation of TNF for apoptosis
                   /*21*/ double Ks_LT_m_TNF_,
                   /// 9) Percentages of cell expressing receptor
                   /*22*/ double LTns_expressing_receptor_,
                   /// 10) Apoptosis rate for TNF
                   /*23*/ double u_LT_TNF_,
                   /// 12) apoptosis related parameters
                   /*24*/ double t_apop_meas_,
                   /*25*/ double t_duration_apoptosis_,
                   /*26*/ double LT_Ab_

                   ):
    /*1*/ LTns_d(ratio_init_LTns_),
    /*2*/ LT0_d(ratio_initLTspecific_),
    /*3*/ LTbo_d(0),
    /*4*/ LTbl_d(0),
    /*5*/ LT_TymTr_incorporated_d(0),
    /*6*/ Total_cells_in_apoptosis_d(0),
    /*7*/ IFN_LTns_prod_rate_d(IFN_LTns_prod_rate_),
    /*8*/ IFN_LTbo_prod_rate_d(IFN_LTbo_prod_rate_),
    /*9*/ IFN_LTbl_prod_rate_d(IFN_LTbl_prod_rate_),
    /*10*/ TNF_LTns_prod_rate_d(TNF_LTns_prod_rate_),
    /*11*/ TNF_LTbo_prod_rate_d(TNF_LTbo_prod_rate_),
    /*12*/ TNF_LTbl_prod_rate_d(TNF_LTbl_prod_rate_),
    /*13*/ percentage_IFN_LTns_prod_rate_d(percentage_IFN_LTns_prod_rate_),
    /*14*/ percentage_IFN_LTbo_prod_rate_d( percentage_IFN_LTbo_prod_rate_),
    //*15*/ percentage_IFN_LTbl_prod_rate_d(percentage_IFN_LTbl_prod_rate_),
    /*16*/ percentage_TNF_LTns_prod_rate_d(percentage_TNF_LTns_prod_rate_),
    /*17*/ percentage_TNF_LTbo_prod_rate_d(percentage_TNF_LTbo_prod_rate_),
    //*18*/ percentage_TNF_LTbl_prod_rate_d(percentage_TNF_LTbl_prod_rate_),
    /*19*/ LTns_proliferation_rate_d(LTns_proliferation_rate_),
    /*20*/ LTbo_proliferation_rate_d(LTbo_proliferation_rate_),
    /*21*/ LTbl_proliferation_rate_d(LTbl_proliferation_rate_),
    /*22*/ LTns_apop_rate_d(LTns_apop_rate_),
    /*23*/ LTbo_apop_rate_d(LTbo_apop_rate_),
    /*24*/ LTbl_apop_rate_d(LTbl_apop_rate_),
    /*25*/ Ks_LT_m_TNF_d(Ks_LT_m_TNF_),
    /*26*/ LTns_expressing_receptor_d(LTns_expressing_receptor_),
    /*27*/ u_LT_TNF_d(u_LT_TNF_),
    /*28*/ t_apop_meas_d (t_apop_meas_),
    /*29*/ t_duration_apoptosis_d(t_duration_apoptosis_),
    /*30*/ LT_Ab_d(LT_Ab_)

    {

}




LT_cells::LT_cells(const LT_cells& other):
    /*1*/ LTns_d(other.LTns_d),
    /*2*/ LT0_d(other.LT0_d),
    /*3*/ LTbo_d(other.LTbo_d),
    /*4*/ LTbl_d(other.LTbl_d),
    /*5*/ LT_TymTr_incorporated_d(other.LT_TymTr_incorporated_d),
    /*6*/ Total_cells_in_apoptosis_d(other.Total_cells_in_apoptosis_d),
    /*7*/ IFN_LTns_prod_rate_d(other.IFN_LTns_prod_rate_d),
    /*8*/ IFN_LTbo_prod_rate_d(other.IFN_LTbo_prod_rate_d),
    /*9*/ IFN_LTbl_prod_rate_d(other.IFN_LTbl_prod_rate_d),
    /*10*/ TNF_LTns_prod_rate_d(other.TNF_LTns_prod_rate_d),
    /*11*/ TNF_LTbo_prod_rate_d(other.TNF_LTbo_prod_rate_d),
    /*12*/ TNF_LTbl_prod_rate_d(other.TNF_LTbl_prod_rate_d),
    /*13*/ percentage_IFN_LTns_prod_rate_d(other.percentage_IFN_LTns_prod_rate_d),
    /*14*/ percentage_IFN_LTbo_prod_rate_d(other. percentage_IFN_LTbo_prod_rate_d),
    //*15*/ percentage_IFN_LTbl_prod_rate_d(other.percentage_IFN_LTbl_prod_rate_d),
    /*16*/ percentage_TNF_LTns_prod_rate_d(other.percentage_TNF_LTns_prod_rate_d),
    /*17*/ percentage_TNF_LTbo_prod_rate_d(other.percentage_TNF_LTbo_prod_rate_d),
    //*18*/ percentage_TNF_LTbl_prod_rate_d(other.percentage_TNF_LTbl_prod_rate_d),
    /*19*/ LTns_proliferation_rate_d(other.LTns_proliferation_rate_d),
    /*20*/ LTbo_proliferation_rate_d(other.LTbo_proliferation_rate_d),
    /*21*/ LTbl_proliferation_rate_d(other.LTbl_proliferation_rate_d),
    /*22*/ LTns_apop_rate_d(other.LTns_apop_rate_d),
    /*23*/ LTbo_apop_rate_d(other.LTbo_apop_rate_d),
    /*24*/ LTbl_apop_rate_d(other.LTbl_apop_rate_d),
    /*25*/ Ks_LT_m_TNF_d(other.Ks_LT_m_TNF_d),
    /*26*/ LTns_expressing_receptor_d(other.LTns_expressing_receptor_d),
    /*27*/ u_LT_TNF_d(other.u_LT_TNF_d),
    /*28*/ t_apop_meas_d (other.t_apop_meas_d),
    /*29*/ t_duration_apoptosis_d(other.t_duration_apoptosis_d),
    /*30*/ LT_Ab_d(other.LT_Ab_d)

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
    /*1*/ std::swap(one.LTns_d,other.LTns_d);
    /*2*/ std::swap(one.LT0_d,other.LT0_d);
    /*3*/ std::swap(one.LTbo_d,other.LTbo_d);
    /*4*/ std::swap(one.LTbl_d,other.LTbl_d);
    /*5*/ std::swap(one.LT_TymTr_incorporated_d,other.LT_TymTr_incorporated_d);
    /*6*/ std::swap(one.Total_cells_in_apoptosis_d,other.Total_cells_in_apoptosis_d);
    /*7*/ std::swap(one.IFN_LTns_prod_rate_d,other.IFN_LTns_prod_rate_d);
    /*8*/ std::swap(one.IFN_LTbo_prod_rate_d,other.IFN_LTbo_prod_rate_d);
    /*9*/ std::swap(one.IFN_LTbl_prod_rate_d,other.IFN_LTbl_prod_rate_d);
    /*10*/ std::swap(one.TNF_LTns_prod_rate_d,other.TNF_LTns_prod_rate_d);
    /*11*/ std::swap(one.TNF_LTbo_prod_rate_d,other.TNF_LTbo_prod_rate_d);
    /*12*/ std::swap(one.TNF_LTbl_prod_rate_d,other.TNF_LTbl_prod_rate_d);
    /*13*/ std::swap(one.percentage_IFN_LTns_prod_rate_d,other.percentage_IFN_LTns_prod_rate_d);
    /*14*/ std::swap(one.percentage_IFN_LTbo_prod_rate_d,other. percentage_IFN_LTbo_prod_rate_d);
    //*15*/ std::swap(one.percentage_IFN_LTbl_prod_rate_d,other.percentage_IFN_LTbl_prod_rate_d);
    /*16*/ std::swap(one.percentage_TNF_LTns_prod_rate_d,other.percentage_TNF_LTns_prod_rate_d);
    /*17*/ std::swap(one.percentage_TNF_LTbo_prod_rate_d,other.percentage_TNF_LTbo_prod_rate_d);
    //*18*/ std::swap(one.percentage_TNF_LTbl_prod_rate_d,other.percentage_TNF_LTbl_prod_rate_d);
    /*19*/ std::swap(one.LTns_proliferation_rate_d,other.LTns_proliferation_rate_d);
    /*20*/ std::swap(one.LTbo_proliferation_rate_d,other.LTbo_proliferation_rate_d);
    /*21*/ std::swap(one.LTbl_proliferation_rate_d,other.LTbl_proliferation_rate_d);
    /*22*/ std::swap(one.LTns_apop_rate_d,other.LTns_apop_rate_d);
    /*23*/ std::swap(one.LTbo_apop_rate_d,other.LTbo_apop_rate_d);
    /*24*/ std::swap(one.LTbl_apop_rate_d,other.LTbl_apop_rate_d);
    /*25*/ std::swap(one.Ks_LT_m_TNF_d,other.Ks_LT_m_TNF_d);
    /*26*/ std::swap(one.LTns_expressing_receptor_d,other.LTns_expressing_receptor_d);
    /*27*/ std::swap(one.u_LT_TNF_d,other.u_LT_TNF_d);
    /*28*/ std::swap(one.t_apop_meas_d ,other.t_apop_meas_d);
    /*29*/ std::swap(one.t_duration_apoptosis_d,other.t_duration_apoptosis_d);
    /*30*/ std::swap(one.LT_Ab_d,other.LT_Ab_d);
}



/// Main step for LT
void LT_cells::update(double& time_step, double t_run, const Media& m, const APC_cells& APC, const NK_cells& NK)


{    /// cells not sensitive to the Ag proliferate passively or die
   double LTns_delta=(
               -LTns_apop_rate_d*LTns_d+
               LTns_proliferation_rate_d*LTns_d*m.prol_ratio()
               )*time_step;
   LTns_d+=LTns_delta;


    /// Ag specific cells proliferate and some of them interact with APC and get activated and express the receptor.
    /// During interaction LT receptor can be block or bind to ligand
   double LT0_delta=(LT0_d*(
                         -LTns_apop_rate_d
                         +LTns_proliferation_rate_d*m.prol_ratio()
                         -APC.APC_LT_1()*(APC.APCa())
                         -APC.APC_LT_1()*(APC.APCbl())
                         -APC.APC_LT_1()*(APC.APCbo())
                         -APC.APC_LT_1()*(APC.APCbo_Ab())
                         ))*time_step;
   LT0_d+=LT0_delta;
    /// Cells interact only once with APC and can recieve signaling by CD137 or not.
   double LTbo_delta=(
               APC.APC_LT_1()*LT_Ab_d/(LT_Ab_d+m.Ab())*LT0_d*(APC.APCa())
              +APC.APC_LT_1()*LT_Ab_d/(LT_Ab_d+m.Ab())*LT0_d*(APC.APCbl())
              +APC.APC_LT_1()*LT_Ab_d/(LT_Ab_d+m.Ab())*LT0_d*(APC.APCbo())
              +APC.APC_LT_1()*LT_Ab_d/(LT_Ab_d+m.Ab())*LT0_d*(APC.APCbo_Ab())
              +LTbo_proliferation_rate_d*m.prol_ratio()*LTbo_d
              -LTbo_apop_rate_d*LTbo_d
              -LTbo_d*u_LT_TNF_d*(m.TNF()/(m.TNF()+Ks_LT_m_TNF_d))
               )*time_step;

    LTbo_d+=LTbo_delta;

//        double LTbl_delta=(
//                             APC.APC_LT_1()*m.Ab()/(LT_Ab_d+m.Ab())*LT0_d*(APC.APCa())
//                            +APC.APC_LT_1()*m.Ab()/(LT_Ab_d+m.Ab())*LT0_d*(APC.APCbl())
//                            +APC.APC_LT_1()*m.Ab()/(LT_Ab_d+m.Ab())*LT0_d*(APC.APCbo())
//                            +APC.APC_LT_1()*m.Ab()/(LT_Ab_d+m.Ab())*LT0_d*(APC.APCbo_Ab())
//                            +LTbo_proliferation_rate_d/LTbl_proliferation_rate_d*m.prol_ratio()*LTbl_d
//                            -LTbo_apop_rate_d*LTbl_apop_rate_d*LTbl_d
//                            -LTbl_d*u_LT_TNF_d*(m.TNF()/(m.TNF()+Ks_LT_m_TNF_d))
//                  )*time_step;

    double LTbl_delta=(
                           APC.APC_LT_1()*m.Ab()/(LT_Ab_d+m.Ab())*LT0_d*(APC.APCa())
                          +APC.APC_LT_1()*m.Ab()/(LT_Ab_d+m.Ab())*LT0_d*(APC.APCbl())
                          +APC.APC_LT_1()*m.Ab()/(LT_Ab_d+m.Ab())*LT0_d*(APC.APCbo())
                          +APC.APC_LT_1()*m.Ab()/(LT_Ab_d+m.Ab())*LT0_d*(APC.APCbo_Ab())
                          +LTbo_proliferation_rate_d*LTbl_proliferation_rate_d*m.prol_ratio()*LTbl_d
                          -LTbo_apop_rate_d*LTbl_apop_rate_d*LTbl_d
                          -LTbl_d*u_LT_TNF_d*(m.TNF()/(m.TNF()+Ks_LT_m_TNF_d))
                )*time_step;
    LTbl_d+=LTbl_delta;


   double LT_TymTr_incorporated_delta;
   if (m.TymidineTriteate()>0){
       LT_TymTr_incorporated_delta=(
                   (LTns_proliferation_rate_d*LTns_d*m.prol_ratio()+
                    LTns_proliferation_rate_d*m.prol_ratio()*LT0_d+
                    LTbo_proliferation_rate_d*LTbo_d*m.prol_ratio()
                                  +LTbo_proliferation_rate_d/LTbl_proliferation_rate_d*m.prol_ratio()*LTbl_d)*m.Prol_TymTr()
                   )*time_step;
       LT_TymTr_incorporated_d+=LT_TymTr_incorporated_delta;


   }

//       double LT_TymTr_incorporated_delta;
//       if (m.TymidineTriteate()>0){
//           LT_TymTr_incorporated_delta=(
//                       (LTns_proliferation_rate_d*LTns_d*m.prol_ratio()+
//                        LTns_proliferation_rate_d*m.prol_ratio()*LT0_d+
//                        LTbo_proliferation_rate_d*LTbo_d*m.prol_ratio()
//                                      +LTbo_proliferation_rate_d*m.prol_ratio()*LTbl_d)*m.Prol_TymTr()
//                       )*time_step;
//           LT_TymTr_incorporated_d+=LT_TymTr_incorporated_delta;

//   double  Total_cells_in_apoptosis_delta;
//   if ((t_run>t_apop_meas_d-t_duration_apoptosis_d)&&(t_run<=t_apop_meas_d)){
//        Total_cells_in_apoptosis_delta=(LTns_apop_rate_d*LTns_d+LTns_apop_rate_d*LT0_d+
//                                     LTbo_apop_rate_d*LTbo_d+LTbo_d*u_LT_TNF_d*(m.TNF()/(m.TNF()+Ks_LT_m_TNF_d))+
//                                     LTbo_apop_rate_d*LTbl_apop_rate_d*LTbl_d+LTbl_d*u_LT_TNF_d*(m.TNF()/(m.TNF()+Ks_LT_m_TNF_d)))*time_step;
//        Total_cells_in_apoptosis_d+= Total_cells_in_apoptosis_delta;

        double  Total_cells_in_apoptosis_delta;
        if ((t_run>t_apop_meas_d-t_duration_apoptosis_d)&&(t_run<=t_apop_meas_d)){
             Total_cells_in_apoptosis_delta=(LTns_apop_rate_d*LTns_d+LTns_apop_rate_d*LT0_d+
                                          LTbo_apop_rate_d*LTbo_d+LTbo_d*u_LT_TNF_d*(m.TNF()/(m.TNF()+Ks_LT_m_TNF_d))+
                                          LTbo_apop_rate_d*LTbl_d+LTbl_d*u_LT_TNF_d*(m.TNF()/(m.TNF()+Ks_LT_m_TNF_d)))*time_step;
             Total_cells_in_apoptosis_d+= Total_cells_in_apoptosis_delta;
}
}


/// 1) Total number of LT
double LT_cells::num_LT() const
    {
      double sum=LTns_d+LT0_d+LTbo_d+LTbl_d;
      return sum;
    }


/// 2) números de células (4)
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

/// 3) Percentage of cells expressing
double LT_cells::LT_percentage_cell_expressing_receptor() const
{
    double sum=100.0*(LTns_expressing_receptor_d*(LTns_d+LT0_d)+LTbo_d+LTbl_d)/num_LT();
    return sum;
}

/// 4) Total cytokines production
double LT_cells::LT_IFNgamma_production_rate() const
   {
//      double sum=LTns_d*percentage_IFN_LTns_prod_rate_d*IFN_LTns_prod_rate_d+
//               LT0_d*percentage_IFN_LTns_prod_rate_d*IFN_LTns_prod_rate_d+
//               LTbo_d*percentage_IFN_LTbo_prod_rate_d*IFN_LTbo_prod_rate_d+
//               LTbl_d*percentage_IFN_LTbo_prod_rate_d*IFN_LTbo_prod_rate_d/IFN_LTbl_prod_rate_d;
//      return sum;

      double sum=LTns_d*percentage_IFN_LTns_prod_rate_d*IFN_LTns_prod_rate_d+
               LT0_d*percentage_IFN_LTns_prod_rate_d*IFN_LTns_prod_rate_d+
               LTbo_d*percentage_IFN_LTbo_prod_rate_d*IFN_LTbo_prod_rate_d+
               LTbl_d*percentage_IFN_LTbo_prod_rate_d*IFN_LTbo_prod_rate_d*IFN_LTbl_prod_rate_d;
      return sum;
  }

double LT_cells::percentage_LT_IFN_production() const
   {
      double sum=(LTns_d*percentage_IFN_LTns_prod_rate_d+
               LT0_d*percentage_IFN_LTns_prod_rate_d+
               LTbo_d*percentage_IFN_LTbo_prod_rate_d+
                  LTbl_d*percentage_IFN_LTbo_prod_rate_d)*100/num_LT();
      return sum;
  }


double LT_cells::TNF_production_rate() const
   {
//      double sum=LTns_d*percentage_TNF_LTns_prod_rate_d*TNF_LTns_prod_rate_d+
//                 LT0_d*percentage_TNF_LTns_prod_rate_d*TNF_LTns_prod_rate_d+
//                 LTbo_d*percentage_TNF_LTbo_prod_rate_d*TNF_LTbo_prod_rate_d+
//                 LTbl_d*percentage_TNF_LTbo_prod_rate_d*TNF_LTbo_prod_rate_d/TNF_LTbl_prod_rate_d;
//      return sum;

      double sum=LTns_d*percentage_TNF_LTns_prod_rate_d*TNF_LTns_prod_rate_d+
                 LT0_d*percentage_TNF_LTns_prod_rate_d*TNF_LTns_prod_rate_d+
                 LTbo_d*percentage_TNF_LTbo_prod_rate_d*TNF_LTbo_prod_rate_d+
                 LTbl_d*percentage_TNF_LTbo_prod_rate_d*TNF_LTbo_prod_rate_d*TNF_LTbl_prod_rate_d;
      return sum;
  }


double LT_cells::percentage_LT_TNF_production() const
   {
      double sum=(LTns_d*percentage_TNF_LTns_prod_rate_d+
               LT0_d*percentage_TNF_LTns_prod_rate_d+
               LTbo_d*percentage_TNF_LTbo_prod_rate_d+
                  LTbl_d*percentage_TNF_LTbo_prod_rate_d)*100/num_LT();
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
       double LT_cells::percentage_apoptotic_LT_cells() const
       { double sum=100.0*Total_cells_in_apoptosis_d/num_LT();
         return sum;
       }

       double& LT_cells::LT_Ab()
       {
           return LT_Ab_d;
       }

       const double& LT_cells::LT_Ab() const
       {
           return LT_Ab_d;
       }




std::ostream& operator<<(std::ostream& s, const LT_cells& c)
{   
   /*1*/ s<<"\n LT ns \t"<<c.LTns_d;
   /*2*/ s<<"\n LT 0 \t"<<c.LT0_d;
   /*3*/ s<<"\n LT bound \t"<<c.LTbo_d;
   /*4*/ s<<"\n LT blocked \t"<<c.LTbl_d;
   /*5*/ s<<"\n Tym incorporated by LT \t"<<c.LT_TymTr_incorporated_d;
   /*6*/ s<<"\n Total cells in apoptosis \t"<<c.Total_cells_in_apoptosis_d;
   if (0)
   {
   s<<"///\n----------------------------------\n";
   s<<"those are parameters that do not vary\n\n";

   /*7*/ s<<"\n IFN_LTns_prod_rate_d \t"<<c.IFN_LTns_prod_rate_d;
   /*8*/ s<<"\n IFN_LTbo_prod_rate_d \t"<<c.IFN_LTbo_prod_rate_d;
   /*9*/ s<<"\n IFN_LTbl_prod_rate_d \t"<<c.IFN_LTbl_prod_rate_d;
   /*10*/ s<<"\n TNF_LTns_prod_rate_d \t"<<c.TNF_LTns_prod_rate_d;
   /*11*/ s<<"\n TNF_LTbo_prod_rate_d \t"<<c.TNF_LTbo_prod_rate_d;
   /*12*/ s<<"\n TNF_LTbl_prod_rate_d \t"<<c.TNF_LTbl_prod_rate_d;
   /*13*/ s<<"\n percentage_IFN_LTns_prod_rate_d \t"<<c.percentage_IFN_LTns_prod_rate_d;
   /*14*/ s<<"\n percentage_IFN_LTbo_prod_rate_d \t"<<c. percentage_IFN_LTbo_prod_rate_d;
   //*15*/ s<<"\n percentage_IFN_LTbl_prod_rate_d \t"<<c.percentage_IFN_LTbl_prod_rate_d;
   /*16*/ s<<"\n percentage_TNF_LTns_prod_rate_d \t"<<c.percentage_TNF_LTns_prod_rate_d;
   //*17*/ s<<"\n percentage_TNF_LTbo_prod_rate_d \t"<<c.percentage_TNF_LTbo_prod_rate_d;
   //*18*/ s<<"\n percentage_TNF_LTbl_prod_rate_d \t"<<c.percentage_TNF_LTbl_prod_rate_d;
   /*19*/ s<<"\n LTns_proliferation_rate_d \t"<<c.LTns_proliferation_rate_d;
   /*20*/ s<<"\n LTbo_proliferation_rate_d \t"<<c.LTbo_proliferation_rate_d;
   /*21*/ s<<"\n LTbl_proliferation_rate_d \t"<<c.LTbl_proliferation_rate_d;
   /*22*/ s<<"\n LTns_apop_rate_d \t"<<c.LTns_apop_rate_d;
   /*23*/ s<<"\n LTbo_apop_rate_d \t"<<c.LTbo_apop_rate_d;
   /*24*/ s<<"\n LTbl_apop_rate_d \t"<<c.LTbl_apop_rate_d;
   /*25*/ s<<"\n Ks_LT_m_TNF_d \t"<<c.Ks_LT_m_TNF_d;
   /*26*/ s<<"\n LTns_expressing_receptor_d \t"<<c.LTns_expressing_receptor_d;
   /*27*/ s<<"\n u_LT_TNF_d \t"<<c.u_LT_TNF_d;
   /*28*/ s<<"\n LT apoptosis measure \t"<<c.t_apop_meas_d;
   /*29*/ s<<"\n LT apoptosis duration \t"<<c.t_duration_apoptosis_d;
   /*30*/ s<<"\n LT_Ab \t"<<c.LT_Ab_d;
}
   return s;
}


LT_cells::LT_cells(const Parameters& p, const Treatment& t):
    /// number of non Ag specific cells
    LTns_d(
        p.mean_ratio("init_K_ratio_LT")*
        (1.0-p.mean_ratio("Kratio_initLTspecific"))*
        t.init_cells),

    /// number of naive Ag specific cells
     LT0_d(
        p.mean_ratio("init_K_ratio_LT")*
        p.mean_ratio("Kratio_initLTspecific")*
        t.init_cells
        ),
    /// number of Ag specific cells that have recieve receptor signaling during sinapsis
     LTbo_d(0.0),
    /// number of Ag specific cells that have not recieve receptor singaling during sinapsis
     LTbl_d(0.0),
    /// Tymidine incorporated by APC cells
     LT_TymTr_incorporated_d(0.0),
    /// Total LT cell undergoing apoptosis
     Total_cells_in_apoptosis_d(0.0),

    /// Parameters 25
    /// 1) Init number of LT
        /*1*/  ratio_init_LTns_d(p.mean_ratio("Kratio_initLTspecific")),
        /*2*/  ratio_initLTspecific_d(p.mean_ratio("Kratio_initLTspecific")),
    /// 2) IFN Poductions rates of each type of LT
        /*3*/  IFN_LTns_prod_rate_d(p.mean("IFN_LTns_prod_rate")),
        /*4*/  IFN_LTbo_prod_rate_d(p.mean("IFN_LTbo_prod_rate")),
        /*5*/  IFN_LTbl_prod_rate_d(p.mean("IFN_LTbl_prod_rate")),
    /// 3) TNF Poductions rates of each type of LT
        /*6*/  TNF_LTns_prod_rate_d(p.mean("TNF_LTns_prod_rate")),
        /*7*/  TNF_LTbo_prod_rate_d(p.mean("TNF_LTbo_prod_rate")),
        /*8*/  TNF_LTbl_prod_rate_d(p.mean("TNF_LTbl_prod_rate")),
    /// 4) Percentages of IFN productions of each type of LT
        /*9*/  percentage_IFN_LTns_prod_rate_d(p.mean_ratio("Kpercentage_IFN_LTns_prod_rate")),
        /*10*/  percentage_IFN_LTbo_prod_rate_d(p.mean_ratio("Kpercentage_IFN_LTbo_prod_rate")),
        //*11*/  percentage_IFN_LTbl_prod_rate_d(p.mean_ratio("Kpercentage_IFN_LTbl_prod_rate")),
    /// 5)Percentages of TNF productions of each type of LT
        /*12*/  percentage_TNF_LTns_prod_rate_d(p.mean_ratio("Kpercentage_TNF_LTns_prod_rate")),
        /*13*/  percentage_TNF_LTbo_prod_rate_d(p.mean_ratio("Kpercentage_TNF_LTbo_prod_rate")),
        //*14*/  percentage_TNF_LTbl_prod_rate_d(p.mean_ratio("Kpercentage_TNF_LTbl_prod_rate")),
    /// 6) Proliferation rates
        /*15*/  LTns_proliferation_rate_d(p.mean("LTns_proliferation_rate")),
        /*16*/  LTbo_proliferation_rate_d(p.mean("LTbo_proliferation_rate")),
        /*17*/  LTbl_proliferation_rate_d(p.mean("LTbl_proliferation_rate")),
    /// 7) Apoptosis rates
        /*18*/  LTns_apop_rate_d(p.mean("LTns_apop_rate")),
        /*19*/  LTbo_apop_rate_d(p.mean("LTbo_apop_rate")),
        /*20*/  LTbl_apop_rate_d(p.mean("LTbl_apop_rate")),
    /// 8) constant saturation of TNF for apoptosis
        /*21*/  Ks_LT_m_TNF_d(p.mean("Ks_LT_m_TNF")),
    /// 9) Percentages of cell expressing receptor
        /*22*/  LTns_expressing_receptor_d(p.mean_ratio("LTns_Kratio_expressing_receptor")),
    /// 10) Apoptosis rate for TNF
        /*23*/  u_LT_TNF_d(p.mean("u_LT_TNF")),
    /// 11) LT exh rate
        /*24*/  t_apop_meas_d(t.t_apop_meas_d),
    /// 12) apoptosis related parameters
        /*25*/ t_duration_apoptosis_d(p.mean("t_duration_apoptosis")),
        /*26*/ LT_Ab_d(p.mean("LT_Ab"))



{}


std::vector<double> LT_cells::Derivative(double t_run, const Media& m, const APC_cells& APC)const
{
    std::vector<double> D;

    /// cells not sensitive to the Ag proliferate passively or die
   double LTns_delta=(
               -LTns_apop_rate_d*LTns_d+
               LTns_proliferation_rate_d*m.prol_ratio()*LTns_d
               );

   D.push_back(LTns_delta);


   /// Ag specific cells proliferate and some of them interact with APC and get activated and express the receptor.
   /// During interaction LT receptor can be block or bind to ligand
   double LT0_delta=LT0_d*(
               -LTns_apop_rate_d
               +LTns_proliferation_rate_d*m.prol_ratio()
               -APC.APC_LT_1()*(APC.APCa()/*/(APC.APCa()+APC.KsAPC_LT())*/)
               -APC.APC_LT_1()*(APC.APCbl()/*/(APC.APCbl()+APC.KsAPC_LT())*/)
               -APC.APC_LT_1()*(APC.APCbo()/*/(APC.APCbo()+APC.KsAPC_LT())*/)
               -APC.APC_LT_1()*(APC.APCbo_Ab()/*/(APC.APCbo_Ab()+APC.KsAPC_LT())*/)
               );
   D.push_back(LT0_delta);
    /// Cells interact only once with APC and can recieve signaling by CD137 or not.
   double LTbo_delta=
                 APC.APC_LT_1()*LT_Ab_d/(LT_Ab_d+m.Ab())*LT0_d*(APC.APCa())
                +APC.APC_LT_1()*LT_Ab_d/(LT_Ab_d+m.Ab())*LT0_d*(APC.APCbl())
                +APC.APC_LT_1()*LT_Ab_d/(LT_Ab_d+m.Ab())*LT0_d*(APC.APCbo())
                +APC.APC_LT_1()*LT_Ab_d/(LT_Ab_d+m.Ab())*LT0_d*(APC.APCbo_Ab()/*/(APC.APCbo_Ab()+APC.KsAPC_LT())*/)
                +LTbo_proliferation_rate_d*m.prol_ratio()*LTbo_d
                -LTbo_apop_rate_d*LTbo_d
                -LTbo_d*u_LT_TNF_d*(m.TNF()/(m.TNF()+Ks_LT_m_TNF_d));

    D.push_back(LTbo_delta);
    double LTbl_delta=(
                APC.APC_LT_1()*m.Ab()/(LT_Ab_d+m.Ab())*LT0_d*(APC.APCa()/*/(APC.APCa()+APC.KsAPC_LT())*/)
               +APC.APC_LT_1()*m.Ab()/(LT_Ab_d+m.Ab())*LT0_d*(APC.APCbl()/*/(APC.APCbl()+APC.KsAPC_LT())*/)
               +APC.APC_LT_1()*m.Ab()/(LT_Ab_d+m.Ab())*LT0_d*(APC.APCbo()/*/(APC.APCbo()+APC.KsAPC_LT())*/)
               +APC.APC_LT_1()*m.Ab()/(LT_Ab_d+m.Ab())*LT0_d*(APC.APCbo_Ab()/*/(APC.APCbo_Ab()+APC.KsAPC_LT())*/)
               +LTbl_proliferation_rate_d*m.prol_ratio()*LTbl_d
               -LTbl_apop_rate_d*LTbl_d
               -LTbl_d*u_LT_TNF_d*(m.TNF()/(m.TNF()+Ks_LT_m_TNF_d))
                );
    D.push_back(LTbl_delta);


   double LT_TymTr_incorporated_delta;
   if (m.TymidineTriteate()>0){
       LT_TymTr_incorporated_delta=((LTns_proliferation_rate_d*m.prol_ratio()*LTns_d+
                                     LTns_proliferation_rate_d*m.prol_ratio()*LT0_d+
                                     LTbo_proliferation_rate_d*m.prol_ratio()*LTbo_d+
                                  LTbl_proliferation_rate_d*m.prol_ratio()*LTbl_d)*m.Prol_TymTr());
       }
   else
       LT_TymTr_incorporated_delta=0;

   D.push_back(LT_TymTr_incorporated_delta);
   double  Total_cells_in_apoptosis_delta;
   if ((t_run>t_apop_meas_d-t_duration_apoptosis_d)&&(t_run<=t_apop_meas_d)){
       Total_cells_in_apoptosis_delta=(LTns_apop_rate_d*LTns_d+LTns_apop_rate_d*LT0_d+
                                       LTbo_apop_rate_d*LTbo_d+LTbo_d*u_LT_TNF_d*(m.TNF()/(m.TNF()+Ks_LT_m_TNF_d))+
                                       LTbl_apop_rate_d*LTbl_d+LTbl_d*u_LT_TNF_d*(m.TNF()/(m.TNF()+Ks_LT_m_TNF_d))
                                       );
       }
   else
       Total_cells_in_apoptosis_delta=0;
   D.push_back(Total_cells_in_apoptosis_delta);


   return D;
}


std::vector<double> LT_cells::getState()const
{

    std::vector<double> S;

   S.push_back(LTns_d);
   S.push_back(LT0_d);
   S.push_back(LTbo_d);
   S.push_back(LTbl_d);
   S.push_back(LT_TymTr_incorporated_d);
   S.push_back(Total_cells_in_apoptosis_d);
        return S;

}
void LT_cells::setState(const std::vector<double>& y)
{

   LTns_d=y[0];
   LT0_d=y[1];
   LTbo_d=y[2];
   LTbl_d=y[3];
   LT_TymTr_incorporated_d=y[4];
   Total_cells_in_apoptosis_d=y[5];
}
