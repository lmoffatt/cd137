#include "Includes/Media.h"
#include "Includes/APC.h"
#include "Includes/LT.h"
#include "Includes/NK.h"


APC_cells::APC_cells(/// 1) init number of APC
                     /*1*/ double init_APC_,

                     /// 2) IFN Poductions rates of each type of APC
                     /*2*/ double IFN_APC0_prod_rate_,
                     /*3*/ double IFN_APCa_prod_rate_,
                     /*4*/ double IFN_APCbo_prod_rate_,


                     /// 3) TNF Poductions rates of each type of APC
                     /*5*/ double TNF_APC0_prod_rate_,
                     /*6*/ double TNF_APCa_prod_rate_,
                     /*7*/ double TNF_APCbo_prod_rate_,


                     /// 4) Percentages of IFN productions of each type of APC
                     /*8*/ double percentage_IFN_APC0_prod_rate_,
                     /*9*/ double percentage_IFN_APCa_prod_rate_,
                     /*10*/ double percentage_IFN_APCbo_prod_rate_,


                     /// 5)Percentages of TNF productions of each type of APC
                     /*11*/ double percentage_TNF_APC0_prod_rate_,
                     /*12*/ double percentage_TNF_APCa_prod_rate_,
                     /*13*/ double percentage_TNF_APCbo_prod_rate_,


                     /// 6) Proliferation rates
                     /*14*/ double APC_bound_proliferation_rate_,

                     /// 7) Apoptosis rates
                     /*15*/ double APC0_apop_rate_,
                     /*16*/ double APCa_apop_rate_,
                     /*17*/ double APCbo_apop_rate_,
                     /*18*/ double APCbl_apop_rate_,
                     /*19*/ double APCexh_apop_rate_,

                     /// 8) constant saturation of TNF for apoptosis
                     /*20*/ double Ks_APC_m_TNF_,

                     /// 9) conversion rates
                     /*21*/ double APC_Ag_,
                     /*22*/ double APC_APC_,
                     /*23*/ double APC_NK_,
                     /*24*/ double APC_LT_1_,
                     /*25*/ double APC_LT_2_,
                     /*26*/ double APC_Ab_,
                     /*27*/ double APC_exh_,

                     /// 10)Saturation constant of IFN and TNF for activation
                     /*28*/ double KsAPC_LT_,

                     /// 11)Saturation constant of APC_LT interaction
                     /*29*/ double Ksi_,
                     /*30*/ double Kst_,

                     /// 12) Percentages of cell expressing receptor
                     /*31*/ double APC0_expressing_receptor_,
                     /*32*/ double APCa_expressing_receptor_,
                     /// 13) Apoptosis rate for TNF
                     /*33*/ double u_APC_TNF_


          ):

          APC0_d(init_APC_),
          APCa_d(0),
          APCbo_d(0),
          APCbo_Ab_d(0),
          APCbl_d (0),
          APCexh_d(0),
          APC_TymTr_incorporated_d(0),

          IFN_APC0_prod_rate_d(IFN_APC0_prod_rate_),
          IFN_APCa_prod_rate_d(IFN_APCa_prod_rate_),
          IFN_APCbo_prod_rate_d(IFN_APCbo_prod_rate_),


          TNF_APC0_prod_rate_d(TNF_APC0_prod_rate_),
          TNF_APCa_prod_rate_d(TNF_APCa_prod_rate_),
          TNF_APCbo_prod_rate_d(TNF_APCbo_prod_rate_),

          percentage_IFN_APC0_prod_rate_d (percentage_IFN_APC0_prod_rate_),
          percentage_IFN_APCa_prod_rate_d (percentage_IFN_APCa_prod_rate_),
          percentage_IFN_APCbo_prod_rate_d (percentage_IFN_APCbo_prod_rate_),

          percentage_TNF_APC0_prod_rate_d (percentage_TNF_APC0_prod_rate_),
          percentage_TNF_APCa_prod_rate_d (percentage_TNF_APCa_prod_rate_),
          percentage_TNF_APCbo_prod_rate_d (percentage_TNF_APCbo_prod_rate_),


          APC_bound_proliferation_rate_d (APC_bound_proliferation_rate_),

          APC0_apop_rate_d(APC0_apop_rate_),
          APCa_apop_rate_d (APCa_apop_rate_),
          APCbo_apop_rate_d (APCbo_apop_rate_),
          APCbl_apop_rate_d (APCbl_apop_rate_),
          APCexh_apop_rate_d (APCexh_apop_rate_),

          Ks_APC_m_TNF_d (Ks_APC_m_TNF_),

          APC_Ag_d (APC_Ag_),
          APC_APC_d(APC_APC_),
          APC_NK_d(APC_NK_),
          APC_LT_1_d(APC_LT_1_),
          APC_LT_2_d(APC_LT_2_),
          APC_Ab_d(APC_Ab_),
          APC_exh_d(APC_exh_),

          KsAPC_LT_d(KsAPC_LT_),

          Ksi_d(Ksi_),
          Kst_d(Kst_),

          APC0_expressing_receptor_d(APC0_expressing_receptor_),
          APCa_expressing_receptor_d(APCa_expressing_receptor_),
          u_APC_TNF_d (u_APC_TNF_)


          {}

// not used yet
/*
APC_cells::APC_cells(const SimParameters& sp,
          const Treatment& tr):

    APC0_d(sp.init_ratio_APC_*tr.init_cells),
    APCa_d(0),
    APCbo_d(0),
    APCbo_Ab_d(0),
    APCbl_d (0),
    APCexh_d(0),
    APC_TymTr_incorporated_d(0),

    IFN_APC0_prod_rate_d(sp.IFN_APC0_prod_rate_),
    IFN_APCa_prod_rate_d(sp.IFN_APCa_prod_rate_),
    IFN_APCbo_prod_rate_d(sp.IFN_APCbo_prod_rate_),


    TNF_APC0_prod_rate_d(sp.TNF_APC0_prod_rate_),
    TNF_APCa_prod_rate_d(sp.TNF_APCa_prod_rate_),
    TNF_APCbo_prod_rate_d(sp.TNF_APCbo_prod_rate_),

    percentage_IFN_APC0_prod_rate_d (sp.percentage_IFN_APC0_prod_rate_),
    percentage_IFN_APCa_prod_rate_d (sp.percentage_IFN_APCa_prod_rate_),
    percentage_IFN_APCbo_prod_rate_d (sp.percentage_IFN_APCbo_prod_rate_),

    percentage_TNF_APC0_prod_rate_d (sp.percentage_TNF_APC0_prod_rate_),
    percentage_TNF_APCa_prod_rate_d (sp.percentage_TNF_APCa_prod_rate_),
    percentage_TNF_APCbo_prod_rate_d (sp.percentage_TNF_APCbo_prod_rate_),


    APC_bound_proliferation_rate_d (sp.APC_bound_proliferation_rate_),

    APC0_apop_rate_d(sp.APC0_apop_rate_),
    APCa_apop_rate_d (sp.APCa_apop_rate_),
    APCbo_apop_rate_d (sp.APCbo_apop_rate_),
    APCbl_apop_rate_d (sp.APCbl_apop_rate_),
    APCexh_apop_rate_d (sp.APCexh_apop_rate_),

    Ks_APC_m_TNF_d (sp.Ks_APC_m_TNF_),

    APC_Ag_d (sp.APC_Ag_),
    APC_APC_d(sp.APC_APC_),
    APC_NK_d(sp.APC_NK_),
    APC_LT_1_d(sp.APC_LT_1_),
    APC_LT_2_d(sp.APC_LT_2_),
    APC_Ab_d(sp.APC_Ab_),
    APC_exh_d(sp.APC_exh_),

    KsAPC_LT_d(sp.KsAPC_LT_),

    Ksi_d(sp.APC_Ksi_),
    Kst_d(sp.APC_Kst_),

    APC0_expressing_receptor_d(sp.APC0_expressing_receptor_),
    APCa_expressing_receptor_d(sp.APCa_expressing_receptor_),
    u_APC_TNF_d (sp.u_APC_TNF_)

          {}
*/
APC_cells::APC_cells(){}

APC_cells::APC_cells(const APC_cells& other):
    APC0_d(other.APC0_d),
    APCa_d(other.APCa_d),
    APCbo_d(other.APCbo_d),
    APCbo_Ab_d(other.APC_Ab_d),
    APCbl_d (other.APCbl_d),
    APCexh_d(other.APCexh_d),
    APC_TymTr_incorporated_d(other.APC_TymTr_incorporated_d),


    IFN_APC0_prod_rate_d(other.IFN_APC0_prod_rate_d),
    IFN_APCa_prod_rate_d(other.IFN_APCa_prod_rate_d),
    IFN_APCbo_prod_rate_d(other.IFN_APCbo_prod_rate_d),


    TNF_APC0_prod_rate_d(other.TNF_APC0_prod_rate_d),
    TNF_APCa_prod_rate_d(other.TNF_APCa_prod_rate_d),
    TNF_APCbo_prod_rate_d(other.TNF_APCbo_prod_rate_d),

    percentage_IFN_APC0_prod_rate_d (other.percentage_IFN_APC0_prod_rate_d),
    percentage_IFN_APCa_prod_rate_d (other.percentage_IFN_APCa_prod_rate_d),
    percentage_IFN_APCbo_prod_rate_d (other.percentage_IFN_APCbo_prod_rate_d),

    percentage_TNF_APC0_prod_rate_d (other.percentage_TNF_APC0_prod_rate_d),
    percentage_TNF_APCa_prod_rate_d (other.percentage_TNF_APCa_prod_rate_d),
    percentage_TNF_APCbo_prod_rate_d (other.percentage_TNF_APCbo_prod_rate_d),


    APC_bound_proliferation_rate_d (other.APC_bound_proliferation_rate_d),

    APC0_apop_rate_d(other.APC0_apop_rate_d),
    APCa_apop_rate_d (other.APCa_apop_rate_d),
    APCbo_apop_rate_d (other.APCbo_apop_rate_d),
    APCbl_apop_rate_d (other.APCbl_apop_rate_d),
    APCexh_apop_rate_d (other.APCexh_apop_rate_d),

    Ks_APC_m_TNF_d (other.Ks_APC_m_TNF_d),

    APC_Ag_d (other.APC_Ag_d),
    APC_APC_d(other.APC_APC_d),
    APC_NK_d(other.APC_NK_d),
    APC_LT_1_d(other.APC_LT_1_d),
    APC_LT_2_d(other.APC_LT_2_d),
    APC_Ab_d(other.APC_Ab_d),
    APC_exh_d(other.APC_exh_d),

    KsAPC_LT_d(other.KsAPC_LT_d),

    Ksi_d(other.Ksi_d),
    Kst_d(other.Kst_d),

    APC0_expressing_receptor_d(other.APC0_expressing_receptor_d),
    APCa_expressing_receptor_d(other.APCa_expressing_receptor_d),
    u_APC_TNF_d (other.u_APC_TNF_d)

    {}


APC_cells&
APC_cells::operator=(const APC_cells& other)
{
    if (this!=&other)
    {
        APC_cells tmp(other);
        swap(*this,tmp);
    }
    return *this;
}

void swap(APC_cells& one, APC_cells& other)
{
   std::swap(one.APC0_d,other.APC0_d);
   std::swap(one.APCa_d,other.APCa_d);
   std::swap(one.APCbo_d,other.APCbo_d);
   std::swap(one.APCbo_Ab_d,other.APCbo_Ab_d);
   std::swap(one.APCbl_d ,other.APCbl_d);
   std::swap(one.APCexh_d,other.APCexh_d);
   std::swap (one.APC_TymTr_incorporated_d,other.APC_TymTr_incorporated_d);

   std::swap(one.IFN_APC0_prod_rate_d,other.IFN_APC0_prod_rate_d);
   std::swap(one.IFN_APCa_prod_rate_d,other.IFN_APCa_prod_rate_d);
   std::swap(one.IFN_APCbo_prod_rate_d,other.IFN_APCbo_prod_rate_d);

   std::swap(one.TNF_APC0_prod_rate_d,other.TNF_APC0_prod_rate_d);
   std::swap(one.TNF_APCa_prod_rate_d,other.TNF_APCa_prod_rate_d);
   std::swap(one.TNF_APCbo_prod_rate_d,other.TNF_APCbo_prod_rate_d);

   std::swap(one.percentage_IFN_APC0_prod_rate_d ,other.percentage_IFN_APC0_prod_rate_d);
   std::swap(one.percentage_IFN_APCa_prod_rate_d ,other.percentage_IFN_APCa_prod_rate_d);
   std::swap(one.percentage_IFN_APCbo_prod_rate_d ,other.percentage_IFN_APCbo_prod_rate_d);

   std::swap(one.percentage_TNF_APC0_prod_rate_d ,other.percentage_TNF_APC0_prod_rate_d);
   std::swap(one.percentage_TNF_APCa_prod_rate_d ,other.percentage_TNF_APCa_prod_rate_d);
   std::swap(one.percentage_TNF_APCbo_prod_rate_d ,other.percentage_TNF_APCbo_prod_rate_d);


   std::swap(one.APC_bound_proliferation_rate_d ,other.APC_bound_proliferation_rate_d);

   std::swap(one.APC0_apop_rate_d,other.APC0_apop_rate_d);
   std::swap(one.APCa_apop_rate_d ,other.APCa_apop_rate_d);
   std::swap(one.APCbo_apop_rate_d ,other.APCbo_apop_rate_d);
   std::swap(one.APCbl_apop_rate_d ,other.APCbl_apop_rate_d);
   std::swap(one.APCexh_apop_rate_d ,other.APCexh_apop_rate_d);

   std::swap(one.Ks_APC_m_TNF_d ,other.Ks_APC_m_TNF_d);

   std::swap(one.APC_Ag_d ,other.APC_Ag_d);
   std::swap(one.APC_APC_d,other.APC_APC_d);
   std::swap(one.APC_NK_d,other.APC_NK_d);
   std::swap(one.APC_LT_1_d,other.APC_LT_1_d);
   std::swap(one.APC_LT_2_d,other.APC_LT_2_d);
   std::swap(one.APC_Ab_d,other.APC_Ab_d);
   std::swap(one.APC_exh_d,other.APC_exh_d);

   std::swap(one.KsAPC_LT_d,other.KsAPC_LT_d);

   std::swap(one.Ksi_d,other.Ksi_d);
   std::swap(one.Kst_d,other.Kst_d);

   std::swap(one.APC0_expressing_receptor_d,other.APC0_expressing_receptor_d);
   std::swap(one.APCa_expressing_receptor_d,other.APCa_expressing_receptor_d);
   std::swap(one.u_APC_TNF_d ,other.u_APC_TNF_d);


}




/// Main step for APC (APC dynamics)
void APC_cells::update(double time_step,const Media& m, const NK_cells& NK, const LT_cells& LT)
{
        /// we update each subpopulation of cells independently and we take into account the transition from one state to the other

    /// the number of free cells (no Ag) they proliferate according to the cell concentration (factor proliferation ratio)
    /// and some of them are "lost" since they internalize the Ag
    APC0_d+=(-APC0_apop_rate_d*APC0_d-
             APC0_d*APC_Ag_d*(m.IFNgamma()/(m.IFNgamma()+ Ksi_d))*(m.TNF()/(m.TNF()+Kst_d))*m.Ag()-
             APC_Ag_d*m.Ag()*APC0_d)*time_step;

    /** the cells that have internalize the Ag proliferate in the same way than the free
    they grow also by the free cells that internalize the Ag
    they shrink by the cells that interact with the LT cells*/
    APCa_d+=(APC0_d*APC_Ag_d*(m.IFNgamma()/(m.IFNgamma()+ Ksi_d))*(m.TNF()/(m.TNF()+Kst_d))*m.Ag() +
            APC_Ag_d*m.Ag()*APC0_d -
            APCa_apop_rate_d*APCa_d -
            u_APC_TNF_d*APCa_d*(m.TNF()/(m.TNF()+ Ks_APC_m_TNF_d))-
            APC_APC_d*APCa_d*APCa_expressing_receptor_d*APCa_d*APCa_expressing_receptor_d-
            APC_APC_d*APCa_d*APCa_expressing_receptor_d*APCbo_d-
            APC_NK_d*APCa_d*APCa_expressing_receptor_d*NK.NKa()*NK.NKa_expressing_receptor()-
            APC_NK_d*APCa_d*APCa_expressing_receptor_d*NK.NKbo()-
            APC_LT_1_d*LT.LT0()*APCa_d*APCa_expressing_receptor_d/(APCa_d*APCa_expressing_receptor_d+KsAPC_LT_d)-
            APC_Ab_d*APCa_d*APCa_expressing_receptor_d*m.Ab()*APCa_d-APCa_d*APC_exh_d)*time_step;


    /// the cells that have interacted with LT grow accordingly with the number of cells that have internalized the Ag and the
    /// number of APC cells expressing the ligand and receptor and bound with monocytes, NK or LT cells with the same rate (Monocytes, NK or LT can
    /// interact only with one cell). We are supposing that all activated cells express receptor and ligand. We are supposing that probabiliities of
    /// interactionts between cells are similar.
    APCbo_d+=(APC_APC_d*APCa_d*APCa_expressing_receptor_d*APCa_d*APCa_expressing_receptor_d+
            APC_APC_d*APCa_d*APCa_expressing_receptor_d*APCbo_d+
            APC_NK_d*APCa_d*APCa_expressing_receptor_d*NK.NKa()*NK.NKa_expressing_receptor()+
            APC_NK_d*APCa_d*APCa_expressing_receptor_d*NK.NKbo()+
            APC_LT_1_d*LT.LT0()*APCa_d*APCa_expressing_receptor_d/(APCa_d*APCa_expressing_receptor_d+KsAPC_LT_d)+
            APCbo_d*APC_bound_proliferation_rate_d-
            APCbo_d*APCbo_apop_rate_d-
            u_APC_TNF_d*APCbo_d*(m.TNF()/(m.TNF()+ Ks_APC_m_TNF_d))-
            APCbo_d*APC_exh_d - APC_Ab_d*APCbo_d*m.Ab())*time_step;


    /// the cells that were signalizeb by receptor and binds the Ab
    APCbo_Ab_d+= (APC_Ab_d*APCbo_d*m.Ab() +
                 APCbo_Ab_d*APC_bound_proliferation_rate_d-
                 APCbo_Ab_d*APCbo_apop_rate_d-
                 u_APC_TNF_d*APCbo_Ab_d*(m.TNF()/(m.TNF()+ Ks_APC_m_TNF_d))-
                 APCbo_Ab_d*APC_exh_d)*time_step;



    /// the cells that are blocked
    APCbl_d+= (APC_Ab_d*APCa_d*APCa_expressing_receptor_d*m.Ab()-
              APCbo_apop_rate_d*APCbl_d-
              u_APC_TNF_d*APCbl_d*(m.TNF()/(m.TNF()+ Ks_APC_m_TNF_d))-
              APCbl_d*APC_exh_d)*time_step;

    /// the cells that are exhausted
    APCexh_d+= (APCa_d*APC_exh_d +
                APCbo_d*APC_exh_d +
                APCbo_Ab_d*APC_exh_d +
                APCbl_d*APC_exh_d -
                APCexh_apop_rate_d-
                u_APC_TNF_d*APCbo_d*(m.TNF()/(m.TNF()+ Ks_APC_m_TNF_d)))*time_step;
    if (m.TymidineTriteate()>0)
        APC_TymTr_incorporated_d+=(APCbo_d+APCbo_Ab_d)*APC_bound_proliferation_rate_d*m.Prol_TymTr();


}
/// Reset number of APC
/*
void APC_cells::reset(const SimParameters& sp,
                      const Treatment& tr)
{
    APC0_d=sp.init_ratio_APC_*tr.init_cells;
    APCa_d=0;
    APCbo_d=0;
    APCbo_Ab_d=0;
    APCbl_d=0;
    APCexh_d=0;
    APC_TymTr_incorporated_d=0;
 }
 */
/// 1) Total number of cells
double& APC_cells::num_APC()
    {
        double sum=APC0_d+APCa_d+APCbo_d+APCbo_Ab_d+APCbl_d+APCexh_d;
        return sum;
    }

const double& APC_cells::num_APC() const
    {
        double sum=APC0_d+APCa_d+APCbo_d+APCbo_Ab_d+APCbl_d+APCexh_d;
        return sum;
    }
/// 2) number of cells (6)
double& APC_cells::APC0()
    {
        return APC0_d;
    }

const double& APC_cells::APC0()const
    {
        return APC0_d;
    }

double& APC_cells::APCa()
    {
        return APCa_d;
    }

const double& APC_cells::APCa()const
    {
        return APCa_d;
    }

double& APC_cells::APCbo()
    {
        return APCbo_d;
    }

const double& APC_cells::APCbo()const
    {
        return APCbo_d;
    }

double& APC_cells::APCbo_Ab()
    {
        return APCbo_d;
    }

const double& APC_cells::APCbo_Ab()const
    {
        return APCbo_Ab_d;
    }


double& APC_cells::APCbl()
    {
        return APCbl_d;
    }

const double& APC_cells::APCbl() const
    {
        return APCbl_d;
    }

double& APC_cells::APCexh()
    {
        return APCexh_d;
    }

const double& APC_cells::APCexh() const
    {
        return APCexh_d;
    }



/// 3) Percentage of cells expressing receptor
double& APC_cells::percentage_cell_expressing_receptor()
    {
       double sum=100*(APC0_d*APC0_expressing_receptor_d+APCa_d*APCa_expressing_receptor_d+APCbo_d+APCbo_Ab_d+APCbl_d)/num_APC();
       return sum;
    }
const double& APC_cells::percentage_cell_expressing_receptor()const
    {
       double sum=100*(APC0_d*APC0_expressing_receptor_d+APCa_d*APCa_expressing_receptor_d+APCbo_d+APCbo_Ab_d+APCbl_d)/num_APC();
       return sum;
    }

/// 4) Cytokines production rate and producing cells (6)

double& APC_cells::APC_IFNgamma_production_rate()
{
        double sum=APC0_d*IFN_APC0_prod_rate_d*percentage_IFN_APC0_prod_rate_d+
                  (APCa_d+APCbl_d)*IFN_APCa_prod_rate_d*percentage_IFN_APCa_prod_rate_d+
                  (APCbo_d+APCbo_Ab_d)*IFN_APCbo_prod_rate_d*percentage_IFN_APCbo_prod_rate_d;
        return sum;
    }
const double& APC_cells::APC_IFNgamma_production_rate() const
    {
        double sum=APC0_d*IFN_APC0_prod_rate_d*percentage_IFN_APC0_prod_rate_d+
                  (APCa_d+APCbl_d)*IFN_APCa_prod_rate_d*percentage_IFN_APCa_prod_rate_d+
                  (APCbo_d+APCbo_Ab_d)*IFN_APCbo_prod_rate_d*percentage_IFN_APCbo_prod_rate_d;
        return sum;
    }

double& APC_cells::APC_TNF_production_rate()
    {
        double sum=APC0_d*TNF_APC0_prod_rate_d*percentage_TNF_APC0_prod_rate_d+
                  (APCa_d+APCbl_d)*TNF_APCa_prod_rate_d*percentage_TNF_APCa_prod_rate_d+
                  (APCbo_d+APCbo_Ab_d)*TNF_APCbo_prod_rate_d*percentage_TNF_APCbo_prod_rate_d;
        return sum;
    }

const double& APC_cells::APC_TNF_production_rate() const
    {
        double sum=APC0_d*TNF_APC0_prod_rate_d*percentage_TNF_APC0_prod_rate_d+
                  (APCa_d+APCbl_d)*TNF_APCa_prod_rate_d*percentage_TNF_APCa_prod_rate_d+
                  (APCbo_d+APCbo_Ab_d)*TNF_APCbo_prod_rate_d*percentage_TNF_APCbo_prod_rate_d;
        return sum;
    }


double& APC_cells::percentage_APC_producing_IFN()
    {
       double sum=100*(percentage_IFN_APC0_prod_rate_d*APC0_d+percentage_IFN_APCa_prod_rate_d*(APCa_d+APCbl_d)+percentage_IFN_APCbo_prod_rate_d*(APCbo_d+APCbo_Ab_d))/num_APC();
       return sum;
    }
const double& APC_cells::percentage_APC_producing_IFN() const
    {
       double sum=100*(percentage_IFN_APC0_prod_rate_d*APC0_d+percentage_IFN_APCa_prod_rate_d*(APCa_d+APCbl_d)+percentage_IFN_APCbo_prod_rate_d*(APCbo_d+APCbo_Ab_d))/num_APC();
       return sum;
    }

double& APC_cells::percentage_APC_producing_TNF()
    {
       double sum=100*(percentage_TNF_APC0_prod_rate_d*APC0_d+percentage_TNF_APCa_prod_rate_d*(APCa_d+APCbl_d)+percentage_TNF_APCbo_prod_rate_d*(APCbo_d+APCbo_Ab_d))/num_APC();
       return sum;
    }
const double& APC_cells::percentage_APC_producing_TNF() const
    {
       double sum=100*(percentage_TNF_APC0_prod_rate_d*APC0_d+percentage_TNF_APCa_prod_rate_d*(APCa_d+APCbl_d)+percentage_TNF_APCbo_prod_rate_d*(APCbo_d+APCbo_Ab_d))/num_APC();
       return sum;
    }

/// APCa TNF production rate
double& APC_cells::APCa_TNF_production_rate()
    {
       return TNF_APCa_prod_rate_d;
    }
const double& APC_cells::APCa_TNF_production_rate() const
    {
       return TNF_APCa_prod_rate_d;
    }

/// APCbo TNF production rate
double& APC_cells::APCbo_TNF_production_rate()
   {
   return TNF_APCbo_prod_rate_d;

   }
const double& APC_cells::APCbo_TNF_production_rate() const
   {
   return TNF_APCa_prod_rate_d;
   }


/// 5) Percentage of activated cells expressing receptor
double& APC_cells::APCa_expressing_receptor ()
   {
        return APCa_expressing_receptor_d;
   }
const double& APC_cells::APCa_expressing_receptor () const
   {
        return APCa_expressing_receptor_d;
   }
/// 6) Union rates of APC (5)
        /// Ag internalization rate

        double& APC_cells::APC_Ag()
        {
             return APC_Ag_d;
        }
        const double& APC_cells::APC_Ag()const
        {
             return APC_Ag_d;
        }

        /// Receptor binding rate to  NK
        const double& APC_cells::APC_NK()const
        {
             return APC_NK_d;
        }
        double& APC_cells::APC_NK()
        {
             return APC_NK_d;
        }

        /// Receptor binding rate to LT
        const double&APC_cells::APC_LT_1() const
        {
             return APC_LT_1_d;
        }
        double& APC_cells::APC_LT_1(){
            return APC_LT_1_d;
       }

        const double&APC_cells:: APC_LT_2() const
        {
             return APC_LT_2_d;
        }
        double& APC_cells::APC_LT_2()
        {
             return APC_LT_2_d;
        }

        /// Ab binding rate
        const double& APC_cells::APC_Ab()const
        {
             return APC_Ab_d;
        }
        double& APC_cells::APC_Ab(){
            return APC_Ab_d;
       }

        /// Constante de saturación Unión APC_LT
        const double& APC_cells::KsAPC_LT () const{
            return KsAPC_LT_d;
       }
        double& APC_cells::KsAPC_LT(){
            return KsAPC_LT_d;
       }

        /// Tymidine incorporated by APC cells
        double& APC_cells::APC_TymTr_incorporated()
        { return APC_TymTr_incorporated_d;}

        const double& APC_cells::APC_TymTr_incorporated()const
         { return APC_TymTr_incorporated_d;}






 std::ostream& operator<<(std::ostream& s, const APC_cells& c)
{
       s<<"\n total APC \t"<<c.num_APC();
       s<<"\n APC0 \t"<<c.APC0_d;
       s<<"\n APCa \t"<<c.APCa_d;
       s<<"\n APCbo \t"<<c.APCbo_d;
       s<<"\n APCbo_Ab \t"<<c.APCbo_Ab_d;
       s<<"\n APCbl \t"<<c.APCbl_d;
       s<<"\n APCexh \t"<<c.APCexh_d;
       s<<"\n APC tymTR incorporate \t"<<c.APC_TymTr_incorporated_d;

       if (0)
       {
       s<<"\n///TNF and INF Poductions rates of each type of APC\n";
       s<<"\n IFN_APC0_prod_rate \t"<<c.IFN_APC0_prod_rate_d;
       s<<"\n IFN_APCa_prod_rate \t"<<c.IFN_APCa_prod_rate_d;
       s<<"\n IFN_APCbo_prod_rate \t"<<c.IFN_APCbo_prod_rate_d;
       s<<"\n TNF_APC0_prod_rate \t"<<c.TNF_APC0_prod_rate_d;
       s<<"\n TNF_APCa_prod_rate \t"<<c.TNF_APCa_prod_rate_d;
       s<<"\n TNF_APCbo_prod_rate \t"<<c.TNF_APCbo_prod_rate_d;


       s<<"\n/// those are parameters that do not vary\n";
       s<<"\n percentage_IFN_APC0_prod_rate_d \t"<<c.percentage_IFN_APC0_prod_rate_d;
       s<<"\n percentage_IFN_APCa_prod_rate_d \t"<<c.percentage_IFN_APCa_prod_rate_d;
       s<<"\n percentage_IFN_APCbo_prod_rate_d \t"<<c.percentage_IFN_APCbo_prod_rate_d;

       s<<"\n percentage_TNF_APC0_prod_rate_d \t"<<c.percentage_TNF_APC0_prod_rate_d;
       s<<"\n percentage_TNF_APCa_prod_rate_d \t"<<c.percentage_TNF_APCa_prod_rate_d;
       s<<"\n percentage_TNF_APCbo_prod_rate_d \t"<<c.percentage_TNF_APCbo_prod_rate_d;


       s<<"\n APC0_apop_rate \t"<<c.APC0_apop_rate_d;
       s<<"\n APCa_apop_rate \t"<<c.APCa_apop_rate_d;
       s<<"\n APCbo_apop_rate \t"<<c.APCbo_apop_rate_d;
       s<<"\n APCbl_apop_rate \t"<<c.APCbl_apop_rate_d;
       s<<"\n  APCexh_apop_rate \t"<<c.APCexh_apop_rate_d;

       s<<"\n APC_bound_proliferation_rate \t"<<c.APC_bound_proliferation_rate_d;
       s<<"\n APC_Ag \t"<<c.APC_Ag_d;
       s<<"\n APC_APC \t"<<c.APC_APC_d;
       s<<"\n APC_NK_d \t"<<c.APC_NK_d;
       s<<"\n APC_LT_1 \t"<<c.APC_LT_1_d;
       s<<"\n APC_LT_1 \t"<<c.APC_LT_2_d;
       s<<"\n APC_Ab \t"<<c.APC_Ab_d;
       s<<"\n APC_exh rate \t"<<c.APC_exh_d;
       s<<"\n Ks_APC_m_TNF_d \t"<<c.Ks_APC_m_TNF_d;
       s<<"\n KsAPC_LT \t"<<c.KsAPC_LT_d;
       s<<"\n Ksi_d \t"<<c.Ksi_d;
       s<<"\n Kst_d \t"<<c.Kst_d;
       s<<"\n APC0_expressing_receptor_d \t"<<c.APC0_expressing_receptor_d;
       s<<"\n APCa_expressing_receptor_d \t"<<c.APCa_expressing_receptor_d;
       s<<"\n u_APC_TNF_d \t"<<c.u_APC_TNF_d;
       }
       return s;


}

APC_cells::APC_cells(const Parameters& p, const Treatment& t):
    /// Variables 7
    /// number of native cells
    APC0_d(
            (1.0-p.mean_ratio("init_K_ratio_LT"))*
            p.mean_ratio("init_K_ratio_APC_NK")*t.init_cells
            ),
    /// number of cells that have internalized the antigen (and therefore express the ligand and receptor)
     APCa_d(0.0),
    /// number of cells that have that have been signaled by receptor or ligand
     APCbo_d(0.0),
    /// number of cells that have that have been signaled by receptor or ligand and bound to the blocking Ab
     APCbo_Ab_d(0.0),
    /// number of cells that binds the blocking mAb
     APCbl_d(0.0),
    /// number of cells that are exhausted
     APCexh_d(0.0),
    /// Timidina incorporada
     APC_TymTr_incorporated_d(0.0),

    /// Parámetros 33
    /// 1) Init ratio of cells
    /*1*/  init_ratio_APC_d(p.mean_ratio("init_K_ratio_APC_NK")),

    /// 2) IFN Poductions rates of each type of APC
    /*2*/  IFN_APC0_prod_rate_d(p.mean("IFN_APC0_prod_rate")),
    /*3*/  IFN_APCa_prod_rate_d(p.mean("IFN_APCa_prod_rate")),
    /*4*/  IFN_APCbo_prod_rate_d(p.mean("IFN_APCbo_prod_rate")),


    /// 3) TNF Poductions rates of each type of APC
    /*5*/  TNF_APC0_prod_rate_d(p.mean("TNF_APC0_prod_rate")),
    /*6*/  TNF_APCa_prod_rate_d(p.mean("TNF_APCa_prod_rate")),
    /*7*/  TNF_APCbo_prod_rate_d(p.mean("TNF_APCbo_prod_rate")),


    /// 4) Percentages of IFN productions of each type of APC
    /*8*/  percentage_IFN_APC0_prod_rate_d(p.mean_ratio("Kpercentage_IFN_APC0_prod_rate")),
    /*9*/  percentage_IFN_APCa_prod_rate_d(p.mean_ratio("Kpercentage_IFN_APCa_prod_rate")),
    /*10*/  percentage_IFN_APCbo_prod_rate_d(p.mean_ratio("Kpercentage_IFN_APCbo_prod_rate")),

    /// 5)Percentages of TNF productions of each type of APC
    /*11*/  percentage_TNF_APC0_prod_rate_d(p.mean_ratio("Kpercentage_TNF_APC0_prod_rate")),
    /*12*/  percentage_TNF_APCa_prod_rate_d(p.mean_ratio("Kpercentage_TNF_APCa_prod_rate")),
    /*13*/  percentage_TNF_APCbo_prod_rate_d(p.mean("Kpercentage_TNF_APCbo_prod_rate")),


    /// 6) Proliferation rates
    /*14*/  APC_bound_proliferation_rate_d(p.mean("APC_bound_proliferation_rate")),

    /// 7) Apoptosis rates
    /*15*/  APC0_apop_rate_d(p.mean("APC0_apop_rate")),
    /*16*/  APCa_apop_rate_d(p.mean("APCa_apop_rate")),
    /*17*/  APCbo_apop_rate_d(p.mean("APCbo_apop_rate")),
    /*18*/  APCbl_apop_rate_d(p.mean("APCbl_apop_rate")),
    /*19*/  APCexh_apop_rate_d(p.mean("APCexh_apop_rate")),

    /// 8) constant saturation of TNF for apoptosis
    /*20*/  Ks_APC_m_TNF_d(p.mean("Ks_APC_m_TNF")),

    /// 9) conversion rates
    /*21*/  APC_Ag_d(p.mean("APC_Ag")),
    /*22*/  APC_APC_d(p.mean("APC_APC")),
    /*23*/  APC_NK_d(p.mean("APC_NK")),
    /*24*/  APC_LT_1_d(p.mean("APC_LT_1")),
    /*25*/  APC_LT_2_d(p.mean("APC_LT_2")),
    /*26*/  APC_Ab_d(p.mean("APC_Ab")),
    /*27*/  APC_exh_d(p.mean("APC_exh")),

    /// 10)Saturation constant of IFN and TNF for activation
    /*28*/  KsAPC_LT_d(p.mean("KsAPC_LT")),

    /// 11)Saturation constant of APC_LT interaction
    /*29*/  Ksi_d(p.mean("APC_Ksi")),
    /*30*/  Kst_d(p.mean("APC_Kst")),

    /// 12) Percentages of cell expressing receptor
    /*31*/  APC0_expressing_receptor_d(p.mean("APC0_Kratio_expressing_receptor")),
    /*32*/  APCa_expressing_receptor_d(p.mean("APCa_Kratio_expressing_receptor")),

    /// 13) Apoptosis rate for TNF
    /*33*/  u_APC_TNF_d(p.mean("u_APC_TNF"))
{}











