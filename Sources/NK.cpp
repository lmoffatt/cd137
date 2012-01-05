#include "Includes/Media.h"
#include "Includes/APC.h"
#include "Includes/NK.h"
#include "Includes/LT.h"
#include <cmath>

NK_cells::NK_cells(/// 1) Init number of NK
                   /*1*/ double init_NK_,

                   /// 2) IFN Poductions rates of each type of NK
                   /*2*/ double IFN_NK0_prod_rate_,
                   /*3*/ double IFN_NKa_prod_rate_,
                   /*4*/ double IFN_NKbo_prod_rate_,


                   /// 3) TNF Poductions rates of each type of NK
                   /*5*/ double TNF_NK0_prod_rate_,
                   /*6*/ double TNF_NKa_prod_rate_,
                   /*7*/ double TNF_NKbo_prod_rate_,


                   /// 4) Percentages of IFN productions of each type of NK
                   /*8*/ double percentage_IFN_NK0_prod_rate_,
                   /*9*/ double percentage_IFN_AgNKa_prod_rate_,
                   /*10*/ double percentage_IFN_NKbo_prod_rate_,


                   /// 5)Percentages of TNF productions of each type of NK
                   /*11*/ double percentage_TNF_NK0_prod_rate_,
                   /*12*/ double percentage_TNF_NKa_prod_rate_,
                   /*13*/ double percentage_TNF_NKbo_prod_rate_,

                   /// 6) Proliferation rates
                   /*13.5*/ double NK0_proliferation_rate_,
                   /*14*/ double NKa_proliferation_rate_,
                   /*15*/ double NKbo_proliferation_rate_,
                   /*16*/ double NKbl_proliferation_rate_,

                   /// 7) Apoptosis rates
                   /*17*/ double NK0_apop_rate_,
                   /*18*/ double NKa_apop_rate_,
                   /*19*/ double NKbo_apop_rate_,
                   /*20*/ double NKbl_apop_rate_,
//                   /*21*/ double NKexh_apop_rate_,

                   /// 8) constant saturation of TNF for apoptosis
                   /*22*/ double Ks_NK_m_TNF_,

                   /// 9) conversion rates
                   /*23*/ double KaNK_,
                   /*24*/ double NK_NK_,
                   /*25*/ double NK_Ab_,
//                   /*26*/ double NK_exh_,

                   /// 10)Saturation constant of APC interaction for activation
                   /*27*/ double KsAPC_NK_,

                   /// 11)Saturation constant of NK_LT interaction
                   /*28*/ double Ksi_,
                   /*29*/ double Kst_,

                   /// 12) Percentages of cell expressing receptor
                   /*30*/ double NK0_expressing_receptor_,
                   /*31*/ double NKa_expressing_receptor_,

                   /// 13) Apoptosis rate for TNF
                   /*32*/ double u_NK_TNF_
          ):

          NK0_d(init_NK_),
          NKa_d(0),
          NKbo_d(0),
          NKbo_Ab_d(0),
          NKbl_d (0),
//          NKexh_d(0),
          NK_TymTr_incorporated_d(0),
    /// 2) IFN Poductions rates of each type of NK
    /*2*/ IFN_NK0_prod_rate_d (IFN_NK0_prod_rate_),
    /*3*/ IFN_NKa_prod_rate_d (IFN_NKa_prod_rate_),
    /*4*/ IFN_NKbo_prod_rate_d (IFN_NKbo_prod_rate_),


    /// 3) TNF Poductions rates of each type of NK
    /*5*/ TNF_NK0_prod_rate_d (TNF_NK0_prod_rate_),
    /*6*/ TNF_NKa_prod_rate_d(TNF_NKa_prod_rate_),
    /*7*/ TNF_NKbo_prod_rate_d (TNF_NKbo_prod_rate_),


    /// 4) Percentages of IFN productions of each type of NK
    /*8*/ percentage_IFN_NK0_prod_rate_d (percentage_IFN_NK0_prod_rate_),
    /*9*/ percentage_IFN_AgNKa_prod_rate_d (percentage_IFN_AgNKa_prod_rate_),
    /*10*/ percentage_IFN_NKbo_prod_rate_d (percentage_IFN_NKbo_prod_rate_),


    /// 5)Percentages of TNF productions of each type of NK
    /*11*/ percentage_TNF_NK0_prod_rate_d (percentage_TNF_NK0_prod_rate_),
    /*12*/ percentage_TNF_NKa_prod_rate_d (percentage_TNF_NKa_prod_rate_),
    /*13*/ percentage_TNF_NKbo_prod_rate_d (percentage_TNF_NKbo_prod_rate_),

    /// 6) Proliferation rates
    /*13.5*/ NK0_proliferation_rate_d (NK0_proliferation_rate_),
    /*14*/ NKa_proliferation_rate_d (NKa_proliferation_rate_),
    /*15*/ NKbo_proliferation_rate_d (NKbo_proliferation_rate_),
    /*16*/ NKbl_proliferation_rate_d (NKbl_proliferation_rate_),

    /// 7) Apoptosis rates
    /*17*/ NK0_apop_rate_d (NK0_apop_rate_),
    /*18*/ NKa_apop_rate_d (NKa_apop_rate_),
    /*19*/ NKbo_apop_rate_d (NKbo_apop_rate_),
    /*20*/ NKbl_apop_rate_d (NKbl_apop_rate_),
//    /*21*/ NKexh_apop_rate_d (NKexh_apop_rate_),

    /// 8) constant saturation of TNF for apoptosis
    /*22*/ Ks_NK_m_TNF_d (Ks_NK_m_TNF_),

    /// 9) conversion rates
    /*23*/ KaNK_d (KaNK_),
    /*24*/ NK_NK_d (NK_NK_),
    /*25*/ NK_Ab_d (NK_Ab_),
//    /*26*/ NK_exh_d (NK_exh_),

    /// 10)Saturation constant of NK interaction for activation
    /*27*/ KsAPC_NK_d (KsAPC_NK_),

    /// 11)Saturation constant of NK_LT interaction
    /*28*/ Ksi_d (Ksi_),
    /*29*/ Kst_d (Kst_),

    /// 12) Percentages of cell expressing receptor
    /*30*/ NK0_expressing_receptor_d (NK0_expressing_receptor_),
    /*31*/ NKa_expressing_receptor_d (NKa_expressing_receptor_),
    /// 13) Apoptosis rate for TNF
    /*32*/ u_NK_TNF_d (u_NK_TNF_)

          {}

// not used yet
/*
NK_cells::NK_cells(const SimParameters& sp,
          const Treatment& tr):

     NK0_d(sp.init_ratio_NK_*tr.init_cells),
     NKa_d(0),
     NKbo_d(0),
     NKbo_Ab_d (0),
     NKbl_d (0),
     NKexh_d(0),
     NK_TymTr_incorporated_d(0),


     IFN_NK0_prod_rate_d (sp.IFN_NK0_prod_rate_),
     IFN_NKa_prod_rate_d (sp.IFN_NKa_prod_rate_),
     IFN_NKbo_prod_rate_d (sp.IFN_NKbo_prod_rate_),


     TNF_NK0_prod_rate_d (sp.TNF_NK0_prod_rate_),
     TNF_NKa_prod_rate_d(sp.TNF_NKa_prod_rate_),
     TNF_NKbo_prod_rate_d (sp.TNF_NKbo_prod_rate_),



     percentage_IFN_NK0_prod_rate_d (sp.percentage_IFN_NK0_prod_rate_),
     percentage_IFN_AgNKa_prod_rate_d (sp.percentage_IFN_AgNKa_prod_rate_),
     percentage_IFN_NKbo_prod_rate_d (sp.percentage_IFN_NKbo_prod_rate_),



     percentage_TNF_NK0_prod_rate_d (sp.percentage_TNF_NK0_prod_rate_),
     percentage_TNF_NKa_prod_rate_d (sp.percentage_TNF_NKa_prod_rate_),
     percentage_TNF_NKbo_prod_rate_d (sp.percentage_TNF_NKbo_prod_rate_),


     NK0_proliferation_rate_d (sp.NK0_proliferation_rate_),
     NKa_proliferation_rate_d (sp.NKa_proliferation_rate_),
     NKbo_proliferation_rate_d (sp.NKbo_proliferation_rate_),
     NKbl_proliferation_rate_d (sp.NKbl_proliferation_rate_),


     NK0_apop_rate_d (sp.NK0_apop_rate_),
     NKa_apop_rate_d (sp.NKa_apop_rate_),
     NKbo_apop_rate_d (sp.NKbo_apop_rate_),
     NKbl_apop_rate_d (sp.NKbl_apop_rate_),
     NKexh_apop_rate_d (sp.NKexh_apop_rate_),


     Ks_NK_m_TNF_d (sp.Ks_NK_m_TNF_),


     KaNK_d (sp.KaNK_),
     NK_NK_d (sp.NK_NK_),
     NK_Ab_d (sp.NK_Ab_),
     NK_exh_d (sp.NK_exh_),

     KsAPC_NK_d (sp.KsAPC_NK_),

     Ksi_d (sp.NK_Ksi_),
     Kst_d (sp.NK_Kst_),

     NK0_expressing_receptor_d (sp.NK0_expressing_receptor_),
     NKa_expressing_receptor_d (sp.NKa_expressing_receptor_),

     u_NK_TNF_d (sp.u_NK_TNF_)

          {}

*/

NK_cells::NK_cells(){}




NK_cells::NK_cells(const NK_cells& other):

   NK0_d(other.NK0_d),
   NKa_d(other.NKa_d),
   NKbo_d(other.NKbo_d),
   NKbo_Ab_d(other.NKbo_Ab_d),
   NKbl_d (other.NKbl_d),
//   NKexh_d(other.NKexh_d),
   NK_TymTr_incorporated_d(other.NK_TymTr_incorporated_d),
   init_ratio_NK_d(other.init_ratio_NK_d),




     IFN_NK0_prod_rate_d(other.IFN_NK0_prod_rate_d),
     IFN_NKa_prod_rate_d(other.IFN_NKa_prod_rate_d),
     IFN_NKbo_prod_rate_d(other.IFN_NKbo_prod_rate_d),

     TNF_NK0_prod_rate_d(other.TNF_NK0_prod_rate_d),
     TNF_NKa_prod_rate_d(other.TNF_NKa_prod_rate_d),
     TNF_NKbo_prod_rate_d(other.TNF_NKbo_prod_rate_d),

     percentage_IFN_NK0_prod_rate_d(other.percentage_IFN_NK0_prod_rate_d),
     percentage_IFN_AgNKa_prod_rate_d(other.percentage_IFN_AgNKa_prod_rate_d),
     percentage_IFN_NKbo_prod_rate_d(other.percentage_IFN_NKbo_prod_rate_d),

     percentage_TNF_NK0_prod_rate_d(other.percentage_TNF_NK0_prod_rate_d),
     percentage_TNF_NKa_prod_rate_d(other.percentage_TNF_NKa_prod_rate_d),
     percentage_TNF_NKbo_prod_rate_d(other.percentage_TNF_NKbo_prod_rate_d),

     NK0_proliferation_rate_d (other.NK0_proliferation_rate_d),
     NKa_proliferation_rate_d(other.NKa_proliferation_rate_d),
     NKbo_proliferation_rate_d(other.NKbo_proliferation_rate_d),
     NKbl_proliferation_rate_d(other.NKbl_proliferation_rate_d),

     NK0_apop_rate_d(other.NK0_apop_rate_d),
     NKa_apop_rate_d(other.NKa_apop_rate_d),
     NKbo_apop_rate_d(other.NKbo_apop_rate_d),
     NKbl_apop_rate_d(other.NKbl_apop_rate_d),
//     NKexh_apop_rate_d(other.NKexh_apop_rate_d),

     Ks_NK_m_TNF_d(other.Ks_NK_m_TNF_d),

     KaNK_d(other.KaNK_d),
     NK_NK_d(other.NK_NK_d),
     NK_Ab_d(other.NK_Ab_d),
//     NK_exh_d(other.NK_exh_d),

     KsAPC_NK_d(other.KsAPC_NK_d),

     Ksi_d(other.Ksi_d),
     Kst_d(other.Kst_d),


     NK0_expressing_receptor_d(other.NK0_expressing_receptor_d),
     NKa_expressing_receptor_d(other.NKa_expressing_receptor_d),


     u_NK_TNF_d(other.u_NK_TNF_d)
    {}

NK_cells&
NK_cells::operator=(const NK_cells& other)
{
    if (this!=&other)
    {
        NK_cells tmp(other);
        swap(*this,tmp);
    }
    return *this;
}

 void swap(NK_cells& one, NK_cells& other)
{

     std::swap(one.NK0_d,other.NK0_d);
     std::swap(one.NKa_d,other.NKa_d);
     std::swap(one.NKbo_d,other.NKbo_d);
     std::swap(one.NKbo_Ab_d,other.NKbo_Ab_d);
     std::swap(one.NKbl_d ,other.NKbl_d);
//     std::swap(one.NKexh_d,other.NKexh_d);
     std::swap (one.NK_TymTr_incorporated_d,other.NK_TymTr_incorporated_d);
     std::swap(one.init_ratio_NK_d,other.init_ratio_NK_d);

     std::swap(one.IFN_NK0_prod_rate_d,other.IFN_NK0_prod_rate_d);
     std::swap(one.IFN_NKa_prod_rate_d,other.IFN_NKa_prod_rate_d);
     std::swap(one.IFN_NKbo_prod_rate_d,other.IFN_NKbo_prod_rate_d);

     std::swap(one.TNF_NK0_prod_rate_d,other.TNF_NK0_prod_rate_d);
     std::swap(one.TNF_NKa_prod_rate_d,other.TNF_NKa_prod_rate_d);
     std::swap(one.TNF_NKbo_prod_rate_d,other.TNF_NKbo_prod_rate_d);

     std::swap(one.percentage_IFN_NK0_prod_rate_d,other.percentage_IFN_NK0_prod_rate_d);
     std::swap(one.percentage_IFN_AgNKa_prod_rate_d,other.percentage_IFN_AgNKa_prod_rate_d);
     std::swap(one.percentage_IFN_NKbo_prod_rate_d,other.percentage_IFN_NKbo_prod_rate_d);

     std::swap(one.percentage_TNF_NK0_prod_rate_d,other.percentage_TNF_NK0_prod_rate_d);
     std::swap(one.percentage_TNF_NKa_prod_rate_d,other.percentage_TNF_NKa_prod_rate_d);
     std::swap(one.percentage_TNF_NKbo_prod_rate_d,other.percentage_TNF_NKbo_prod_rate_d);

     std::swap(one.NK0_proliferation_rate_d ,other.NK0_proliferation_rate_d),
     std::swap(one.NKa_proliferation_rate_d,other.NKa_proliferation_rate_d);
     std::swap(one.NKbo_proliferation_rate_d,other.NKbo_proliferation_rate_d);
     std::swap(one.NKbl_proliferation_rate_d,other.NKbl_proliferation_rate_d);

     std::swap(one.NK0_apop_rate_d,other.NK0_apop_rate_d);
     std::swap(one.NKa_apop_rate_d,other.NKa_apop_rate_d);
     std::swap(one.NKbo_apop_rate_d,other.NKbo_apop_rate_d);
     std::swap(one.NKbl_apop_rate_d,other.NKbl_apop_rate_d);
//     std::swap(one.NKexh_apop_rate_d,other.NKexh_apop_rate_d);

     std::swap(one.Ks_NK_m_TNF_d,other.Ks_NK_m_TNF_d);

     std::swap(one.KaNK_d,other.KaNK_d);
     std::swap(one.NK_NK_d,other.NK_NK_d);
     std::swap(one.NK_Ab_d,other.NK_Ab_d);
//     std::swap(one.NK_exh_d,other.NK_exh_d);


     std::swap(one.KsAPC_NK_d,other.KsAPC_NK_d);

     std::swap(one.Ksi_d,other.Ksi_d);
     std::swap(one.Kst_d,other.Kst_d);

     std::swap(one.NK0_expressing_receptor_d,other.NK0_expressing_receptor_d);
     std::swap(one.NKa_expressing_receptor_d,other.NKa_expressing_receptor_d);

     std::swap(one.u_NK_TNF_d,other.u_NK_TNF_d);

}



/// main step for the NK cells(NK cell dynamics)
void NK_cells::update(double& time_step,const Media& m, const APC_cells& APC,const LT_cells& LT)
{
    /// we update each subpopulation of cells independently and we take into account the transition from one state to the other

    /// the number of free cells (no Ag) they proliferate according and some of them are "lost" since they activate

    double NK0_delta=(NK0_d*NK0_proliferation_rate_d-NK0_d*NK0_apop_rate_d-
            NK0_d*m.Ag()*KaNK_d*((APC.APCa()+(APC.APCbo_TNF_production_rate()/APC.APCa_TNF_production_rate())+APC.APCbl())/
                                (APC.APCa()+(APC.APCbo_TNF_production_rate()/APC.APCa_TNF_production_rate())+APC.APCbl()+KsAPC_NK_d)))*time_step;
    NK0_d+=NK0_delta;
    /** activated NK cells dynamics*/
    double NKa_delta=(NKa_proliferation_rate_d*NKa_d+
            NK0_d*m.Ag()*KaNK_d*(APC.APCa()+ APC.APCbl()+ APC.APCbo_Ab() + (APC.APCbo_TNF_production_rate()/APC.APCa_TNF_production_rate())*APC.APCbl())/
            (APC.APCa()+ APC.APCbl()+ APC.APCbo_Ab() + (APC.APCbo_TNF_production_rate()/APC.APCa_TNF_production_rate())*APC.APCbl()+KsAPC_NK_d)-
           NK_NK_d*NKa_d*NKa_expressing_receptor_d*NKa_d*NKa_expressing_receptor_d-
           NK_NK_d*NKa_d*NKa_expressing_receptor_d*NKbo_d-
           APC.APC_NK()*NKa_d*NKa_expressing_receptor_d*APC.APCa()*APC.APCa_expressing_receptor()-
           APC.APC_NK()*NKa_d*NKa_expressing_receptor_d*APC.APCbo()-
           NK_Ab_d*NKa_d*NKa_expressing_receptor_d*m.Ab()-NKa_apop_rate_d*NKa_d-
           u_NK_TNF_d*NKa_d*(m.TNF()/(m.TNF()+Ks_NK_m_TNF_d))/*-NKa_d*NK_exh_d*/)*time_step;
    NKa_d+=NKa_delta;
/// El porcentaje de células productoras de IL-12 está dentro de la constante

    /// the cells that have interacted with LT grow accordingly with the number of cells that have internalized the Ag and the
    /// number of LT cells, monocytes and NK expressing the receptor with the same affinity (Monocytes, NK or LT can
    /// interact only with one cell). We are supposing that all activated cells express receptor and ligand.
    double NKbo_delta= (NKbo_d*NKbo_proliferation_rate_d+
            NK_NK_d*NKa_d*NKa_expressing_receptor_d*NKa_d*NKa_expressing_receptor_d+
            NK_NK_d*NKa_d*NKa_expressing_receptor_d*NKbo_d+
            APC.APC_NK()*NKa_d*NKa_expressing_receptor_d*APC.APCa()*APC.APCa_expressing_receptor()+
            APC.APC_NK()*NKa_d*NKa_expressing_receptor_d*APC.APCbo()-
            NKbo_apop_rate_d*NKbo_d-u_NK_TNF_d*NKbo_d*(m.TNF()/(m.TNF()+Ks_NK_m_TNF_d))-
            NKbo_d*NK_Ab_d*m.Ab()/*-NKbo_d*NK_exh_d*/)*time_step;

     NKbo_d+=NKbo_delta;

    /// NK cells can interact with other cells or with the blocking mAb

    double NKbo_Ab_delta=(NKbo_Ab_d*NKbo_proliferation_rate_d+NKbo_d*NK_Ab_d*m.Ab()-
            NKbo_apop_rate_d*NKbo_Ab_d-u_NK_TNF_d*NKbo_Ab_d*(m.TNF()/(m.TNF()+Ks_NK_m_TNF_d))/*-
            NKbo_Ab_d*NK_exh_d*/)*time_step;
    NKbo_Ab_d+=NKbo_Ab_delta;
    /// the cells that have interacted with LT get exhausted
    double NKbl_delta=(NKbl_d*NKbl_proliferation_rate_d+NK_Ab_d*NKa_d*NKa_expressing_receptor_d*m.Ab()-
            NKbl_apop_rate_d*NKbl_d-u_NK_TNF_d*NKbl_d*(m.TNF()/(m.TNF()+Ks_NK_m_TNF_d))/*-NKbl_d*NK_exh_d*/)*time_step;

   NKbl_d+=NKbl_delta;
//   double  NKexh_delta=(NKa_d*NK_exh_d+NKbo_d*NK_exh_d+NKbo_Ab_d*NK_exh_d+NKbl_d*NK_exh_d-
//             NKexh_apop_rate_d*NKexh_d-u_NK_TNF_d*NKexh_d*(m.TNF()/(m.TNF()+Ks_NK_m_TNF_d)))*time_step;
//    NKexh_d+=NKexh_delta;
    double NK_TymTr_incorporated_delta;
    if (m.TymidineTriteate()>0){
        NK_TymTr_incorporated_delta=
                ((NK0_d*NK0_proliferation_rate_d+(NKa_d+NKbl_d)*NKa_proliferation_rate_d+
                                  (NKbo_d+NKbo_Ab_d)*NKbo_proliferation_rate_d)*m.Prol_TymTr()
                 )*time_step;
        NK_TymTr_incorporated_d+=NK_TymTr_incorporated_delta;
    }

}
/*

void NK_cells::reset(const SimParameters& sp,
                      const Treatment& tr)
{
    NK0_d=sp.init_ratio_NK_*tr.init_cells;
    NKa_d=0;
    NKbo_d=0;
    NKbo_Ab_d=0;
    NKbl_d=0;
    NKexh_d=0;
    NK_TymTr_incorporated_d=0;
 }
*/
double NK_cells::num_NK() const
    {
        double sum=NK0_d+NKa_d+NKbo_d+NKbo_Ab_d+NKbl_d;
        return sum;
    }


double NK_cells::NK_IFNgamma_production_rate() const
    {
        double sum=IFN_NK0_prod_rate_d*NK0_d*percentage_IFN_NK0_prod_rate_d+
                   IFN_NKa_prod_rate_d*(NKa_d+NKbl_d)*percentage_IFN_AgNKa_prod_rate_d+
                   IFN_NKbo_prod_rate_d*(NKbo_d+NKbo_Ab_d)*percentage_IFN_NKbo_prod_rate_d;
        return sum;
    }



double NK_cells::NK_TNF_production_rate() const
    {
        double sum=TNF_NK0_prod_rate_d*NK0_d*percentage_TNF_NK0_prod_rate_d+
                   TNF_NKa_prod_rate_d*(NKa_d+NKbl_d)*percentage_TNF_NKa_prod_rate_d+
                   TNF_NKbo_prod_rate_d*(NKbo_d+NKbo_Ab_d)*percentage_TNF_NKbo_prod_rate_d;
        return sum;
    }



double NK_cells::percentage_NK_expressing_receptor() const
    {
    double sum=(NK0_expressing_receptor_d*NK0_d+NKa_expressing_receptor_d+NKbo_d+NKbo_Ab_d+NKbl_d)*100/num_NK();
    return sum;
    }


double NK_cells::percentage_NK_producing_IFN() const
   {
    double sum=100*(percentage_IFN_NK0_prod_rate_d*NK0_d
                    +percentage_IFN_AgNKa_prod_rate_d*(NKa_d+NKbl_d)
                    +percentage_IFN_NKbo_prod_rate_d*(NKbo_d+NKbl_d))/num_NK();
    return sum;
   }



double NK_cells::percentage_NK_producing_TNF() const
   {
    double sum=100*(percentage_TNF_NK0_prod_rate_d*NK0_d
                    +percentage_TNF_NKa_prod_rate_d*(NKa_d+NKbl_d)
                    +percentage_TNF_NKbo_prod_rate_d*(NKbo_d+NKbl_d))/num_NK();
    return sum;
   }


double& NK_cells::NK0()
    {
        return NK0_d;
    }

const double& NK_cells::NK0() const
    {
        return NK0_d;
    }

double& NK_cells::NKa()
    {
        return NKa_d;
    }

const double& NK_cells::NKa()const
    {
        return NKa_d;
    }


double& NK_cells::NKbo()
    {
        return NKbo_d;
    }

const double& NK_cells::NKbo() const
    {
        return NKbo_d;
    }

double& NK_cells::NKbo_Ab()
    {
        return NKbo_Ab_d;
    }

const double& NK_cells::NKbo_Ab()const
    {
        return NKbo_Ab_d;
    }


double& NK_cells::NKbl()
    {
        return NKbl_d;
    }

const double& NK_cells::NKbl() const
    {
        return NKbl_d;
    }

//double& NK_cells::NKexh()
//    {
//        return NKexh_d;
//    }

//const double& NK_cells::NKexh()const
//    {
//        return NKexh_d;
//    }

double& NK_cells::NK_TymTr_incorporated()
{ return NK_TymTr_incorporated_d;}

const double& NK_cells::NK_TymTr_incorporated()const
 { return NK_TymTr_incorporated_d;}

double& NK_cells::NKa_expressing_receptor()
{ return NKa_expressing_receptor_d;}

const double& NK_cells::NKa_expressing_receptor() const
{ return NKa_expressing_receptor_d;}



 std::ostream& operator<<(std::ostream& s, const NK_cells& c)
{

    s<<"\n NK0 \t"<<c.NK0_d;
    s<<"\n NKa \t"<<c.NKa_d;
    s<<"\n NKbo_d \t"<<c.NKbo_d;
    s<<"\n NKbo_Ab \t"<<c.NKbo_Ab_d;
    s<<"\n NKbl \t"<<c.NKbl_d;
//    s<<"\n NKexh \t"<<c.NKexh_d;
    s<<"\n NK tymTR incorporate \t"<<c.NK_TymTr_incorporated_d;


    if (0)
    {

       s<<"\n init_ratio_NK_ \t"<<c.init_ratio_NK_d;


       s<<"\n IFN_NK0_prod_rate_ \t"<<c.IFN_NK0_prod_rate_d;
       s<<"\n IFN_NK0_prod_rate_ \t"<<c.IFN_NK0_prod_rate_d;
       s<<"\n IFN_NKbo_prod_rate_ \t"<<c.IFN_NKbo_prod_rate_d;

       s<<"\n TNF_NK0_prod_rate_ \t"<<c.TNF_NK0_prod_rate_d;
       s<<"\n TNF_NKa_prod_rate_ \t"<<c.TNF_NKa_prod_rate_d;
       s<<"\n TNF_NKbo_prod_rate_ \t"<<c.TNF_NKbo_prod_rate_d;

       s<<"\n percentage_IFN_NK0_prod_rate_ \t"<<c.percentage_IFN_NK0_prod_rate_d;
       s<<"\n percentage_IFN_AgNKa_prod_rate_ \t"<<c.percentage_IFN_AgNKa_prod_rate_d;
       s<<"\n percentage_IFN_NKbo_prod_rate_ \t"<<c.percentage_IFN_NKbo_prod_rate_d;

       s<<"\n percentage_TNF_NK0_prod_rate_ \t"<<c.percentage_TNF_NK0_prod_rate_d;
       s<<"\n percentage_TNF_NKa_prod_rate_ \t"<<c.percentage_TNF_NKa_prod_rate_d;
       s<<"\n percentage_TNF_NKbo_prod_rate_ \t"<<c.percentage_TNF_NKbo_prod_rate_d;

       s<<"\n NK0_proliferation_rate_d \t"<<c.NK0_proliferation_rate_d;
       s<<"\n NKa_proliferation_rate_ \t"<<c.NKa_proliferation_rate_d;
       s<<"\n NKbo_proliferation_rate_ \t"<<c.NKbo_proliferation_rate_d;
       s<<"\n NKbl_proliferation_rate_ \t"<<c.NKbl_proliferation_rate_d;

       s<<"\n NK0_apop_rate_ \t"<<c.NK0_apop_rate_d;
       s<<"\n NKa_apop_rate_ \t"<<c.NKa_apop_rate_d;
       s<<"\n NKbo_apop_rate_ \t"<<c.NKbo_apop_rate_d;
       s<<"\n NKbl_apop_rate_ \t"<<c.NKbl_apop_rate_d;
//       s<<"\n NKexh_apop_rate_ \t"<<c.NKexh_apop_rate_d;

       s<<"\n Ks_NK_m_TNF_ \t"<<c.Ks_NK_m_TNF_d;

       s<<"\n KaNK_ \t"<<c.KaNK_d;
       s<<"\n NK_NK_ \t"<<c.NK_NK_d;
       s<<"\n NK_Ab_ \t"<<c.NK_Ab_d;
//       s<<"\n NK_exh_ \t"<<c.NK_exh_d;

       s<<"\n KsAPC_NK_ \t"<<c.KsAPC_NK_d;

       s<<"\n Ksi_ \t"<<c.Ksi_d;
       s<<"\n Kst_ \t"<<c.Kst_d;

       s<<"\n NK0_expressing_receptor_ \t"<<c.NK0_expressing_receptor_d;
       s<<"\n NKa_expressing_receptor_ \t"<<c.NKa_expressing_receptor_d;

       s<<"\n u_NK_TNF_ \t"<<c.u_NK_TNF_d;


    }
    return s;
}

 NK_cells::NK_cells(const Parameters& p, const Treatment& t):
 /// Variables
     NK0_d(
         (1.0-p.mean_ratio("init_K_ratio_LT"))*
         (1.0-p.mean_ratio("init_K_ratio_APC_NK"))
         *t.init_cells
         ),
     NKa_d(0.0),
     NKbo_d(0.0),
     NKbo_Ab_d(0.0),
     NKbl_d (0.0),
//     NKexh_d(0.0),
     NK_TymTr_incorporated_d(0.0),
/// 2) IFN Poductions rates of each type of NK
/*2*/ IFN_NK0_prod_rate_d (p.mean("IFN_NK0_prod_rate")),
/*3*/ IFN_NKa_prod_rate_d (p.mean("IFN_NKa_prod_rate")),
/*4*/ IFN_NKbo_prod_rate_d (p.mean("IFN_NKbo_prod_rate")),


/// 3) TNF Poductions rates of each type of NK
/*5*/ TNF_NK0_prod_rate_d (p.mean("TNF_NK0_prod_rate")),
/*6*/ TNF_NKa_prod_rate_d(p.mean("TNF_NKa_prod_rate")),
/*7*/ TNF_NKbo_prod_rate_d (p.mean("TNF_NKbo_prod_rate")),


/// 4) Percentages of IFN productions of each type of NK
/*8*/ percentage_IFN_NK0_prod_rate_d (p.mean_ratio("Kpercentage_IFN_NK0_prod_rate")),
/*9*/ percentage_IFN_AgNKa_prod_rate_d (p.mean_ratio("Kpercentage_IFN_AgNKa_prod_rate")),
/*10*/ percentage_IFN_NKbo_prod_rate_d (p.mean_ratio("Kpercentage_IFN_NKbo_prod_rate")),


/// 5)Percentages of TNF productions of each type of NK
/*11*/ percentage_TNF_NK0_prod_rate_d (p.mean_ratio("Kpercentage_TNF_NK0_prod_rate")),
/*12*/ percentage_TNF_NKa_prod_rate_d (p.mean_ratio("Kpercentage_TNF_NKa_prod_rate")),
/*13*/ percentage_TNF_NKbo_prod_rate_d (p.mean_ratio("Kpercentage_TNF_NKbo_prod_rate")),

/// 6) Proliferation rates
/*13.5*/ NK0_proliferation_rate_d (p.mean("NK0_proliferation_rate")),
/*14*/ NKa_proliferation_rate_d (p.mean("NKa_proliferation_rate")),
/*15*/ NKbo_proliferation_rate_d (p.mean("NKbo_proliferation_rate")),
/*16*/ NKbl_proliferation_rate_d (p.mean("NKbl_proliferation_rate")),

/// 7) Apoptosis rates
/*17*/ NK0_apop_rate_d (p.mean("NK0_apop_rate")),
/*18*/ NKa_apop_rate_d (p.mean("NKa_apop_rate")),
/*19*/ NKbo_apop_rate_d (p.mean("NKbo_apop_rate")),
/*20*/ NKbl_apop_rate_d (p.mean("NKbl_apop_rate")),
///*21*/ NKexh_apop_rate_d (p.mean("NKexh_apop_rate")),

/// 8) constant saturation of TNF for apoptosis
/*22*/ Ks_NK_m_TNF_d (p.mean("Ks_NK_m_TNF")),

/// 9) conversion rates
/*23*/ KaNK_d (p.mean("KaNK")),
/*24*/ NK_NK_d (p.mean("NK_NK")),
/*25*/ NK_Ab_d (p.mean("NK_Ab")),
///*26*/ NK_exh_d (p.mean("NK_exh")),

/// 10)Saturation constant of NK interaction for activation
/*27*/ KsAPC_NK_d (p.mean("KsAPC_NK")),

/// 11)Saturation constant of NK_LT interaction
/*28*/ Ksi_d (p.mean("NK_Ksi")),
/*29*/ Kst_d (p.mean("NK_Kst")),

/// 12) Percentages of cell expressing receptor
/*30*/ NK0_expressing_receptor_d (p.mean_ratio("NK0_Kratio_expressing_receptor")),
/*31*/ NKa_expressing_receptor_d (p.mean_ratio("NKa_Kratio_expressing_receptor")),
/// 13) Apoptosis rate for TNF
/*32*/ u_NK_TNF_d (p.mean("u_NK_TNF"))
 {

 }




 std::vector<double> NK_cells::Derivative(const Media& m, const APC_cells& APC)const
 {

std::vector<double> D;
     /// we update each subpopulation of cells independently and we take into account the transition from one state to the other

     /// the number of free cells (no Ag) they proliferate according and some of them are "lost" since they activate

     double NK0_delta=(NK0_d*NK0_proliferation_rate_d-NK0_d*NK0_apop_rate_d-
             NK0_d*m.Ag()*KaNK_d*((APC.APCa()+(APC.APCbo_TNF_production_rate()/APC.APCa_TNF_production_rate())+APC.APCbl())/
                                 (APC.APCa()+(APC.APCbo_TNF_production_rate()/APC.APCa_TNF_production_rate())+APC.APCbl()+KsAPC_NK_d)));
     D.push_back(NK0_delta);
     /** activated NK cells dynamics*/
     double NKa_delta=(NKa_proliferation_rate_d*NKa_d+
             NK0_d*m.Ag()*KaNK_d*(APC.APCa()+ APC.APCbl()+ APC.APCbo_Ab() + (APC.APCbo_TNF_production_rate()/APC.APCa_TNF_production_rate())*APC.APCbl())/
             (APC.APCa()+ APC.APCbl()+ APC.APCbo_Ab() + (APC.APCbo_TNF_production_rate()/APC.APCa_TNF_production_rate())*APC.APCbl()+KsAPC_NK_d)-
            NK_NK_d*NKa_d*NKa_expressing_receptor_d*NKa_d*NKa_expressing_receptor_d-
            NK_NK_d*NKa_d*NKa_expressing_receptor_d*NKbo_d-
            APC.APC_NK()*NKa_d*NKa_expressing_receptor_d*APC.APCa()*APC.APCa_expressing_receptor()-
            APC.APC_NK()*NKa_d*NKa_expressing_receptor_d*APC.APCbo()-
            NK_Ab_d*NKa_d*NKa_expressing_receptor_d*m.Ab()-NKa_apop_rate_d*NKa_d-
            u_NK_TNF_d*NKa_d*(m.TNF()/(m.TNF()+Ks_NK_m_TNF_d))/*-NKa_d*NK_exh_d*/);
     D.push_back(NKa_delta);
 /// El porcentaje de células productoras de IL-12 está dentro de la constante

     /// the cells that have interacted with LT grow accordingly with the number of cells that have internalized the Ag and the
     /// number of LT cells, monocytes and NK expressing the receptor with the same affinity (Monocytes, NK or LT can
     /// interact only with one cell). We are supposing that all activated cells express receptor and ligand.
     double NKbo_delta= (NKbo_d*NKbo_proliferation_rate_d+
             NK_NK_d*NKa_d*NKa_expressing_receptor_d*NKa_d*NKa_expressing_receptor_d+
             NK_NK_d*NKa_d*NKa_expressing_receptor_d*NKbo_d+
             APC.APC_NK()*NKa_d*NKa_expressing_receptor_d*APC.APCa()*APC.APCa_expressing_receptor()+
             APC.APC_NK()*NKa_d*NKa_expressing_receptor_d*APC.APCbo()-
             NKbo_apop_rate_d*NKbo_d-u_NK_TNF_d*NKbo_d*(m.TNF()/(m.TNF()+Ks_NK_m_TNF_d))-
             NKbo_d*NK_Ab_d*m.Ab()/*-NKbo_d*NK_exh_d*/);

      D.push_back(NKbo_delta);

     /// NK cells can interact with other cells or with the blocking mAb

     double NKbo_Ab_delta=(NKbo_Ab_d*NKbo_proliferation_rate_d+NKbo_d*NK_Ab_d*m.Ab()-
             NKbo_apop_rate_d*NKbo_Ab_d-u_NK_TNF_d*NKbo_Ab_d*(m.TNF()/(m.TNF()+Ks_NK_m_TNF_d))/*-
             NKbo_Ab_d*NK_exh_d*/);
     D.push_back(NKbo_Ab_delta);
     /// the cells that have interacted with LT get exhausted
     double NKbl_delta=(NKbl_d*NKbl_proliferation_rate_d+NK_Ab_d*NKa_d*NKa_expressing_receptor_d*m.Ab()-
             NKbl_apop_rate_d*NKbl_d-u_NK_TNF_d*NKbl_d*(m.TNF()/(m.TNF()+Ks_NK_m_TNF_d))/*-NKbl_d*NK_exh_d*/);

    D.push_back(NKbl_delta);
//    double  NKexh_delta=(NKa_d*NK_exh_d+NKbo_d*NK_exh_d+NKbo_Ab_d*NK_exh_d+NKbl_d*NK_exh_d-
//              NKexh_apop_rate_d*NKexh_d-u_NK_TNF_d*NKexh_d*(m.TNF()/(m.TNF()+Ks_NK_m_TNF_d)));
//     D.push_back(NKexh_delta);
     double NK_TymTr_incorporated_delta=0;
     if (m.TymidineTriteate()>0){
         NK_TymTr_incorporated_delta=((NK0_d*NK0_proliferation_rate_d+(NKa_d+NKbl_d)*NKa_proliferation_rate_d+
                                   (NKbo_d+NKbo_Ab_d)*NKbo_proliferation_rate_d)*m.Prol_TymTr());
     }
     D.push_back(NK_TymTr_incorporated_delta);

     return D;
 }



 std::vector<double> NK_cells::getState()const
  {
     std::vector<double> S;
     /// we update each subpopulation of cells independently and we take into account the transition from one state to the other

     /// the number of free cells (no Ag) they proliferate according and some of them are "lost" since they activate

     S.push_back(NK0_d);
     /** activated NK cells dynamics*/
     S.push_back(NKa_d);
 /// El porcentaje de células productoras de IL-12 está dentro de la constante

     /// the cells that have interacted with LT grow accordingly with the number of cells that have internalized the Ag and the
     /// number of LT cells, monocytes and NK expressing the receptor with the same affinity (Monocytes, NK or LT can
     /// interact only with one cell). We are supposing that all activated cells express receptor and ligand.

      S.push_back(NKbo_d);

     /// NK cells can interact with other cells or with the blocking mAb

     S.push_back(NKbo_Ab_d);
     /// the cells that have interacted with LT get exhausted

    S.push_back(NKbl_d);
//     S.push_back(NKexh_d);
     S.push_back(NK_TymTr_incorporated_d);

     return S;

 }




 void NK_cells::setState(std::vector<double> y)
 {
     /// we update each subpopulation of cells independently and we take into account the transition from one state to the other

     /// the number of free cells (no Ag) they proliferate according and some of them are "lost" since they activate

     NK0_d=y[0];
     /** activated NK cells dynamics*/
     NKa_d=y[1];
 /// El porcentaje de células productoras de IL-12 está dentro de la constante

     /// the cells that have interacted with LT grow accordingly with the number of cells that have internalized the Ag and the
     /// number of LT cells, monocytes and NK expressing the receptor with the same affinity (Monocytes, NK or LT can
     /// interact only with one cell). We are supposing that all activated cells express receptor and ligand.
      NKbo_d=y[2];

     /// NK cells can interact with other cells or with the blocking mAb

     NKbo_Ab_d=y[3];
     /// the cells that have interacted with LT get exhausted
     NKbl_d=y[4];
//     NKexh_d=y[5];
     NK_TymTr_incorporated_d=y[6];
  }
