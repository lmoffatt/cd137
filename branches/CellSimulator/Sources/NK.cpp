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
                   /// 5)Percentages of TNF productions of each type of NK
                   /*9*/ double percentage_TNF_NK0_prod_rate_,
                   /*10*/ double percentage_TNF_NKa_prod_rate_,
                   //*11*/ double percentage_TNF_NKbo_prod_rate_,
                   /// 6) Proliferation rates
                   /*12*/ double NK0_proliferation_rate_,
                   /*13*/ double NKa_proliferation_rate_,
                   //*14*/ double NKbo_proliferation_rate_,
                   /// 7) Apoptosis rates
                   /*15*/ double NK0_apop_rate_,
                   /*16*/ double NKa_apop_rate_,
                   //*17*/ double NKbo_apop_rate_,
                   /// 8) constant saturation of TNF for apoptosis
                   /*18*/ double Ks_NK_m_TNF_,
                   /// 9) conversion rates
                   /*19*/ double KaNK_,
                   /*20*/ double NK_NK_,
                   /*21*/ double NK_Ab_,
                   /// 10)Saturation constant of APC interaction for activation
                   /*22*/ double KsAPC_NK_,
                   /// 11)Saturation constant of NK_LT interaction
                   /*23*/ double Ksi_,
                   /*24*/ double Kst_,
                   /// 12) Percentages of cell expressing receptor
                   /*25*/ double NK0_expressing_receptor_,
                   /*26*/ double NKa_expressing_receptor_,
                   /// 13) Apoptosis rate for TNF
                   /*27*/ double u_NK_TNF_
                   ):

    /*1*/ NK0_d(init_NK_),
    NKa_d(0),
    NKbo_d(0),
    NKbo_Ab_d(0),
    NKbl_d (0),
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
    /// 5)Percentages of TNF productions of each type of NK
    /*9*/ percentage_TNF_NK0_prod_rate_d (percentage_TNF_NK0_prod_rate_),
    /*10*/ percentage_TNF_NKa_prod_rate_d (percentage_TNF_NKa_prod_rate_),
    //*11*/ percentage_TNF_NKbo_prod_rate_d (percentage_TNF_NKbo_prod_rate_),
    /// 6) Proliferation rates
    /*12*/ NK0_proliferation_rate_d (NK0_proliferation_rate_),
    /*13*/ NKa_proliferation_rate_d (NKa_proliferation_rate_),
    //*14*/ NKbo_proliferation_rate_d (NKbo_proliferation_rate_),
    /// 7) Apoptosis rates
    /*15*/ NK0_apop_rate_d (NK0_apop_rate_),
    /*16*/ NKa_apop_rate_d (NKa_apop_rate_),
    //*17*/ NKbo_apop_rate_d (NKbo_apop_rate_),
    /// 8) constant saturation of TNF for apoptosis
    /*18*/ Ks_NK_m_TNF_d (Ks_NK_m_TNF_),
    /// 9) conversion rates
    /*19*/ KaNK_d (KaNK_),
    /*20*/ NK_NK_d (NK_NK_),
    /*21*/ NK_Ab_d (NK_Ab_),
    /// 10)Saturation constant of NK interaction for activation
    /*22*/ KsAPC_NK_d (KsAPC_NK_),
    /// 11)Saturation constant of NK_LT interaction
    /*23*/ Ksi_d (Ksi_),
    /*24*/ Kst_d (Kst_),
    /// 12) Percentages of cell expressing receptor
    /*25*/ NK0_expressing_receptor_d (NK0_expressing_receptor_),
    /*26*/ NKa_expressing_receptor_d (NKa_expressing_receptor_),
    /// 13) Apoptosis rate for TNF
    /*27*/ u_NK_TNF_d (u_NK_TNF_)

{}

NK_cells::NK_cells(){}

NK_cells::NK_cells(const NK_cells& other):

     /*1*/ NK0_d(other.NK0_d),
     /*2*/ NKa_d(other.NKa_d),
     /*3*/ NKbo_d(other.NKbo_d),
     /*4*/ NKbo_Ab_d(other.NKbo_Ab_d),
     /*5*/ NKbl_d (other.NKbl_d),
     /*6*/ NK_TymTr_incorporated_d(other.NK_TymTr_incorporated_d),
     /*7*/ init_ratio_NK_d(other.init_ratio_NK_d),
     /*8*/ IFN_NK0_prod_rate_d(other.IFN_NK0_prod_rate_d),
     /*9*/ IFN_NKa_prod_rate_d(other.IFN_NKa_prod_rate_d),
     /*10*/ IFN_NKbo_prod_rate_d(other.IFN_NKbo_prod_rate_d),
     /*11*/ TNF_NK0_prod_rate_d(other.TNF_NK0_prod_rate_d),
     /*12*/ TNF_NKa_prod_rate_d(other.TNF_NKa_prod_rate_d),
     /*13*/ TNF_NKbo_prod_rate_d(other.TNF_NKbo_prod_rate_d),
     /*14*/ percentage_IFN_NK0_prod_rate_d(other.percentage_IFN_NK0_prod_rate_d),
     /*15*/ percentage_TNF_NK0_prod_rate_d(other.percentage_TNF_NK0_prod_rate_d),
     /*16*/ percentage_TNF_NKa_prod_rate_d(other.percentage_TNF_NKa_prod_rate_d),
     //*17*/ percentage_TNF_NKbo_prod_rate_d(other.percentage_TNF_NKbo_prod_rate_d),
     /*18*/ NK0_proliferation_rate_d (other.NK0_proliferation_rate_d),
     /*19*/ NKa_proliferation_rate_d(other.NKa_proliferation_rate_d),
     //*20*/ NKbo_proliferation_rate_d(other.NKbo_proliferation_rate_d),
     /*21*/ NK0_apop_rate_d(other.NK0_apop_rate_d),
     /*22*/ NKa_apop_rate_d(other.NKa_apop_rate_d),
     //*23*/ NKbo_apop_rate_d(other.NKbo_apop_rate_d),
     /*24*/ Ks_NK_m_TNF_d(other.Ks_NK_m_TNF_d),
     /*25*/ KaNK_d(other.KaNK_d),
     /*26*/ NK_NK_d(other.NK_NK_d),
     /*27*/ NK_Ab_d(other.NK_Ab_d),
     /*28*/ KsAPC_NK_d(other.KsAPC_NK_d),
     /*29*/ Ksi_d(other.Ksi_d),
     /*30*/ Kst_d(other.Kst_d),
     /*31*/ NK0_expressing_receptor_d(other.NK0_expressing_receptor_d),
     /*32*/ NKa_expressing_receptor_d(other.NKa_expressing_receptor_d),
     /*33*/ u_NK_TNF_d(other.u_NK_TNF_d)
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

    /*1*/ std::swap(one.NK0_d,other.NK0_d);
    /*2*/ std::swap(one.NKa_d,other.NKa_d);
    /*3*/ std::swap(one.NKbo_d,other.NKbo_d);
    /*4*/ std::swap(one.NKbo_Ab_d,other.NKbo_Ab_d);
    /*5*/ std::swap(one.NKbl_d ,other.NKbl_d);
    /*6*/ std::swap (one.NK_TymTr_incorporated_d,other.NK_TymTr_incorporated_d);
    /*7*/ std::swap(one.init_ratio_NK_d,other.init_ratio_NK_d);
    /*8*/ std::swap(one.IFN_NK0_prod_rate_d,other.IFN_NK0_prod_rate_d);
    /*9*/ std::swap(one.IFN_NKa_prod_rate_d,other.IFN_NKa_prod_rate_d);
    /*10*/ std::swap(one.IFN_NKbo_prod_rate_d,other.IFN_NKbo_prod_rate_d);
    /*11*/ std::swap(one.TNF_NK0_prod_rate_d,other.TNF_NK0_prod_rate_d);
    /*12*/ std::swap(one.TNF_NKa_prod_rate_d,other.TNF_NKa_prod_rate_d);
    /*13*/ std::swap(one.TNF_NKbo_prod_rate_d,other.TNF_NKbo_prod_rate_d);
    /*14*/ std::swap(one.percentage_IFN_NK0_prod_rate_d,other.percentage_IFN_NK0_prod_rate_d);
    /*15*/ std::swap(one.percentage_TNF_NK0_prod_rate_d,other.percentage_TNF_NK0_prod_rate_d);
    /*16*/ std::swap(one.percentage_TNF_NKa_prod_rate_d,other.percentage_TNF_NKa_prod_rate_d);
    //*17*/ std::swap(one.percentage_TNF_NKbo_prod_rate_d,other.percentage_TNF_NKbo_prod_rate_d);
    /*18*/ std::swap(one.NK0_proliferation_rate_d ,other.NK0_proliferation_rate_d),
    /*19*/ std::swap(one.NKa_proliferation_rate_d,other.NKa_proliferation_rate_d);
    //*20*/ std::swap(one.NKbo_proliferation_rate_d,other.NKbo_proliferation_rate_d);
    /*21*/ std::swap(one.NK0_apop_rate_d,other.NK0_apop_rate_d);
    /*22*/ std::swap(one.NKa_apop_rate_d,other.NKa_apop_rate_d);
    //*23*/ std::swap(one.NKbo_apop_rate_d,other.NKbo_apop_rate_d);
    /*24*/ std::swap(one.Ks_NK_m_TNF_d,other.Ks_NK_m_TNF_d);
    /*25*/ std::swap(one.KaNK_d,other.KaNK_d);
    /*26*/ std::swap(one.NK_NK_d,other.NK_NK_d);
    /*27*/ std::swap(one.NK_Ab_d,other.NK_Ab_d);
    /*28*/ std::swap(one.KsAPC_NK_d,other.KsAPC_NK_d);
    /*29*/ std::swap(one.Ksi_d,other.Ksi_d);
    /*30*/ std::swap(one.Kst_d,other.Kst_d);
    /*31*/ std::swap(one.NK0_expressing_receptor_d,other.NK0_expressing_receptor_d);
    /*32*/ std::swap(one.NKa_expressing_receptor_d,other.NKa_expressing_receptor_d);
    /*33*/ std::swap(one.u_NK_TNF_d,other.u_NK_TNF_d);

}



//*********************** main step for the NK cells(NK cell dynamics)*************************************************/
void NK_cells::update(double& time_step,const Media& m, const APC_cells& APC,const LT_cells& LT)
{
    /// we update each subpopulation of cells independently and we take into account the transition from one state to the other

    /// naive cell, die, proliferate and some of them are "lost" since they activate
    double NK0_delta=(NK0_d*NK0_proliferation_rate_d*m.prol_ratio()-
                      NK0_d*NK0_apop_rate_d-
                      NK0_d*m.Ag()*KaNK_d*
                      (
                          (APC.APCa()+
                           APC.APCbo()+
                           (APC.APCbo_TNF_production_rate())*(APC.APCbl()+APC.APCbo_Ab())
                           )/
                          ((APC.APCa()+
                            APC.APCbo()+
                            (APC.APCbo_TNF_production_rate())*(APC.APCbl()+APC.APCbo_Ab())
                            )+
                           KsAPC_NK_d
                           )
                       )*
                      ((APC.APCa()+APC.APCbl()+APC.APCbo()+APC.APCbo_Ab())
                       ))
                      *time_step;
    NK0_d+=NK0_delta;

    /// activated NK "came from NK0" cell die, proliferate and some of them are "lost" since they bound ligand or mAb
    double NKa_delta=(NKa_proliferation_rate_d*NKa_d*m.prol_ratio()+
                      NK0_d*m.Ag()*KaNK_d*
                      (
                          (APC.APCa()+
                           APC.APCbo()+
                           (APC.APCbo_TNF_production_rate())*(APC.APCbl()+APC.APCbo_Ab())
                           )/
                          ((APC.APCa()+
                            APC.APCbo()+
                            (APC.APCbo_TNF_production_rate())*(APC.APCbl()+APC.APCbo_Ab())
                            )+
                           KsAPC_NK_d
                           )
                       )*
                      ((APC.APCa()+APC.APCbl()+APC.APCbo()+APC.APCbo_Ab())
                       )-
                      NKa_d*NK_NK_d*NKa_expressing_receptor_d*(2*NKa_d*NKa_expressing_receptor_d+NKbo_d+NKbl_d)-
                      NKa_d*APC.APC_NK()*NKa_expressing_receptor_d*(APC.APCa()+APC.APCbo()+APC.APCbl()+APC.APCbo_Ab())-
                      NKa_d*NK_Ab_d*NKa_expressing_receptor_d*m.Ab()-
                      NKa_d*NKa_apop_rate_d-
                      NKa_d*u_NK_TNF_d*(m.TNF()/(m.TNF()+Ks_NK_m_TNF_d))
                      )*
                      time_step;
    NKa_d+=NKa_delta;
    /// El porcentaje de células productoras de IL-12 está dentro de la constante

    /// signalized NK cell "came" from NKa, proliferate, die and some of them are "lost" since they bound mAb
    double NKbo_delta= (NKbo_d*NKa_proliferation_rate_d*m.prol_ratio()+
                        NK_NK_d*NKa_d*NKa_expressing_receptor_d*(2*NKa_d*NKa_expressing_receptor_d+NKbo_d+NKbl_d)+
                        APC.APC_NK()*NKa_d*NKa_expressing_receptor_d*(APC.APCa()+APC.APCbo()+APC.APCbl()+APC.APCbo_Ab())-
                        NKa_apop_rate_d*NKbo_d-
                        u_NK_TNF_d*NKbo_d*(m.TNF()/(m.TNF()+Ks_NK_m_TNF_d))-
                        NKbo_d*NK_Ab_d*m.Ab())*time_step;

    NKbo_d+=NKbo_delta;

    /// signalized NK cells are blocke if they bind with the blocking mAb, proliferate and die

    double NKbo_Ab_delta=(NKbo_Ab_d*NKa_proliferation_rate_d*m.prol_ratio()+
                          NKbo_d*NK_Ab_d*m.Ab()-
                          NKa_apop_rate_d*NKbo_Ab_d-
                          u_NK_TNF_d*NKbo_Ab_d*(m.TNF()/(m.TNF()+Ks_NK_m_TNF_d)))
            *time_step;
    NKbo_Ab_d+=NKbo_Ab_delta;

    /// block NK cells are activated cells that binds to mAb, proliferate and die
    double NKbl_delta=(NKbl_d*NKa_proliferation_rate_d*m.prol_ratio()+
                       NK_Ab_d*NKa_d*NKa_expressing_receptor_d*m.Ab()-
                       NKa_apop_rate_d*NKbl_d-
                       u_NK_TNF_d*NKbl_d*(m.TNF()/(m.TNF()+Ks_NK_m_TNF_d)))*time_step;

    NKbl_d+=NKbl_delta;

    double NK_TymTr_incorporated_delta;
    if (m.TymidineTriteate()>0){
        NK_TymTr_incorporated_delta=
                (((NK0_d*NK0_proliferation_rate_d*m.prol_ratio()+(NKa_d+NKbl_d)*NKa_proliferation_rate_d*m.prol_ratio()+
                  (NKbo_d+NKbo_Ab_d)*NKa_proliferation_rate_d*m.prol_ratio()))*m.Prol_TymTr()
                 )*time_step;
        NK_TymTr_incorporated_d+=NK_TymTr_incorporated_delta;
    }

}

double NK_cells::num_NK() const
{
    double sum=NK0_d+NKa_d+NKbo_d+NKbo_Ab_d+NKbl_d;
    return sum;
}

double NK_cells::NK_IFNgamma_production_rate() const
{
    double sum=IFN_NK0_prod_rate_d*NK0_d*percentage_IFN_NK0_prod_rate_d+
            IFN_NKa_prod_rate_d*(NKa_d+NKbl_d)+
            (NKbo_d+NKbo_Ab_d)*IFN_NKa_prod_rate_d/IFN_NKbo_prod_rate_d;
    return sum;
}

double NK_cells::NK_TNF_production_rate() const
{
    double sum=TNF_NK0_prod_rate_d*NK0_d*percentage_TNF_NK0_prod_rate_d+
            TNF_NKa_prod_rate_d*(NKa_d+NKbl_d)*percentage_TNF_NKa_prod_rate_d+
            (NKbo_d+NKbo_Ab_d)*percentage_TNF_NKa_prod_rate_d*TNF_NKa_prod_rate_d/TNF_NKbo_prod_rate_d;
    return sum;
}

double NK_cells::percentage_NK_expressing_receptor() const
{
    double sum=(NK0_expressing_receptor_d*NK0_d+NKa_expressing_receptor_d*NKa_d+NKbo_d+NKbo_Ab_d+NKbl_d)*100/num_NK();
    return sum;
}

double NK_cells::percentage_NK_producing_IFN() const
{
    double sum=100*((percentage_IFN_NK0_prod_rate_d*NK0_d)+NKa_d+NKbl_d+NKbo_d+NKbl_d)/num_NK();
    return sum;
}

double NK_cells::percentage_NK_producing_TNF() const
{
    double sum=100*(percentage_TNF_NK0_prod_rate_d*NK0_d
                    +percentage_TNF_NKa_prod_rate_d*(NKa_d+NKbl_d)
                    +percentage_TNF_NKa_prod_rate_d*(NKbo_d+NKbo_Ab_d))/num_NK();
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
    s<<"\n NK tymTR incorporate \t"<<c.NK_TymTr_incorporated_d;


    if (0)
    {

        /*1*/ s<<"\n init_ratio_NK_ \t"<<c.init_ratio_NK_d;
        /*2*/ s<<"\n IFN_NK0_prod_rate_ \t"<<c.IFN_NK0_prod_rate_d;
        /*3*/ s<<"\n IFN_NK0_prod_rate_ \t"<<c.IFN_NK0_prod_rate_d;
        /*4*/ s<<"\n IFN_NKbo_prod_rate_ \t"<<c.IFN_NKbo_prod_rate_d;
        /*5*/ s<<"\n TNF_NK0_prod_rate_ \t"<<c.TNF_NK0_prod_rate_d;
        /*6*/ s<<"\n TNF_NKa_prod_rate_ \t"<<c.TNF_NKa_prod_rate_d;
        /*7*/ s<<"\n TNF_NKbo_prod_rate_ \t"<<c.TNF_NKbo_prod_rate_d;
        /*8*/ s<<"\n percentage_IFN_NK0_prod_rate_ \t"<<c.percentage_IFN_NK0_prod_rate_d;
        /*9*/ s<<"\n percentage_TNF_NK0_prod_rate_ \t"<<c.percentage_TNF_NK0_prod_rate_d;
        /*10*/ s<<"\n percentage_TNF_NKa_prod_rate_ \t"<<c.percentage_TNF_NKa_prod_rate_d;
        //*11*/ s<<"\n percentage_TNF_NKbo_prod_rate_ \t"<<c.percentage_TNF_NKbo_prod_rate_d;
        /*12*/ s<<"\n NK0_proliferation_rate_d \t"<<c.NK0_proliferation_rate_d;
        /*13*/ s<<"\n NKa_proliferation_rate_ \t"<<c.NKa_proliferation_rate_d;
        //*14*/ s<<"\n NKbo_proliferation_rate_ \t"<<c.NKbo_proliferation_rate_d;
        /*15*/ s<<"\n NK0_apop_rate_ \t"<<c.NK0_apop_rate_d;
        //*16*/ s<<"\n NKa_apop_rate_ \t"<<c.NKa_apop_rate_d;
        //*17*/ s<<"\n NKbo_apop_rate_ \t"<<c.NKbo_apop_rate_d;
        /*18*/ s<<"\n Ks_NK_m_TNF_ \t"<<c.Ks_NK_m_TNF_d;
        /*19*/ s<<"\n KaNK_ \t"<<c.KaNK_d;
        /*20*/ s<<"\n NK_NK_ \t"<<c.NK_NK_d;
        /*21*/ s<<"\n NK_Ab_ \t"<<c.NK_Ab_d;
        /*22*/ s<<"\n KsAPC_NK_ \t"<<c.KsAPC_NK_d;
        /*23*/ s<<"\n Ksi_ \t"<<c.Ksi_d;
        /*24*/ s<<"\n Kst_ \t"<<c.Kst_d;
        /*25*/ s<<"\n NK0_expressing_receptor_ \t"<<c.NK0_expressing_receptor_d;
        /*26*/ s<<"\n NKa_expressing_receptor_ \t"<<c.NKa_expressing_receptor_d;
        /*27*/ s<<"\n u_NK_TNF_ \t"<<c.u_NK_TNF_d;
        }
    return s;
}

NK_cells::NK_cells(const Parameters& p, const Treatment& t):
    /// Variables
   NK0_d(
        (1.0-p.mean_ratio("init_K_ratio_LT"))*
        (p.mean_ratio("init_K_ratio_NK_APC"))
        *t.init_cells
        ),
    NKa_d(0.0),
    NKbo_d(0.0),
    NKbo_Ab_d(0.0),
    NKbl_d (0.0),
    NK_TymTr_incorporated_d(0.0),
    init_ratio_NK_d(p.mean_ratio("init_K_ratio_NK_APC")),
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

    /// 5)Percentages of TNF productions of each type of NK
    /*9*/ percentage_TNF_NK0_prod_rate_d (p.mean_ratio("Kpercentage_TNF_NK0_prod_rate")),
    /*10*/ percentage_TNF_NKa_prod_rate_d (p.mean_ratio("Kpercentage_TNF_NKa_prod_rate")),
    //*11*/ percentage_TNF_NKbo_prod_rate_d (p.mean_ratio("Kpercentage_TNF_NKbo_prod_rate")),

    /// 6) Proliferation rates
    /*12*/ NK0_proliferation_rate_d (p.mean("NK0_proliferation_rate")),
    /*13*/ NKa_proliferation_rate_d (p.mean("NKa_proliferation_rate")),
    //*14*/ NKbo_proliferation_rate_d (p.mean("NKbo_proliferation_rate")),

    /// 7) Apoptosis rates
    /*15*/ NK0_apop_rate_d (p.mean("NK0_apop_rate")),
    /*16*/ NKa_apop_rate_d (p.mean("NKa_apop_rate")),
    //*17*/ NKbo_apop_rate_d (p.mean("NKbo_apop_rate")),

    /// 8) constant saturation of TNF for apoptosis
    /*18*/ Ks_NK_m_TNF_d (p.mean("Ks_NK_m_TNF")),

    /// 9) conversion rates
    /*19*/ KaNK_d (p.mean("KaNK")),
    /*20*/ NK_NK_d (p.mean("NK_NK")),
    /*21*/ NK_Ab_d (p.mean("NK_Ab")),

    /// 10)Saturation constant of NK interaction for activation
    /*22*/ KsAPC_NK_d (p.mean("KsAPC_NK")),

    /// 11)Saturation constant of NK_LT interaction
    /*23*/ Ksi_d (p.mean("NK_Ksi")),
    /*24*/ Kst_d (p.mean("NK_Kst")),

    /// 12) Percentages of cell expressing receptor
    /*25*/ NK0_expressing_receptor_d (p.mean_ratio("NK0_Kratio_expressing_receptor")),
    /*26*/ NKa_expressing_receptor_d (p.mean_ratio("NKa_Kratio_expressing_receptor")),
    /// 13) Apoptosis rate for TNF
    /*27*/ u_NK_TNF_d (p.mean("u_NK_TNF"))
{

}




std::vector<double> NK_cells::Derivative(const Media& m, const APC_cells& APC)const
{

    std::vector<double> D;
    /// we update each subpopulation of cells independently and we take into account the transition from one state to the other

    /// naive cell, die, proliferate and some of them are "lost" since they activate

    double NK0_delta=(NK0_d*NK0_proliferation_rate_d*m.prol_ratio()-
                      NK0_d*NK0_apop_rate_d-
                      NK0_d*m.Ag()*KaNK_d*
                      ((APC.APCa()+
                       APC.APCbl()+
                      (APC.APCbo_TNF_production_rate()/APC.APCa_TNF_production_rate())*(APC.APCbl()+APC.APCbo_Ab()))/
                      ((APC.APCa()+
                       APC.APCbl()+
                      (APC.APCbo_TNF_production_rate()/APC.APCa_TNF_production_rate())*(APC.APCbl()+APC.APCbo_Ab()))+ KsAPC_NK_d))*
                      ((APC.APCa()+APC.APCbl()+APC.APCbl()+APC.APCbo_Ab())/
                      ((APC.APCa()+APC.APCbl()+APC.APCbl()+APC.APCbo_Ab())+ KsAPC_NK_d)));
    D.push_back(NK0_delta);

    /// activated NK "came from NK0" cell die, proliferate and some of them are "lost" since they bound ligand or mAb
    double NKa_delta=(NK0_d*m.Ag()*KaNK_d*
                      ((APC.APCa()+
                       APC.APCbl()+
                      (APC.APCbo_TNF_production_rate()/APC.APCa_TNF_production_rate())*(APC.APCbl()+APC.APCbo_Ab()))/
                      ((APC.APCa()+
                       APC.APCbl()+
                      (APC.APCbo_TNF_production_rate()/APC.APCa_TNF_production_rate())*(APC.APCbl()+APC.APCbo_Ab()))+ KsAPC_NK_d))*
                      ((APC.APCa()+APC.APCbl()+APC.APCbl()+APC.APCbo_Ab())/
                      ((APC.APCa()+APC.APCbl()+APC.APCbl()+APC.APCbo_Ab())+ KsAPC_NK_d))-
                      NK_NK_d*NKa_d*NKa_expressing_receptor_d*(2*NKa_d*NKa_expressing_receptor_d+NKbo_d+NKbl_d)-
                      APC.APC_NK()*NKa_d*NKa_expressing_receptor_d*(APC.APCa()+APC.APCbo()+APC.APCbl())-
                      NK_Ab_d*NKa_d*NKa_expressing_receptor_d*m.Ab()-
                      NKa_apop_rate_d*NKa_d-
                      u_NK_TNF_d*NKa_d*(m.TNF()/(m.TNF()+Ks_NK_m_TNF_d)));
    D.push_back(NKa_delta);
    /// El porcentaje de células productoras de IL-12 está dentro de la constante

    /// signalized NK cell "came" from NKa, proliferate, die and some of them are "lost" since they bound mAb
    double NKbo_delta= (NKbo_d*NKa_proliferation_rate_d*m.prol_ratio()+
                        NK_NK_d*NKa_d*NKa_expressing_receptor_d*NKa_d*NKa_expressing_receptor_d+
                        NK_NK_d*NKa_d*NKa_expressing_receptor_d*NKbo_d+
                        APC.APC_NK()*NKa_d*NKa_expressing_receptor_d*APC.APCa()+
                        APC.APC_NK()*NKa_d*NKa_expressing_receptor_d*APC.APCbo()-
                        NKa_apop_rate_d*NKbo_d-u_NK_TNF_d*NKbo_d*(m.TNF()/(m.TNF()+Ks_NK_m_TNF_d))-
                        NKbo_d*NK_Ab_d*m.Ab()/*-NKbo_d*NK_exh_d*/);

    D.push_back(NKbo_delta);

    /// signalized NK cells are blocke if they bind with the blocking mAb, proliferate and die

    double NKbo_Ab_delta=(NKbo_Ab_d*NKa_proliferation_rate_d*m.prol_ratio()+
                          NKbo_d*NK_Ab_d*m.Ab()-
                          NKa_apop_rate_d*NKbo_Ab_d-
                          u_NK_TNF_d*NKbo_Ab_d*(m.TNF()/(m.TNF()+Ks_NK_m_TNF_d)));
    D.push_back(NKbo_Ab_delta);
    /// block NK cells are activated cells that binds to mAb, proliferate and die
    double NKbl_delta=(NKbo_Ab_d*NKa_proliferation_rate_d*m.prol_ratio()+
                       NKbo_d*NK_Ab_d*m.Ab()-
                       NKa_apop_rate_d*NKbo_Ab_d-
                       u_NK_TNF_d*NKbo_Ab_d*(m.TNF()/(m.TNF()+Ks_NK_m_TNF_d)));

    D.push_back(NKbl_delta);

    double NK_TymTr_incorporated_delta=0;
    if (m.TymidineTriteate()>0){
        NK_TymTr_incorporated_delta=((NK0_d*NK0_proliferation_rate_d*m.prol_ratio()+
                                      (NKa_d+NKbl_d)*NKa_proliferation_rate_d*m.prol_ratio()+
                                      (NKbo_d+NKbo_Ab_d)*NKa_proliferation_rate_d*m.prol_ratio())*m.Prol_TymTr());
    }
    D.push_back(NK_TymTr_incorporated_delta);

    return D;
}



std::vector<double> NK_cells::getState()const
{
    std::vector<double> S;
    S.push_back(NK0_d);
    S.push_back(NKa_d);
    S.push_back(NKbo_d);
    S.push_back(NKbo_Ab_d);
    S.push_back(NKbl_d);
    S.push_back(NK_TymTr_incorporated_d);

    return S;

}




void NK_cells::setState(std::vector<double> y)
{
    NK0_d=y[0];
    NKa_d=y[1];
    NKbo_d=y[2];
    NKbo_Ab_d=y[3];
    NKbl_d=y[4];
    NK_TymTr_incorporated_d=y[5];
}
