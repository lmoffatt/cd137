#include "Includes/Media.h"
#include "Includes/APC.h"
#include "Includes/LT.h"
#include "Includes/NK.h"
#include <cmath>


APC_cells::APC_cells(/// 1) Init number of APC
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
                     /// 5)Percentages of TNF productions of APC0
                     /*11*/ double percentage_TNF_APC0_prod_rate_,
                     /// 6) Proliferation rates
                     /*12*/ double APC_bound_proliferation_rate_,
                     /// 7) Apoptosis rates
                     /*13*/ double APC0_apop_rate_,
                     /*14*/ double APCa_apop_rate_,
                     /*15*/ double APCbo_apop_rate_,
                     /// 8) constant saturation of TNF for apoptosis
                     /*16*/ double Ks_APC_m_TNF_,
                     /// 9) conversion rates
                     /*17*/ double APC_Ag_,
                     /*18*/ double APC_APC_,
                     /*19*/ double APC_NK_,
                     /*20*/ double APC_LT_1_,
                     /*21*/ double APC_LT_2_,
                     /*22*/ double APC_Ab_,
                     /// 10)Saturation constant of IFN and TNF for activation
                     /*23*/ double KsAPC_LT_,
                     /// 11)Saturation constant of APC_LT interaction
                     /*24*/ double Ksi_,
                     /*25*/ double Kst_,
                     /// 12) Percentages of cell expressing receptor
                     /*26*/ double APC0_expressing_receptor_,
                     /// 13) Apoptosis rate for TNF
                     /*27*/ double u_APC_TNF_


                     ):

    /*1*/ APC0_d(init_APC_),
    APCa_d(0),
    APCbo_d(0),
    APCbo_Ab_d(0),
    APCbl_d (0),
    APC_TymTr_incorporated_d(0),
    /*2*/ IFN_APC0_prod_rate_d(IFN_APC0_prod_rate_),
    /*3*/ IFN_APCa_prod_rate_d(IFN_APCa_prod_rate_),
    /*4*/ IFN_APCbo_prod_rate_d(IFN_APCbo_prod_rate_),
    /*5*/ TNF_APC0_prod_rate_d(TNF_APC0_prod_rate_),
    /*6*/ TNF_APCa_prod_rate_d(TNF_APCa_prod_rate_),
    /*7*/ TNF_APCbo_prod_rate_d(TNF_APCbo_prod_rate_),
    /*8*/ percentage_IFN_APC0_prod_rate_d (percentage_IFN_APC0_prod_rate_),
    /*9*/ percentage_IFN_APCa_prod_rate_d (percentage_IFN_APCa_prod_rate_),
    /*10*/ percentage_IFN_APCbo_prod_rate_d (percentage_IFN_APCbo_prod_rate_),
    /*11*/ percentage_TNF_APC0_prod_rate_d (percentage_TNF_APC0_prod_rate_),
    /*12*/ APC_bound_proliferation_rate_d (APC_bound_proliferation_rate_),
    /*13*/ APC0_apop_rate_d(APC0_apop_rate_),
    /*14*/ APCa_apop_rate_d (APCa_apop_rate_),
    /*15*/ APCbo_apop_rate_d (APCbo_apop_rate_),
    /*16*/ Ks_APC_m_TNF_d (Ks_APC_m_TNF_),
    /*17*/ APC_Ag_d (APC_Ag_),
    /*18*/ APC_APC_d(APC_APC_),
    /*19*/ APC_NK_d(APC_NK_),
    /*20*/ APC_LT_1_d(APC_LT_1_),
    /*21*/ APC_LT_2_d(APC_LT_2_),
    /*22*/ APC_Ab_d(APC_Ab_),
    /*23*/ KsAPC_LT_d(KsAPC_LT_),
    /*24*/ Ksi_d(Ksi_),
    /*25*/ Kst_d(Kst_),
    /*26*/ APC0_expressing_receptor_d(APC0_expressing_receptor_),
    /*27*/ u_APC_TNF_d (u_APC_TNF_)


{}


APC_cells::APC_cells(){}

APC_cells::APC_cells(const APC_cells& other):
    /*1*/ APC0_d(other.APC0_d),
    /*2*/ APCa_d(other.APCa_d),
    /*3*/ APCbo_d(other.APCbo_d),
    /*4*/ APCbo_Ab_d(other.APCbo_Ab_d),
    /*5*/ APCbl_d (other.APCbl_d),
    /*6*/ APC_TymTr_incorporated_d(other.APC_TymTr_incorporated_d),
    /*7*/ init_ratio_APC_d (other.init_ratio_APC_d),
    /*8*/ IFN_APC0_prod_rate_d(other.IFN_APC0_prod_rate_d),
    /*9*/ IFN_APCa_prod_rate_d(other.IFN_APCa_prod_rate_d),
    /*10*/ IFN_APCbo_prod_rate_d(other.IFN_APCbo_prod_rate_d),
    /*11*/ TNF_APC0_prod_rate_d(other.TNF_APC0_prod_rate_d),
    /*12*/ TNF_APCa_prod_rate_d(other.TNF_APCa_prod_rate_d),
    /*13*/ TNF_APCbo_prod_rate_d(other.TNF_APCbo_prod_rate_d),
    /*14*/ percentage_IFN_APC0_prod_rate_d (other.percentage_IFN_APC0_prod_rate_d),
    /*15*/ percentage_IFN_APCa_prod_rate_d (other.percentage_IFN_APCa_prod_rate_d),
    /*16*/ percentage_IFN_APCbo_prod_rate_d (other.percentage_IFN_APCbo_prod_rate_d),
    /*17*/ percentage_TNF_APC0_prod_rate_d (other.percentage_TNF_APC0_prod_rate_d),
    /*18*/ APC_bound_proliferation_rate_d (other.APC_bound_proliferation_rate_d),
    /*19*/ APC0_apop_rate_d(other.APC0_apop_rate_d),
    /*20*/ APCa_apop_rate_d (other.APCa_apop_rate_d),
    /*21*/ APCbo_apop_rate_d (other.APCbo_apop_rate_d),
    /*22*/ Ks_APC_m_TNF_d (other.Ks_APC_m_TNF_d),
    /*23*/ APC_Ag_d (other.APC_Ag_d),
    /*24*/ APC_APC_d(other.APC_APC_d),
    /*25*/ APC_NK_d(other.APC_NK_d),
    /*26*/ APC_LT_1_d(other.APC_LT_1_d),
    /*27*/ APC_LT_2_d(other.APC_LT_2_d),
    /*28*/ APC_Ab_d(other.APC_Ab_d),
    /*29*/ KsAPC_LT_d(other.KsAPC_LT_d),
    /*30*/ Ksi_d(other.Ksi_d),
    /*31*/ Kst_d(other.Kst_d),
    /*32*/ APC0_expressing_receptor_d(other.APC0_expressing_receptor_d),
    /*33*/ u_APC_TNF_d (other.u_APC_TNF_d)

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
    /*1*/ std::swap(one.APC0_d,other.APC0_d);
    /*2*/ std::swap(one.APCa_d,other.APCa_d);
    /*3*/ std::swap(one.APCbo_d,other.APCbo_d);
    /*4*/ std::swap(one.APCbo_Ab_d,other.APCbo_Ab_d);
    /*5*/ std::swap(one.APCbl_d ,other.APCbl_d);
    /*6*/ std::swap (one.APC_TymTr_incorporated_d,other.APC_TymTr_incorporated_d);
    /*7*/ std::swap (one.init_ratio_APC_d,other.init_ratio_APC_d);
    /*8*/ std::swap(one.IFN_APC0_prod_rate_d,other.IFN_APC0_prod_rate_d);
    /*9*/ std::swap(one.IFN_APCa_prod_rate_d,other.IFN_APCa_prod_rate_d);
    /*10*/ std::swap(one.IFN_APCbo_prod_rate_d,other.IFN_APCbo_prod_rate_d);
    /*11*/ std::swap(one.TNF_APC0_prod_rate_d,other.TNF_APC0_prod_rate_d);
    /*12*/ std::swap(one.TNF_APCa_prod_rate_d,other.TNF_APCa_prod_rate_d);
    /*13*/ std::swap(one.TNF_APCbo_prod_rate_d,other.TNF_APCbo_prod_rate_d);
    /*14*/ std::swap(one.percentage_IFN_APC0_prod_rate_d ,other.percentage_IFN_APC0_prod_rate_d);
    /*15*/ std::swap(one.percentage_IFN_APCa_prod_rate_d ,other.percentage_IFN_APCa_prod_rate_d);
    /*16*/ std::swap(one.percentage_IFN_APCbo_prod_rate_d ,other.percentage_IFN_APCbo_prod_rate_d);
    /*17*/ std::swap(one.percentage_TNF_APC0_prod_rate_d ,other.percentage_TNF_APC0_prod_rate_d);
    /*18*/ std::swap(one.APC_bound_proliferation_rate_d ,other.APC_bound_proliferation_rate_d);
    /*19*/ std::swap(one.APC0_apop_rate_d,other.APC0_apop_rate_d);
    /*20*/ std::swap(one.APCa_apop_rate_d ,other.APCa_apop_rate_d);
    /*21*/ std::swap(one.APCbo_apop_rate_d ,other.APCbo_apop_rate_d);
    /*22*/ std::swap(one.Ks_APC_m_TNF_d ,other.Ks_APC_m_TNF_d);
    /*23*/ std::swap(one.APC_Ag_d ,other.APC_Ag_d);
    /*24*/ std::swap(one.APC_APC_d,other.APC_APC_d);
    /*25*/ std::swap(one.APC_NK_d,other.APC_NK_d);
    /*26*/ std::swap(one.APC_LT_1_d,other.APC_LT_1_d);
    /*27*/ std::swap(one.APC_LT_2_d,other.APC_LT_2_d);
    /*28*/ std::swap(one.APC_Ab_d,other.APC_Ab_d);
    /*29*/ std::swap(one.KsAPC_LT_d,other.KsAPC_LT_d);
    /*30*/ std::swap(one.Ksi_d,other.Ksi_d);
    /*31*/ std::swap(one.Kst_d,other.Kst_d);
    /*32*/ std::swap(one.APC0_expressing_receptor_d,other.APC0_expressing_receptor_d);
    /*33*/ std::swap(one.u_APC_TNF_d ,other.u_APC_TNF_d);
}


/// **************** Main step for APC (APC dynamics)**************************************
void APC_cells::update(double& time_step,const Media& m, const NK_cells& NK, const LT_cells& LT)
{
    /// we update each subpopulation of cells independently and we take into account the transition from one state to the other

    /// the number of naive cells (no Ag) die and some of them are "lost" since they internalize the Ag

    double APC0_delta=(-APC0_apop_rate_d*APC0_d //apoptosis
                       // activation in inflamatory context APC0_d ->APCa_d
                       -m.Ag()*APC0_d*APC_Ag_d*(m.IFNgamma()/(m.IFNgamma()+ Ksi_d))*(m.TNF()/(m.TNF()+Kst_d))
                       //activation basal (dendritic cells)
                       -APC_Ag_d*m.Ag()*APC0_d)*time_step;

    APC0_d+=APC0_delta;
    /// activated APC cells die and some of them are "lost" since they bind receptor or mAb, they came from APC0
    double APCa_delta_pos=(m.Ag()*APC0_d*APC_Ag_d*(m.IFNgamma()/(m.IFNgamma()+ Ksi_d))*(m.TNF()/(m.TNF()+Kst_d))+
                           APC_Ag_d*m.Ag()*APC0_d)*time_step;
    double APCa_delta_neg= (-APCa_apop_rate_d*APCa_d
                            -u_APC_TNF_d*APCa_d*(m.TNF()/(m.TNF()+ Ks_APC_m_TNF_d))
                            -APC_APC_d*APCa_d*(2*APCa_d+APCbo_d)
                            -APC_NK_d*APCa_d*(NK.NKa()*NK.NKa_expressing_receptor()+NK.NKbo())
                            -APC_LT_1_d*LT.LT_Ab()/(LT.LT_Ab()+m.Ab())*LT.LT0()*(APCa_d)
                            -APC_Ab_d*m.Ab()*APCa_d
                            )*time_step;

    APCa_d+=APCa_delta_pos+APCa_delta_neg;
    /// only bound APC could proliferate. They also die. They came from APCa. Some of them are "lost" since they bind mAb
    double APCbo_delta=(
                APC_APC_d*APCa_d*(2*APCa_d+APCbo_d)
                +APC_NK_d*APCa_d*(NK.NKa()*NK.NKa_expressing_receptor()+NK.NKbo())
                +APC_LT_1_d*LT.LT_Ab()/(LT.LT_Ab()+m.Ab())*LT.LT0()*(APCa_d)
//                +APCbo_d*APC_bound_proliferation_rate_d*m.prol_ratio()
                -APCbo_d*APCbo_apop_rate_d
                -APCbo_d*APC_Ab_d*m.Ab()
                -u_APC_TNF_d*APCbo_d*(m.TNF()/(m.TNF()+ Ks_APC_m_TNF_d))
                )*time_step;

    APCbo_d+=APCbo_delta;

    /// block APC has ligand free, so they can be signalized and "lost". They came from APCa.
    double APCbl_delta= (APC_Ab_d*m.Ab()*APCa_d
                         -APCa_apop_rate_d*APCbl_d
                         -APCbl_d*APC_APC_d*(APCa_d+APCbo_d)
                         -APCbl_d*APC_NK_d*(NK.NKa()+NK.NKbo())
                         -u_APC_TNF_d*APCbl_d*(m.TNF()/(m.TNF()+ Ks_APC_m_TNF_d)))*time_step;
    APCbl_d+=APCbl_delta;

    /// These cells are signalized and blocked. They came from APCbl or APCbo. They can proliferate or die
    double APCbo_Ab_delta= (APCbo_d*APC_Ab_d*m.Ab()
                            +APCbl_d*APC_APC_d*(APCa_d+APCbo_d)
                            +APCbl_d*APC_APC_d*(NK.NKa()+NK.NKbo())
//                            +APCbo_Ab_d*APC_bound_proliferation_rate_d*m.prol_ratio()
                            -APCbo_apop_rate_d*APCbo_Ab_d
                            -APCbo_Ab_d*u_APC_TNF_d*APCbo_Ab_d*(m.TNF()/(m.TNF()+ Ks_APC_m_TNF_d)))*time_step;
    APCbo_Ab_d+=APCbo_Ab_delta;
    double APC_TymTr_incorporated_delta;
    if (m.TymidineTriteate()>0)
    {
        APC_TymTr_incorporated_delta=(
                    APCbo_d*APC_bound_proliferation_rate_d*m.prol_ratio()*m.Prol_TymTr()
                    )*time_step;
        APC_TymTr_incorporated_d+=APC_TymTr_incorporated_delta;
    }

}

/// 1) Total number of cells
double APC_cells::num_APC() const
{
    double sum=APC0_d+APCa_d+APCbo_d+APCbl_d+APCbo_Ab_d;
    return sum;
}


/// 2) number of cells (5)
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
    return APCbo_Ab_d;
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

/// 3) Percentage of cells expressing receptor
double APC_cells::percentage_cell_expressing_receptor() const
{
    double sum=100*(APC0_d*APC0_expressing_receptor_d+APCa_d+APCbo_d+APCbl_d+APCbo_Ab_d)/num_APC();
    return sum;
}

/// 4) Cytokines production rate and producing cells (6)

double APC_cells::APC_IFNgamma_production_rate() const
{
    double sum=APC0_d*IFN_APC0_prod_rate_d*percentage_IFN_APC0_prod_rate_d+
            (APCa_d+APCbl_d)*IFN_APCa_prod_rate_d*percentage_IFN_APCa_prod_rate_d+
            (APCbo_d+APCbo_Ab_d)*IFN_APCbo_prod_rate_d*percentage_IFN_APCbo_prod_rate_d;
    return sum;
}

double APC_cells::APC_TNF_production_rate() const
{
    double sum=APC0_d*TNF_APC0_prod_rate_d*percentage_TNF_APC0_prod_rate_d+
            (APCa_d+APCbl_d)*TNF_APCa_prod_rate_d+
            (APCbo_d+APCbo_Ab_d)*TNF_APCbo_prod_rate_d;
    return sum;
}

double APC_cells::percentage_APC_producing_IFN() const
{
    double sum=100*(percentage_IFN_APC0_prod_rate_d*APC0_d+
                    percentage_IFN_APCa_prod_rate_d*(APCa_d+APCbl_d)+
                    percentage_IFN_APCbo_prod_rate_d*(APCbo_d+APCbo_Ab_d))/num_APC();
    return sum;
}


double APC_cells::percentage_APC_producing_TNF() const
{
    double sum=100*(percentage_TNF_APC0_prod_rate_d*APC0_d+APCa_d+APCbl_d+APCbo_d+APCbo_Ab_d)/num_APC();
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
    return TNF_APCbo_prod_rate_d;
}


/// 5) Union rates of APC
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
    s<<"\n APC tymTR incorporate \t"<<c.APC_TymTr_incorporated_d;

    if (0)
    {
        /*1*/ s<<"\n///TNF and INF Poductions rates of each type of APC\n";
        /*2*/ s<<"\n IFN_APC0_prod_rate \t"<<c.IFN_APC0_prod_rate_d;
        /*3*/ s<<"\n IFN_APCa_prod_rate \t"<<c.IFN_APCa_prod_rate_d;
        /*4*/ s<<"\n IFN_APCbo_prod_rate \t"<<c.IFN_APCbo_prod_rate_d;
        /*5*/ s<<"\n TNF_APC0_prod_rate \t"<<c.TNF_APC0_prod_rate_d;
        /*6*/ s<<"\n TNF_APCa_prod_rate \t"<<c.TNF_APCa_prod_rate_d;
        /*7*/ s<<"\n TNF_APCbo_prod_rate \t"<<c.TNF_APCbo_prod_rate_d;
              s<<"\n/// those are parameters that do not vary\n";
        /*8*/ s<<"\n percentage_IFN_APC0_prod_rate_d \t"<<c.percentage_IFN_APC0_prod_rate_d;
        /*9*/ s<<"\n percentage_IFN_APCa_prod_rate_d \t"<<c.percentage_IFN_APCa_prod_rate_d;
        /*10*/ s<<"\n percentage_IFN_APCbo_prod_rate_d \t"<<c.percentage_IFN_APCbo_prod_rate_d;
        /*11*/ s<<"\n percentage_TNF_APC0_prod_rate_d \t"<<c.percentage_TNF_APC0_prod_rate_d;
        /*12*/ s<<"\n APC0_apop_rate \t"<<c.APC0_apop_rate_d;
        /*13*/ s<<"\n APCa_apop_rate \t"<<c.APCa_apop_rate_d;
        /*14*/ s<<"\n APCbo_apop_rate \t"<<c.APCbo_apop_rate_d;
        /*15*/ s<<"\n APC_bound_proliferation_rate \t"<<c.APC_bound_proliferation_rate_d;
        /*16*/ s<<"\n APC_Ag \t"<<c.APC_Ag_d;
        /*17*/ s<<"\n APC_APC \t"<<c.APC_APC_d;
        /*18*/ s<<"\n APC_NK_d \t"<<c.APC_NK_d;
        /*19*/ s<<"\n APC_LT_1 \t"<<c.APC_LT_1_d;
        /*20*/ s<<"\n APC_LT_1 \t"<<c.APC_LT_2_d;
        /*21*/ s<<"\n APC_Ab \t"<<c.APC_Ab_d;
        /*22*/ s<<"\n Ks_APC_m_TNF_d \t"<<c.Ks_APC_m_TNF_d;
        /*23*/ s<<"\n KsAPC_LT \t"<<c.KsAPC_LT_d;
        /*24*/ s<<"\n Ksi_d \t"<<c.Ksi_d;
        /*25*/ s<<"\n Kst_d \t"<<c.Kst_d;
        /*26*/ s<<"\n APC0_expressing_receptor_d \t"<<c.APC0_expressing_receptor_d;
        /*27*/ s<<"\n u_APC_TNF_d \t"<<c.u_APC_TNF_d;
    }
    return s;


}

APC_cells::APC_cells(const Parameters& p, const Treatment& t):
    /// Variables 6
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
    /// Timidina incorporada
    APC_TymTr_incorporated_d(0.0),

    /// Parámetros 27
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
    /// 6) Proliferation rates
    /*12*/  APC_bound_proliferation_rate_d(p.mean("APC_bound_proliferation_rate")),
    /// 7) Apoptosis rates
    /*13*/  APC0_apop_rate_d(p.mean("APC0_apop_rate")),
    /*14*/  APCa_apop_rate_d(p.mean("APCa_apop_rate")),
    /*15*/  APCbo_apop_rate_d(p.mean("APCbo_apop_rate")),
    /// 8) constant saturation of TNF for apoptosis
    /*16*/  Ks_APC_m_TNF_d(p.mean("Ks_APC_m_TNF")),
    /// 9) conversion rates
    /*17*/  APC_Ag_d(p.mean("APC_Ag")),
    /*18*/  APC_APC_d(p.mean("APC_APC")),
    /*19*/  APC_NK_d(p.mean("APC_NK")),
    /*20*/  APC_LT_1_d(p.mean("APC_LT_1")),
    /*21*/  APC_LT_2_d(p.mean("APC_LT_2")),
    /*22*/  APC_Ab_d(p.mean("APC_Ab")),
    /// 10)Saturation constant of IFN and TNF for activation
    /*23*/  KsAPC_LT_d(p.mean("KsAPC_LT")),
    /// 11)Saturation constant of APC_LT interaction
    /*24*/  Ksi_d(p.mean("APC_Ksi")),
    /*25*/  Kst_d(p.mean("APC_Kst")),
    /// 12) Percentages of cell expressing receptor
    /*26*/  APC0_expressing_receptor_d(p.mean_ratio("APC0_Kratio_expressing_receptor")),
    /// 13) Apoptosis rate for TNF
    /*27*/  u_APC_TNF_d(p.mean("u_APC_TNF"))
{
    if (p.mode()=="minimal")
    {
        IFN_APC0_prod_rate_d=p.mean("IFN_APC_prod_rate");
        IFN_APCa_prod_rate_d=p.mean("IFN_APC_prod_rate");
        IFN_APCbo_prod_rate_d=p.mean("IFN_APC_prod_rate");

    }


}

std::vector<double> APC_cells::Derivative(const Media& m, const NK_cells& NK, const LT_cells& LT)
{
    std::vector<double> D;
    /// we update each subpopulation of cells independently and we take into account the transition from one state to the other

    /// the number of free cells (no Ag) they proliferate according to the cell concentration (factor proliferation ratio)
    /// and some of them are "lost" since they internalize the Ag


    double APC0_delta=(-APC0_apop_rate_d*APC0_d-
                       APC0_d*APC_Ag_d*(m.IFNgamma()/(m.IFNgamma()+ Ksi_d))*(m.TNF()/(m.TNF()+Kst_d))*m.Ag()-
                       APC_Ag_d*m.Ag()*APC0_d);

    D.push_back(APC0_delta);

    /** the cells that have internalize the Ag proliferate in the same way than the free
    they grow also by the free cells that internalize the Ag
    they shrink by the cells that interact with the LT cells*/


    double APCa_delta_pos=(APC0_d*APC_Ag_d*(m.IFNgamma()/(m.IFNgamma()+ Ksi_d))*(m.TNF()/(m.TNF()+Kst_d))*m.Ag()+
                           APC_Ag_d*m.Ag()*APC0_d );
    double APCa_delta_neg= (APCa_apop_rate_d*APCa_d -
                            u_APC_TNF_d*APCa_d*(m.TNF()/(m.TNF()+ Ks_APC_m_TNF_d))-
                            APC_APC_d*APCa_d*(2*APCa_d+APCbo_d)-
                            APC_NK_d*APCa_d*(NK.NKa()*NK.NKa_expressing_receptor()+NK.NKbo())-
                            (APC_LT_1_d-(LT.LT_Ab()*m.Ab()))*LT.LT0()*(APCa_d/(APCa_d+KsAPC_LT_d))-
                            APC_Ab_d*m.Ab()*APCa_d);

    D.push_back(APCa_delta_neg+APCa_delta_pos);

    /// the cells that have interacted with LT grow accordingly with the number of cells that have internalized the Ag and the
    /// number of APC cells expressing the ligand and receptor and bound with monocytes, NK or LT cells with the same rate (Monocytes, NK or LT can
    /// interact only with one cell). We are supposing that all activated cells express receptor and ligand. We are supposing that probabiliities of
    /// interactionts between cells are similar.

    double APCbo_delta=(
                APC_APC_d*APCa_d*(2*APCa_d+APCbo_d)+
                APC_NK_d*APCa_d*(NK.NKa()*NK.NKa_expressing_receptor()+NK.NKbo())+
                (APC_LT_1_d-(LT.LT_Ab()*m.Ab()))*LT.LT0()*(APCa_d/(APCa_d+KsAPC_LT_d))+
                APCbo_d*APC_bound_proliferation_rate_d*m.prol_ratio()-
                APCbo_d*APCbo_apop_rate_d-
                APCbo_d*APC_Ab_d*m.Ab()-
                u_APC_TNF_d*APCbo_d*(m.TNF()/(m.TNF()+ Ks_APC_m_TNF_d))
                );

    D.push_back(APCbo_delta);


    /// the cells that are blocked

    double APCbl_delta= (APC_Ab_d*m.Ab()*APCa_d-
                         APCa_apop_rate_d*APCbl_d-
                         APCbl_d*APC_APC_d*(APCa_d+APCbo_d)+
                         APCbl_d*APC_APC_d*(NK.NKa()+NK.NKbo())+
                         u_APC_TNF_d*APCbl_d*(m.TNF()/(m.TNF()+ Ks_APC_m_TNF_d)));

    D.push_back(APCbl_delta);

    double APCbo_Ab_delta= (APCbo_d*APC_Ab_d*m.Ab()+
                            APCbl_d*APC_APC_d*(APCa_d+APCbo_d)+
                            APCbl_d*APC_APC_d*(NK.NKa()+NK.NKbo())+
                            APCbo_Ab_d*APC_bound_proliferation_rate_d*m.prol_ratio()-
                            APCbo_apop_rate_d*APCbo_Ab_d-
                            APCbo_Ab_d*u_APC_TNF_d*APCbo_Ab_d*(m.TNF()/(m.TNF()+ Ks_APC_m_TNF_d)));
    D.push_back(APCbo_Ab_delta);

    double APC_TymTr_incorporated_delta;
    if (m.TymidineTriteate()>0)
    {
        APC_TymTr_incorporated_delta=
                ((APCbo_d+APCbo_Ab_d)*APC_bound_proliferation_rate_d*m.prol_ratio()*m.Prol_TymTr());
    }
    else
        APC_TymTr_incorporated_delta=0;
    D.push_back(APC_TymTr_incorporated_delta);

    return D;

}

std::vector<double> APC_cells::getState()const
{
    std::vector<double> S;
    S.push_back(APC0_d);
    S.push_back(APCa_d);
    S.push_back(APCbo_d);
    S.push_back(APCbo_Ab_d);
    S.push_back(APCbl_d);
    S.push_back( APC_TymTr_incorporated_d);
    return S;

}

void APC_cells::setState(const std::vector<double>& y)
{
    APC0_d=y[0];
    APCa_d=y[1];
    APCbo_d=y[2];
    APCbl_d+=y[3];
    APCbo_Ab_d=y[4];
    APC_TymTr_incorporated_d=y[5];
}










