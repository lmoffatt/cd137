#ifndef APC_H_INCLUDED
#define APC_H_INCLUDED
#include "Includes/SimParameters.h"
#include "Includes/Treatment.h"
#include "Includes/Parameters.h"

class Media;
class NK_cells;
class LT_cells;

class APC_cells
{
public:
    /// 1) Total APC cells
    double num_APC() const;

    /// 2) APC types (5)
    /// Number of naive cells
    double& APC0();
    const double& APC0()const;
    /// #APC cells that have internalized the antigen and express the receptor and ligand
    double& APCa();
    const double& APCa()const;
    /// #APC cells that have internalized the antigen, express receptor and ligand and have been signalized by reverse signaling
    double& APCbo();
    const double& APCbo()const;
    /// #APC cells that are blocked but had been previously signalized
    double& APCbo_Ab();
    const double& APCbo_Ab()const;
    /// #APC cells that express the receptor and bind the mAb
    double& APCbl();
    const double& APCbl()const;

    /// 3) Percentage of APC cells exprssing receptor
    double percentage_cell_expressing_receptor () const;

    /// 4) Cytokines production rate and producing cells (6)
    /// Total production of interpheron gamma
    double APC_IFNgamma_production_rate() const;
    /// Total production of Tumor Necrosis Factor
    double APC_TNF_production_rate() const;
    /// Percentage of cells producing IFN
    double percentage_APC_producing_IFN() const;
    /// Percentage of cells producing TNF
    double percentage_APC_producing_TNF () const;
    /// APCa TNF production rate
    double& APCa_TNF_production_rate ();
    const double& APCa_TNF_production_rate () const;
    /// APCbo TNF production rate
    double& APCbo_TNF_production_rate ();
    const double& APCbo_TNF_production_rate () const;


    bool withinThreshold(const APC_cells& other,double threshold)const;

    /// 5) Union rates of APC (5)
    /// Ag internalization rate
    double& APC_Ag();
    const double& APC_Ag()const;
    /// Receptor binding rate to  NK
    double& APC_NK();
    const double& APC_NK()const;
    /// Receptor binding rate to LT
    double& APC_LT_1();
    const double& APC_LT_1() const;
    double& APC_Ag_2();
    const double& APC_Ag_2() const;
    /// Ab binding rate
    double& APC_Ab();
    const double& APC_Ab()const;
//    /// Constante de saturación Unión APC_LT
//    double& KsAPC_LT();
//    const double& KsAPC_LT() const;

    /// 6) Tymidine incorporated by APC cells
    double& APC_TymTr_incorporated();
    const double& APC_TymTr_incorporated()const;

    void update(double& time_step,const Media& m, const NK_cells& NK,const LT_cells& LT);

    APC_cells(const Parameters& p, const Treatment& t);

    APC_cells(
        /// 1) Init number of APC
        /*1*/ double init_APC_,
        /// 2) IFN Poductions rates of each type of APC
        /*2*/ double IFN_APC0_prod_rate_,
        /*3*/ double IFN_APCa_prod_rate_,
        /*4*/ double IFN_APCbo_prod_rate_,
              double IFN_APC_generic_prod_rate_,
        /// 3) TNF Poductions rates of each type of APC
        /*5*/ double TNF_APC0_prod_rate_,
        /*6*/ double TNF_APCa_prod_rate_,
        /*7*/ double TNF_APCbo_prod_rate_,
              double TNF_APC_generic_prod_rate_,
        /// 4) Percentages of IFN productions of each type of APC
        /*8*/ double percentage_IFN_APC0_prod_rate_,
        /*9*/ double percentage_IFN_APCa_prod_rate_,
        //*10*/ double percentage_IFN_APCbo_prod_rate_,
        /// 5)Percentages of TNF productions of APC0
        /*11*/ double percentage_TNF_APC0_prod_rate_,
        /// 6) Proliferation rates
        /*12*/ double APC_bound_proliferation_rate_,
        /// 7) Apoptosis rates
        /*13*/ double APC0_apop_rate_,
        /*14*/ double APCa_apop_rate_,
        /*15*/ double APCbo_apop_rate_,
               double APC_generic_apop_rate_,
        /// 8) constant saturation of TNF for apoptosis
        /*16*/ double Ks_APC_m_TNF_,
        /// 9) conversion rates
        /*17*/ double APC_Ag_,
        /*18*/ double APC_APC_,
        /*19*/ double APC_NK_,
        /*20*/ double APC_LT_1_,
        /*21*/ double APC_Ag_2_,
        /*22*/ double APC_Ab_,
//        /// 10)Saturation constant of IFN and TNF for activation
//        /*23*/ double KsAPC_LT_,
        /// 11)Saturation constant of APC_LT interaction
        /*24*/ double Ksi_,
        /*25*/ double Kst_,
        /// 12) Percentages of cell expressing receptor
        /*26*/ double APC0_expressing_receptor_,
        /// 13) Apoptosis rate for TNF
        /*27*/ double u_APC_TNF_
        );

    APC_cells(const APC_cells& other);
    friend void swap(APC_cells& one, APC_cells& other);
    APC_cells& operator=(const APC_cells& other);
    APC_cells();
    ~APC_cells(){};


    friend std::ostream& operator<<(std::ostream& s, const APC_cells& c);
    std::vector<double> Derivative(const Media& m, const NK_cells& NK, const LT_cells& LT);
    std::vector<double> getState()const;
    void setState(const std::vector<double>& y);



private:


    /// Variables 6
    /// number of native cells
    double APC0_d;
    /// number of cells that have internalized the antigen (and therefore express the ligand and receptor)
    double APCa_d;
    /// number of cells that have that have been signaled by receptor or ligand
    double APCbo_d;
    /// number of cells that have that have been signaled by receptor or ligand and bound to the blocking Ab
    double APCbo_Ab_d;
    /// number of cells that binds the blocking mAb
    double APCbl_d;
    /// Timidina incorporada
    double APC_TymTr_incorporated_d;

    /// Parámetros 27
    /// 1) Init ratio of cells
    /*1*/ double init_ratio_APC_d;
    /// 2) IFN Poductions rates of each type of APC
    /*2*/ double IFN_APC0_prod_rate_d;
    /*3*/ double IFN_APCa_prod_rate_d;
    /*4*/ double IFN_APCbo_prod_rate_d;
          double IFN_APC_generic_prod_rate_d;
    /// 3) TNF Poductions rates of each type of APC
    /*5*/ double TNF_APC0_prod_rate_d;
    /*6*/ double TNF_APCa_prod_rate_d;
    /*7*/ double TNF_APCbo_prod_rate_d;
          double TNF_APC_generic_prod_rate_d;
    /// 4) Percentages of IFN productions of each type of APC
    /*8*/ double percentage_IFN_APC0_prod_rate_d;
    /*9*/ double percentage_IFN_APCa_prod_rate_d;
    //*10*/ double percentage_IFN_APCbo_prod_rate_d;
    /// 5)Percentages of TNF productions of each type of APC
    /*11*/ double percentage_TNF_APC0_prod_rate_d;
    /// 6) Proliferation rates
    /*12*/ double APC_bound_proliferation_rate_d;
    /// 7) Apoptosis rates
    /*13*/ double APC0_apop_rate_d;
    /*14*/ double APCa_apop_rate_d;
    /*15*/ double APCbo_apop_rate_d;
           double APC_generic_apop_rate_d;
    /// 8) constant saturation of TNF for apoptosis
    /*16*/ double Ks_APC_m_TNF_d;
    /// 9) conversion rates
    /*17*/ double APC_Ag_d;
    /*18*/ double APC_APC_d;
    /*19*/ double APC_NK_d;
    /*20*/ double APC_LT_1_d;
    /*21*/ double APC_Ag_2_d;
    /*22*/ double APC_Ab_d;
    /// 10)Saturation constant of IFN and TNF for activation
    /*23*/ double KsAPC_LT_d;
    /// 11)Saturation constant of APC_LT interaction
    /*24*/ double Ksi_d;
    /*25*/ double Kst_d;
    /// 12) Percentages of cell expressing receptor
    /*26*/ double APC0_expressing_receptor_d;
    /// 13) Apoptosis rate for TNF
    /*27*/ double u_APC_TNF_d;
};



#endif // APC_H_INCLUDED
