#ifndef LT_H_INCLUDED
#define LT_H_INCLUDED
#include "Includes/SimParameters.h"
#include "Includes/Treatment.h"
#include "Includes/Parameters.h"

class Media;
class APC_cells;
class NK_cells;

/// Lymphocytes T cells
class LT_cells
{
public:
    ~LT_cells(){}
    /// 1) Total LT cells
    double num_LT() const;
    /// 2) #LT cells types (4)
    /// # LT non Ag specific
    double& LTns();
    const double& LTns()const;
    /// # naive LT Ag specific cells
    double& LT0();
    const double& LT0()const;
    /// # activated LT cells that have been signaled by CD137
    double& LTbo();
    const double& LTbo()const;
    /// # activated LT cells that have been blocked for signaling by CD137
    double& LTbl();
    const double& LTbl()const;

    /// 3) Percentage of cell expressing the receptor
    double LT_percentage_cell_expressing_receptor () const;

    /// 4) Cytokines production rate and producing cells (4)
    /// Total production of interpheron gamma by LT
    double LT_IFNgamma_production_rate() const;
    /// Total production of Tumor Necrosis alpha
    double TNF_production_rate() const;
    /// percentage of LT cells that produce TNF
    double percentage_LT_IFN_production() const;
    /// percentage of LT cells that produce TNF
    double percentage_LT_TNF_production() const;


    /// 5) Tymidine incorporated by LT cells
    double& LT_TymTr_incorporated();
    const double& LT_TymTr_incorporated()const;

    /// 6) Percentage of LT cells undergoing apoptosis
    double percentage_apoptotic_LT_cells () const;

    /// 7) LT-Ab binding rate
    double& LT_Ab();
    const double& LT_Ab() const;



    void update(double& time_step, double t_run, const Media& m, const APC_cells& APC, const NK_cells& NK);

    LT_cells(const Parameters& p, const Treatment& t);

    LT_cells(/// 1) Init number of LT
             /*1*/ double init_number_LTns_,
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
             /// 13) LT-Ab binding rate
             /*26*/ double LT_Ab_
             );

    LT_cells(){}
    LT_cells(const LT_cells& other);
    LT_cells& operator=(const LT_cells& other);
    friend void swap(LT_cells& one, LT_cells& other);
    friend std::ostream& operator<<(std::ostream& s, const LT_cells& c);
    std::vector<double> Derivative(double t_run, const Media& m, const APC_cells& APC)const;
    std::vector<double> getState()const;
    void setState(const std::vector<double>& y);




private:
    /// Variables 6
    /// number of non Ag specific cells
    double LTns_d;
    /// number of naive Ag specific cells
    double LT0_d;
    /// number of Ag specific cells that have recieve receptor signaling during sinapsis
    double LTbo_d;
    /// number of Ag specific cells that have not recieve receptor singaling during sinapsis
    double LTbl_d;
    /// Tymidine incorporated by APC cells
    double LT_TymTr_incorporated_d;
    /// Total LT cell undergoing apoptosis
    double Total_cells_in_apoptosis_d;

    /// Parameters 26
    /// 1) Init number of LT
    /*1*/ double ratio_init_LTns_d;
    /*2*/ double ratio_initLTspecific_d;
    /// 2) IFN Poductions rates of each type of LT
    /*3*/ double IFN_LTns_prod_rate_d;
    /*4*/ double IFN_LTbo_prod_rate_d;
    /*5*/ double IFN_LTbl_prod_rate_d;
    /// 3) TNF Poductions rates of each type of LT
    /*6*/ double TNF_LTns_prod_rate_d;
    /*7*/ double TNF_LTbo_prod_rate_d;
    /*8*/ double TNF_LTbl_prod_rate_d;
    /// 4) Percentages of IFN productions of each type of LT
    /*9*/ double percentage_IFN_LTns_prod_rate_d;
    /*10*/ double percentage_IFN_LTbo_prod_rate_d;
    //*11*/ double percentage_IFN_LTbl_prod_rate_d;
    /// 5)Percentages of TNF productions of each type of LT
    /*12*/ double percentage_TNF_LTns_prod_rate_d;
    /*13*/ double percentage_TNF_LTbo_prod_rate_d;
    //*14*/ double percentage_TNF_LTbl_prod_rate_d;
    /// 6) Proliferation rates
    /*15*/ double LTns_proliferation_rate_d;
    /*16*/ double LTbo_proliferation_rate_d;
    /*17*/ double LTbl_proliferation_rate_d;
    /// 7) Apoptosis rates
    /*18*/ double LTns_apop_rate_d;
    /*19*/ double LTbo_apop_rate_d;
    /*20*/ double LTbl_apop_rate_d;
    /// 8) constant saturation of TNF for apoptosis
    /*21*/ double Ks_LT_m_TNF_d;
    /// 9) Percentages of cell expressing receptor
    /*22*/ double LTns_expressing_receptor_d;
    /// 10) Apoptosis rate for TNF
    /*23*/ double u_LT_TNF_d;
    /// 12) apoptosis related parameters
    /*24*/ double t_apop_meas_d;
    /*25*/ double t_duration_apoptosis_d;
    /// 13) LT-Ab binding rate
    /*26*/ double LT_Ab_d;


};

#endif // LT_H_INCLUDED
