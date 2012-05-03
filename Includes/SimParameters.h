#ifndef SIMPARAMETERS_H_INCLUDED
#define SIMPARAMETERS_H_INCLUDED

#include <vector>
#include<string>
#include <iostream>

class SimParameters
{
public:
    std::vector<double> getParameters()const;
    std::vector<double> getRandomParameters(double range)const;
    SimParameters& applyParameters(const std::vector<double>& param);
    std::string mode_;
    /// APC
    /// 1) Init ratio of cells
    /*1*/ double init_ratio_APC_;
    /// 2) IFN Poductions rates of each type of APC
    /*2*/ double IFN_APC0_prod_rate_;
    /*3*/ double IFN_APCa_prod_rate_;
    /*4*/ double IFN_APCbo_prod_rate_;
          double IFN_APC_generic_prod_rate_;
    /// 3) TNF Poductions rates of each type of APC
    /*5*/ double TNF_APC0_prod_rate_;
    /*6*/ double TNF_APCa_prod_rate_;
    /*7*/ double APC_TNF_Induction_CD137;
          double TNF_APC_generic_prod_rate_;
    /// 4) Percentages of IFN productions of each type of APC
    /*8*/ double percentage_IFN_APC0_prod_rate_;
    /*9*/ double percentage_IFN_APCa_prod_rate_;
    //*10*/ double percentage_IFN_APCbo_prod_rate_;
    /// 5)Percentages of TNF productions of each type of APC
    /*11*/ double percentage_TNF_APC0_prod_rate_;
    /// 6) Proliferation rates
    /*12*/ double APC_bound_proliferation_rate_;
    /// 7) Apoptosis rates
    /*13*/ double APC0_apop_rate_;
    /*14*/ double APCa_apop_rate_;
    /*15*/ double APCbo_apop_rate_;
           double APC_generic_apop_rate_;
    /// 8) constant saturation of TNF for apoptosis
    /*16*/ double Ks_APC_m_TNF_;
    /// 9) conversion rates
    /*17*/ double APC_Ag_;
    /*18*/ double APC_APC_;
    /*19*/ double APC_NK_;
    /*20*/ double APC_LT_1_;
    /*21*/ double APC_Ag_2_;
    /*22*/ double APC_Ab_;
//    /// 10)Saturation constant of IFN and TNF for activation
//    /*23*/ double KsAPC_LT_;
    /// 11)Saturation constant of APC_LT interaction
    /*24*/ double APC_Ksi_;
    /*25*/ double APC_Kst_;
    /// 12) Percentages of cell expressing receptor
    /*26*/ double APC0_expressing_receptor_;
    /// 13) Apoptosis rate for TNF
    /*27*/ double u_APC_TNF_;

    /// NK
    /// 1) Init ratio of cells
    /*1*/ double init_ratio_NK_;
    /// 2) IFN Poductions rates of each type of NK
    /*2*/ double IFN_NK0_prod_rate_;
    /*3*/ double IFN_NKa_prod_rate_;
    /*4*/ double IFN_NKbo_prod_rate_;
           double IFN_NK_generic_prod_rate_;
    /// 3) TNF Poductions rates of each type of NK
    /*5*/ double TNF_NK0_prod_rate_;
    /*6*/ double TNF_NKa_prod_rate_;
    /*7*/ double TNF_NKbo_prod_rate_;
           double TNF_NK_generic_prod_rate_;
    /// 4) Percentages of IFN productions of each type of NK
    /*8*/ double percentage_IFN_NK0_prod_rate_;
    /// 5)Percentages of TNF productions of each type of NK
    /*9*/ double percentage_TNF_NK0_prod_rate_;
    /*10*/ double percentage_TNF_NKa_prod_rate_;
    //*11*/ double percentage_TNF_NKbo_prod_rate_;
    /// 6) Proliferation rates
    /*12*/ double NK0_proliferation_rate_;
    /*13*/ double NKa_proliferation_rate_;
    /*14*/ double NK_generic_proliferation_rate_;
    /// 7) Apoptosis rates
    /*15*/ double NK0_apop_rate_;
    /*16*/ double NKa_apop_rate_;
    /*17*/ double NK_generic_apop_rate_;
    /// 8) constant saturation of TNF for apoptosis
    /*18*/ double Ks_NK_m_TNF_;
    /// 9) conversion rates
    /*19*/ double KaNK_;
    /*20*/ double NK_NK_;
    /*21*/ double NK_Ab_;
    /// 10)Saturation constant of NK interaction for activation
    /*22*/ double KsAPC_NK_;
    /// 11)Saturation constant of NK_LT interaction
    /*23*/ double NK_Ksi_;
    /*24*/ double NK_Kst_;
    /// 12) Percentages of cell expressing receptor
    /*25*/ double NK0_expressing_receptor_;
    /*26*/ double NKa_expressing_receptor_;
    /// 13) Apoptosis rate for TNF
    /*27*/ double u_NK_TNF_;

    /// 1) Init number of LT
    /*1*/ double ratio_init_LTns_;
    /*2*/ double ratio_initLTspecific_;
    /// 2) IFN Poductions rates of each type of LT
    /*3*/ double IFN_LTns_prod_rate_;
    /*4*/ double IFN_LTbo_prod_rate_;
    /*5*/ double IFN_LTbl_prod_rate_;
           double IFN_LT_generic_prod_rate_;
    /// 3) TNF Poductions rates of each type of LT
    /*6*/ double TNF_LTns_prod_rate_;
    /*7*/ double TNF_LTbo_prod_rate_;
    /*8*/ double TNF_LTbl_prod_rate_;
           double TNF_LT_generic_prod_rate_;
    /// 4) Percentages of IFN productions of each type of LT
    /*9*/ double percentage_IFN_LTns_prod_rate_;
    /*10*/ double percentage_IFN_LTbo_prod_rate_;
    //*11*/ double percentage_IFN_LTbl_prod_rate_;
    /// 5)Percentages of TNF productions of each type of LT
    /*12*/ double percentage_TNF_LTns_prod_rate_;
    /*13*/ double percentage_TNF_LTbo_prod_rate_;
    //*14*/ double percentage_TNF_LTbl_prod_rate_;
    /// 6) Proliferation rates
    /*15*/ double LTns_proliferation_rate_;
    /*16*/ double LTbo_proliferation_rate_;
    /*17*/ double LTbl_proliferation_rate_;
           double LT_generic_proliferation_rate_;
    /// 7) Apoptosis rates
    /*18*/ double LTns_apop_rate_;
    /*19*/ double LTbo_apop_rate_;
    /*20*/ double LTbl_apop_rate_;
           double LT_generic_apop_rate_;
    /// 8) constant saturation of TNF for apoptosis
    /*21*/ double Ks_LT_m_TNF_;
    /// 9) Percentages of cell expressing receptor
    /*22*/ double LTns_expressing_receptor_;
    /// 10) Apoptosis rate for TNF
    /*23*/ double u_LT_TNF_;
    /// 12) apoptosis related parameters
    /*24*/ double t_apop_meas_;
    /*25*/ double t_duration_apoptosis_;
    double LT_Ab_;

    /// Media
    /*1*/ double TNF_deg_;
    /*2*/ double IFN_deg_;
    /*3*/ double Ag_deg_;
    /*4*/ double TymidineTriteate_;// No es parámetro
    /*5*/ double Prol_TymTr_;
    /*6*/ double max_num_cells_;

    SimParameters(const SimParameters& other);
    friend void swap(SimParameters& one, SimParameters& other);
    friend std::ostream& operator<<(std::ostream& s,SimParameters p);
    SimParameters& operator=(const SimParameters& other);
    SimParameters();
    ~SimParameters(){}
};



#endif // SIMPARAMETERS_H_INCLUDED
