#ifndef LT_H_INCLUDED
#define LT_H_INCLUDED
#include "Includes/SimParameters.h"
#include "Includes/Treatment.h"

class Media;
class APC_cells;
class NK_cells;

/// Lymphocytes T cells
class LT_cells
{
    public:

    ~LT_cells(){}
/// 1) Total LT cells
        double& num_LT();
        const double& num_LT()const;

/// 2) números de células (6)
       /// Number of LT non Ag specific
           double& LTns();
           const double& LTns()const;

       /// Number of naive LT Ag specific cells
           double& LT0();
           const double& LT0()const;

       /// Number of activated LT cells that have been signaled by CD137
           double& LTbo();
           const double& LTbo()const;

       /// Number of activated LT cells that have been blocked for signaling by CD137
           double& LTbl();
           const double& LTbl()const;

       /// Number of LT cells that are exhausted
           double& LTexh ();
           const double& LTexh()const;

/// 3) Percentage of cells exprssing receptor
           double& percentage_cell_expressing_receptor ();
           const double& percentage_cell_expressing_receptor() const;

/// 4) Cytokines production rate and producing cells (6)
    /// Total production of interpheron gamma by LT
        double& LT_IFNgamma_production_rate();
        const double& LT_IFNgamma_production_rate();
    /// Total production of Tumor Necrosis alpha
        double& TNF_production_rate();
        const double& TNF_production_rate() const;
    /// percentage of LT cells that produce TNF
        double& percentage_LT_TNF_production();
        const double& percentage_LT_TNF_production();

/// 5) Percentage of cell expressing the receptor
       double& LT_percentage_cell_expressing_receptor ();
       const double& LT_percentage_cell_expressing_receptor() const;
/// 6) Percentage of LT cell undergoing apoptosis
       double& LT_percentage_cell_undergoing_apoptosis();
       const double& LT_percentage_cell_undergoing_apoptosis();

/// 7) Tymidine incorporated by APC cells
       double& LT_TymTr_incorporated();
       const double& LT_TymTr_incorporated()const;

/// 8) Percentage of LT cells undergoing apoptosis
       double& percentage_apopototic_LT_cells ();
       const double& percentage_apoptotic_LT_cells() const;


        void update(double time_step, double t_run, const Media& m, const APC_cells& a, const NK_cells& NK);

        LT_cells(/// 1) Init number of LT
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
                 /*27*/ double t_duration_apoptosis_,
                 );


        LT_cells(const SimParameters& sp,
                 const Treatment& tr):
            LTns_d(sp.ratio_init_LTns_*tr.init_cells),
            LT0_d(sp.ratio_initLTspecific_*tr.init_cells),
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
            percentage_IFN_LTbo_prod_rate_d(percentage_IFN_LTbo_prod_rate_),
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

            LT_exh_rate_d (LT_exh_rate_),

            t_apop_meas_d (t_apop_meas_),
            t_duration_apoptosis_d(t_duration_apoptosis_)
           {}

            LT_cells(){}

            LT_cells(const LT_cells& other);

            LT_cells& operator=(const LT_cells& other);

            friend void swap(LT_cells& one, LT_cells& other);
	    friend std::ostream& operator<<(std::ostream& s, const LT_cells& c);


            void reset(const SimParameters& sp,const Treatment& tr);

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

        /// number of LT exhausted
        double LTexh_d;

        /// Tymidine incorporated by APC cells
        double LT_TymTr_incorporated_d;

        /// Total LT cell undergoing apoptosis
        double Total_cells_in_apoptosis_d;

        /// Parameters 24
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
            /*11*/ double percentage_IFN_LTbl_prod_rate_d;

        /// 5)Percentages of TNF productions of each type of LT
            /*12*/ double percentage_TNF_LTns_prod_rate_d;
            /*13*/ double percentage_TNF_LTbo_prod_rate_d;
            /*14*/ double percentage_TNF_LTbl_prod_rate_d;

        /// 6) Proliferation rates
            /*15*/ double LTns_proliferation_rate_d;
            /*16*/ double LTbo_proliferation_rate_d;
            /*17*/ double LTbl_proliferation_rate_d;

        /// 7) Apoptosis rates
            /*18*/ double LTns_apop_rate_d;
            /*19*/ double LTbo_apop_rate_d;
            /*20*/ double LTbl_apop_rate_d;
            /*21*/ double LTexh_apop_rate_d;

        /// 8) constant saturation of TNF for apoptosis
            /*22*/ double Ks_LT_m_TNF_d;

        /// 9) Percentages of cell expressing receptor
            /*23*/ double LTns_expressing_receptor_d;

        /// 10) Apoptosis rate for TNF
            /*24*/ double u_LT_TNF_d;

        /// 11) LT exh rate
            /*25*/ double LT_exh_rate_d;

        /// 12) apoptosis related parameters
            /*26*/ double t_apop_meas_d;
            /*27*/ double t_duration_apoptosis_d;


};

#endif // LT_H_INCLUDED
