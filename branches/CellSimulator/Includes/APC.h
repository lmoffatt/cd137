#ifndef APC_H_INCLUDED
#define APC_H_INCLUDED
#include "Includes/SimParameters.h"
#include "Includes/Treatment.h"
#include "Includes/Parameters.h"


class Media;
class NK_cells;
class LT_cells;



/// Antigen presentation cells (Monocytes)
class APC_cells
{
    public:
 /// 1) Total APC cells
        double& num_APC();
        const double& num_APC()const;


/// 2) números de células (6)
        /// Number of naive cells
        double& APC0();
        const double& APC0()const;

        /// Number of cells that have bound the antigen and express the receptor and ligand
        double& APCa();
        const double& APCa()const;

        /// Number of cells that have express the receptor and ligand and bound to ligand/receptor
        double& APCbo();
        const double& APCbo()const;

        /// Number of cells that have express the receptor and ligand and bound to ligand/receptor and boundo the antibody
        double& APCbo_Ab();
        const double& APCbo_Ab()const;

        /// Number of cells that express the receptor and bind the mAb
        double& APCbl();
        const double& APCbl()const;

        /// Number of cells that are exhausted
        double& APCexh ();
        const double& APCexh()const;

/// 3) Percentage of cells exprssing receptor
        double& percentage_cell_expressing_receptor ();
        const double& percentage_cell_expressing_receptor() const;

/// 4) Cytokines production rate and producing cells (6)
        /// Total production of interpheron gamma
        double& APC_IFNgamma_production_rate();
        const double& APC_IFNgamma_production_rate()  const;

        /// Total production of Tumor Necrosis Factor
        double& APC_TNF_production_rate();
        const double& APC_TNF_production_rate()const;

        /// Percentage of cells producing IFN
        double& percentage_APC_producing_IFN();
        const double& percentage_APC_producing_IFN() const;

        /// Percentage of cells producing TNF
        double& percentage_APC_producing_TNF ();
        const double& percentage_APC_producing_TNF() const;

        /// APCa TNF production rate
        double& APCa_TNF_production_rate ();
        const double& APCa_TNF_production_rate () const;

        /// APCbo TNF production rate
        double& APCbo_TNF_production_rate ();
        const double& APCbo_TNF_production_rate () const;


/// 5) Percentage of activated cells expressing receptor
        double& APCa_expressing_receptor();
        const double& APCa_expressing_receptor () const;


/// 6) Union rates of APC (5)
        /// Ag internalization rate
        double& APC_Ag();
        const double& APC_Ag()const;

        /// Receptor binding rate to  NK
        double& APC_NK();
        const double& APC_NK()const;

        /// Receptor binding rate to LT
        double& APC_LT_1();
        const double& APC_LT_1() const;

        double& APC_LT_2();
        const double& APC_LT_2() const;

        /// Ab binding rate
        double& APC_Ab();
        const double& APC_Ab()const;


        /// Constante de saturación Unión APC_LT
        double& KsAPC_LT();
        const double& KsAPC_LT () const;

  /// 8) Tymidine incorporated by APC cells
        double& APC_TymTr_incorporated();
        const double& APC_TymTr_incorporated()const;

        void update(double time_step,const Media& m, const NK_cells& NK,const LT_cells& LT);


        APC_cells(const Parameters& p, const Treatment& t);

        APC_cells(
        /// 1) Init number of APC
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
        /*9*/ double percentage_IFN_AgAPCa_prod_rate_,
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

                          );

        APC_cells(const APC_cells& other);

        friend void swap(APC_cells& one, APC_cells& other);

        APC_cells& operator=(const APC_cells& other);

  /*      APC_cells(const SimParameters& sp,
                  const Treatment& tr);*/

        APC_cells();
        ~APC_cells(){};

     //   void reset(const SimParameters& sp,const Treatment& tr);

        /// main step for the APC cells

	friend std::ostream& operator<<(std::ostream& s, const APC_cells& c);

    private:


        /// Variables 7
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
        /// number of cells that are exhausted
        double APCexh_d;
        /// Timidina incorporada
        double APC_TymTr_incorporated_d;

        /// Parámetros 33
        /// 1) Init ratio of cells
        /*1*/ double init_ratio_APC_d;

        /// 2) IFN Poductions rates of each type of APC
        /*2*/ double IFN_APC0_prod_rate_d;
        /*3*/ double IFN_APCa_prod_rate_d;
        /*4*/ double IFN_APCbo_prod_rate_d;


        /// 3) TNF Poductions rates of each type of APC
        /*5*/ double TNF_APC0_prod_rate_d;
        /*6*/ double TNF_APCa_prod_rate_d;
        /*7*/ double TNF_APCbo_prod_rate_d;


        /// 4) Percentages of IFN productions of each type of APC
        /*8*/ double percentage_IFN_APC0_prod_rate_d;
        /*9*/ double percentage_IFN_APCa_prod_rate_d;
        /*10*/ double percentage_IFN_APCbo_prod_rate_d;

        /// 5)Percentages of TNF productions of each type of APC
        /*11*/ double percentage_TNF_APC0_prod_rate_d;
        /*12*/ double percentage_TNF_APCa_prod_rate_d;
        /*13*/ double percentage_TNF_APCbo_prod_rate_d;


        /// 6) Proliferation rates
        /*14*/ double APC_bound_proliferation_rate_d;

        /// 7) Apoptosis rates
        /*15*/ double APC0_apop_rate_d;
        /*16*/ double APCa_apop_rate_d;
        /*17*/ double APCbo_apop_rate_d;
        /*18*/ double APCbl_apop_rate_d;
        /*19*/ double APCexh_apop_rate_d;

        /// 8) constant saturation of TNF for apoptosis
        /*20*/ double Ks_APC_m_TNF_d;

        /// 9) conversion rates
        /*21*/ double APC_Ag_d;
        /*22*/ double APC_APC_d;
        /*23*/ double APC_NK_d;
        /*24*/ double APC_LT_1_d;
        /*25*/ double APC_LT_2_d;
        /*26*/ double APC_Ab_d;
        /*27*/ double APC_exh_d;

        /// 10)Saturation constant of IFN and TNF for activation
        /*28*/ double KsAPC_LT_d;

        /// 11)Saturation constant of APC_LT interaction
        /*29*/ double Ksi_d;
        /*30*/ double Kst_d;

        /// 12) Percentages of cell expressing receptor
        /*31*/ double APC0_expressing_receptor_d;
        /*32*/ double APCa_expressing_receptor_d;

        /// 13) Apoptosis rate for TNF
        /*33*/ double u_APC_TNF_d;



};



#endif // APC_H_INCLUDED
