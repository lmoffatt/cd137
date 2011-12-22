#ifndef NK_H_INCLUDED
#define NK_H_INCLUDED
#include "Includes/SimParameters.h"
#include "Includes/Treatment.h"
#include <iostream>

class Media;
class APC_cells;
class LT_cells;

/// NK cells
class NK_cells
{
    public:
    ~NK_cells(){}
/// 1) Total NK cells
       double num_NK() const;

/// 2) number of cells (6)
       /// Number of native cells
       double& NK0();
       const double& NK0()const;

       /// Number of cells that have bound the antigen and express the receptor
       double& NKa();
       const double& NKa()const;

       /// Number of cells that have express the ligand/receptor and bound
       double& NKbo();
       const double& NKbo() const;

       /// Number of cells that have express the ligand/receptor and bound but then blocked
       double& NKbo_Ab();
       const double& NKbo_Ab() const;

       /// Number of cells that bound to the blocking mAb;
       double& NKbl();
       const double& NKbl()const;

       /// Number of cells that are exhausted
       double& NKexh ();
       const double& NKexh()const;

/// 3) percentage of cells expressing receptor
       double percentage_NK_expressing_receptor () const;

/// 4) Cytokines production rate and producing cells (4)
        /// Total production of interpheron gamma
        double NK_IFNgamma_production_rate() const;

        /// Total production of Tumor Necrosis Factor
        double NK_TNF_production_rate() const;

        /// Percentage of cells producing IFN
        double percentage_NK_producing_IFN () const;

        /// Percentage of cells producing TNF
        double percentage_NK_producing_TNF () const;

        /// 5) Percentage of activated cells expressing receptor
                double& NKa_expressing_receptor();
                const double& NKa_expressing_receptor () const;


        /// 6) Union rates of NK

                /// Ab binding rate
                double& NK_Ab();
                const double& NK_Ab()const;

       /// 7) APC proliferating cells
                double& NK_proliferating_cells();
                const double& NK_proliferating_cells() const;



       /// 8) Tymidine incorporated by APC cells
                 double& NK_TymTr_incorporated();
                 const double& NK_TymTr_incorporated()const;

        NK_cells(const Parameters& p, const Treatment& t);

        void update(double& time_step,const Media& m, const APC_cells& APC,const LT_cells& LT);


        NK_cells ( /// 1) Init number of NK
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
                   /*21*/ double NKexh_apop_rate_,

                   /// 8) constant saturation of TNF for apoptosis
                   /*22*/ double Ks_NK_m_TNF_,

                   /// 9) conversion rates
                   /*23*/ double KaNK_,
                   /*24*/ double NK_NK_,
                   /*25*/ double NK_Ab_,
                   /*26*/ double NK_exh_,

                   /// 10)Saturation constant of NK interaction for activation
                   /*27*/ double KsAPC_NK_,

                   /// 11)Saturation constant of NK_LT interaction
                   /*28*/ double Ksi_,
                   /*29*/ double Kst_,

                   /// 12) Percentages of cell expressing receptor
                   /*30*/ double NK0_expressing_receptor_,
                   /*31*/ double NKa_expressing_receptor_,

                   /// 13) Apoptosis rate for TNF
                   /*32*/ double u_NK_TNF_
                  );
/*
        NK_cells (const SimParameters& sp,
                  const Treatment& tr);
                  */

        NK_cells& operator=(const NK_cells& other);

        friend void swap(NK_cells& one, NK_cells& other);

        NK_cells(const NK_cells& other);

        NK_cells();

     //   void reset(const SimParameters& sp,const Treatment& tr);

        /// main step for the NK cells


	friend std::ostream& operator<<(std::ostream& s, const NK_cells& c);


    private:
      /// Variables 7
         /// number of native cells
             double NK0_d;
         /// number of cells that have internalized the antigen (and therefore express the ligand and receptor)
             double NKa_d;
         /// number of cells that have that have been signaled by receptor or ligand
             double NKbo_d;
         /// number of cells that have that have been signaled by receptor or ligand and bound to the blocking Ab
             double NKbo_Ab_d;
         /// number of cells that binds the blocking mAb
             double NKbl_d;
         /// number of cells that are exhausted
             double NKexh_d;
         /// Timidina incorporada
             double NK_TymTr_incorporated_d;

        /// 1) Init ratio of cells
        /*1*/ double init_ratio_NK_d;

        /// 2) IFN Poductions rates of each type of NK
        /*2*/ double IFN_NK0_prod_rate_d;
        /*3*/ double IFN_NKa_prod_rate_d;
        /*4*/ double IFN_NKbo_prod_rate_d;

        /// 3) TNF Poductions rates of each type of NK
        /*5*/ double TNF_NK0_prod_rate_d;
        /*6*/ double TNF_NKa_prod_rate_d;
        /*7*/ double TNF_NKbo_prod_rate_d;

        /// 4) Percentages of IFN productions of each type of NK
        /*8*/ double percentage_IFN_NK0_prod_rate_d;
        /*9*/ double percentage_IFN_AgNKa_prod_rate_d;
        /*10*/ double percentage_IFN_NKbo_prod_rate_d;

        /// 5)Percentages of TNF productions of each type of NK
        /*11*/ double percentage_TNF_NK0_prod_rate_d;
        /*12*/ double percentage_TNF_NKa_prod_rate_d;
        /*13*/ double percentage_TNF_NKbo_prod_rate_d;

        /// 6) Proliferation rates
        /*13.5*/ double NK0_proliferation_rate_d;
        /*14*/ double NKa_proliferation_rate_d;
        /*15*/ double NKbo_proliferation_rate_d;
        /*16*/ double NKbl_proliferation_rate_d;

        /// 7) Apoptosis rates
        /*17*/ double NK0_apop_rate_d;
        /*18*/ double NKa_apop_rate_d;
        /*19*/ double NKbo_apop_rate_d;
        /*20*/ double NKbl_apop_rate_d;
        /*21*/ double NKexh_apop_rate_d;

        /// 8) constant saturation of TNF for apoptosis
        /*22*/ double Ks_NK_m_TNF_d;

        /// 9) conversion rates
        /*23*/ double KaNK_d;
        /*24*/ double NK_NK_d;
        /*25*/ double NK_Ab_d;
        /*26*/ double NK_exh_d;

        /// 10)Saturation constant of NK interaction for activation
        /*27*/ double KsAPC_NK_d;

        /// 11)Saturation constant of NK_LT interaction
        /*28*/ double Ksi_d;
        /*29*/ double Kst_d;

        /// 12) Percentages of cell expressing receptor
        /*30*/ double NK0_expressing_receptor_d;
        /*31*/ double NKa_expressing_receptor_d;

        /// 13) Apoptosis rate for TNF
        /*32*/ double u_NK_TNF_d;




};
#endif // NK_H_INCLUDED
