#include "Includes/SimParameters.h"
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cstring>

SimParameters::SimParameters():
    mode_("FULL"),
    /// APC
    /// 1) Init ratio of cells
    /*1*/ init_ratio_APC_(1e5/1e6),

    /// 2) IFN Poductions rates of each type of APC
    /*2*/ IFN_APC0_prod_rate_(1e-10),
    /*3*/ IFN_APCa_prod_rate_(1e-10),
    /*4*/ IFN_APCbo_prod_rate_(1e-10),


    /// 3) TNF Poductions rates of each type of APC
    /*5*/ TNF_APC0_prod_rate_(1e-10),
    /*6*/ TNF_APCa_prod_rate_(1e-10),
    /*7*/ TNF_APCbo_prod_rate_(1e-10),


    /// 4) Percentages of IFN productions of each type of APC
    /*8*/ percentage_IFN_APC0_prod_rate_(0.01),
    /*9*/ percentage_IFN_APCa_prod_rate_(0.01),
    /*10*/ percentage_IFN_APCbo_prod_rate_(0.01),

    /// 5)Percentages of TNF productions of each type of APC
    /*11*/ percentage_TNF_APC0_prod_rate_(0.01),
    /*12*/ percentage_TNF_APCa_prod_rate_(0.01),
    /*13*/ percentage_TNF_APCbo_prod_rate_(0.01),


    /// 6) Proliferation rates
    /*14*/ APC_bound_proliferation_rate_(1e-6),

    /// 7) Apoptosis rates
    /*15*/ APC0_apop_rate_(1e-6),
    /*16*/ APCa_apop_rate_(1e-6),
    /*17*/ APCbo_apop_rate_(1e-6),
    /*18*/ APCbl_apop_rate_(1e-6),
    /*19*/ APCexh_apop_rate_(1e-6),

    /// 8) constant saturation of TNF for apoptosis
    /*20*/ Ks_APC_m_TNF_(1e-6),

    /// 9) conversion rates
    /*21*/ APC_Ag_(1e-6),
    /*22*/ APC_APC_(1e-6),
    /*23*/ APC_NK_(1e-6),
    /*24*/ APC_LT_1_(1e-6),
    /*25*/ APC_LT_2_(1e-6),
    /*26*/ APC_Ab_(1e-6),
    /*27*/ APC_exh_(1e-6),

    /// 10)Saturation constant of IFN and TNF for activation
    /*28*/ KsAPC_LT_(1e-6),

    /// 11)Saturation constant of APC_LT interaction
    /*29*/ APC_Ksi_(1e-6),
    /*30*/ APC_Kst_(1e-6),

    /// 12) Percentages of cell expressing receptor
    /*31*/ APC0_expressing_receptor_(0.01),
    /*32*/ APCa_expressing_receptor_(0.01),
    /// 13) Apoptosis rate for TNF
    /*33*/ u_APC_TNF_ (1e-12),

    /// NK
    /// 1) Init ratio of cells
    /*1*/  init_ratio_NK_(1e5/1e6),

    /// 2) IFN Poductions rates of each type of NK
    /*2*/  IFN_NK0_prod_rate_(1e-10),
    /*3*/  IFN_NKa_prod_rate_(1e-10),
    /*4*/  IFN_NKbo_prod_rate_(1e-10),

    /// 3) TNF Poductions rates of each type of NK
    /*5*/  TNF_NK0_prod_rate_(1e-10),
    /*6*/  TNF_NKa_prod_rate_(1e-10),
    /*7*/  TNF_NKbo_prod_rate_(1e-10),

    /// 4) Percentages of IFN productions of each type of NK
    /*8*/  percentage_IFN_NK0_prod_rate_(0.01),
    /*9*/  percentage_IFN_AgNKa_prod_rate_(0.01),
    /*10*/  percentage_IFN_NKbo_prod_rate_(0.01),

    /// 5)Percentages of TNF productions of each type of NK
    /*11*/  percentage_TNF_NK0_prod_rate_(0.01),
    /*12*/  percentage_TNF_NKa_prod_rate_(0.01),
    /*13*/  percentage_TNF_NKbo_prod_rate_(0.01),

    /// 6) Proliferation rates
    /*13.5*/  NK0_proliferation_rate_(1e-6),
    /*14*/  NKa_proliferation_rate_(1e-6),
    /*15*/  NKbo_proliferation_rate_(1e-6),
    /*16*/  NKbl_proliferation_rate_(1e-6),

    /// 7) Apoptosis rates
    /*17*/  NK0_apop_rate_(1e-6),
    /*18*/  NKa_apop_rate_(1e-6),
    /*19*/  NKbo_apop_rate_(1e-6),
    /*20*/  NKbl_apop_rate_(1e-6),
    /*21*/  NKexh_apop_rate_(1e-6),

    /// 8) constant saturation of TNF for apoptosis
    /*22*/  Ks_NK_m_TNF_(1e-6),

    /// 9) conversion rates
    /*23*/  KaNK_(1e-6),
    /*24*/  NK_NK_(1e-6),
    /*25*/  NK_Ab_(1e-6),
    /*26*/  NK_exh_(1e-6),

    /// 10)Saturation constant of APC NK interaction for activation
    /*27*/  KsAPC_NK_(1e-6),

    /// 11)Saturation constant of NK_LT interaction
    /*28*/  NK_Ksi_(1e-6),
    /*29*/  NK_Kst_(1e-6),

    /// 12) Percentages of cell expressing receptor
    /*30*/  NK0_expressing_receptor_(0.01),
    /*31*/  NKa_expressing_receptor_(0.01),

    /// 13) Apoptosis rate for TNF
    /*32*/  u_NK_TNF_(1e-6),

    /// LT
    /// 1) Init number of LT
       /*1*/  ratio_init_LTns_(0.9e6/1e6),
       /*2*/  ratio_initLTspecific_(1e2/1e6),

    /// 2) IFN Poductions rates of each type of LT
       /*3*/  IFN_LTns_prod_rate_(1e-10),
       /*4*/  IFN_LTbo_prod_rate_(1e-10),
       /*5*/  IFN_LTbl_prod_rate_(1e-10),

   /// 3) TNF Poductions rates of each type of LT
       /*6*/  TNF_LTns_prod_rate_(1e-10),
       /*7*/  TNF_LTbo_prod_rate_(1e-10),
       /*8*/  TNF_LTbl_prod_rate_(1e-10),

   /// 4) Percentages of IFN productions of each type of LT
       /*9*/  percentage_IFN_LTns_prod_rate_(0.05),
       /*10*/  percentage_IFN_LTbo_prod_rate_(0.8),
       /*11*/  percentage_IFN_LTbl_prod_rate_(0.8),


   /// 5)Percentages of TNF productions of each type of LT
       /*12*/  percentage_TNF_LTns_prod_rate_(0.05),
       /*13*/  percentage_TNF_LTbo_prod_rate_(0.8),
       /*14*/  percentage_TNF_LTbl_prod_rate_(0.8),

   /// 6) Proliferation rates
       /*15*/  LTns_proliferation_rate_(1e-10),
       /*16*/  LTbo_proliferation_rate_(1e-10),
       /*17*/  LTbl_proliferation_rate_(1e-10),

   /// 7) Apoptosis rates
       /*18*/  LTns_apop_rate_(1e-10),
       /*19*/  LTbo_apop_rate_(1e-10),
       /*20*/  LTbl_apop_rate_(1e-10),
       /*21*/  LTexh_apop_rate_(1e-10),

   /// 8) constant saturation of TNF for apoptosis
       /*22*/  Ks_LT_m_TNF_(1e-10),

   /// 9) Percentages of cell expressing receptor
       /*23*/  LTns_expressing_receptor_(1e-10),

   /// 10) Apoptosis rate for TNF
       /*24*/  u_LT_TNF_(1e-10),

   /// 11) LT exh rate
       /*25*/ LT_exh_rate_(1e-10),

   /// 12) apoptosis related parameters
       /*26*/ t_apop_meas_(120),
       /*27*/ t_duration_apoptosis_(2),

    /// Media
    /*1*/ TNF_deg_(0.5),
    /*2*/ IFN_deg_(0.5),
    /*3*/ TymidineTriteate_(0.5),
    /*4*/ Prol_TymTr_(0.5)

{}

std::vector<double> SimParameters::getParameters()const
{
    std::vector<double> par;

    if (mode_=="FULL")
    {
   /// APC
        // multiplicar si estoy seguro del valor
        /// 1) Init ratio of cells
        /*1*/   par.push_back(log(init_ratio_APC_));
        /// 2) IFN Poductions rates of each type of APC
        /*2*/ par.push_back(log(IFN_APC0_prod_rate_));
        /*3*/ par.push_back(log(IFN_APCa_prod_rate_));
        /*4*/ par.push_back(log(IFN_APCbo_prod_rate_));


        /// 3) TNF Poductions rates of each type of APC
        /*5*/ par.push_back(log(TNF_APC0_prod_rate_));
        /*6*/ par.push_back(log(TNF_APCa_prod_rate_));
        /*7*/ par.push_back(log(TNF_APCbo_prod_rate_));


        /// 4) Percentages of IFN productions of each type of APC
        /*8*/ par.push_back(log(percentage_IFN_APC0_prod_rate_));
        /*9*/ par.push_back(log(percentage_IFN_APCa_prod_rate_));
        /*10*/ par.push_back(log(percentage_IFN_APCbo_prod_rate_));

        /// 5)Percentages of TNF productions of each type of APC
        /*11*/ par.push_back(log(percentage_TNF_APC0_prod_rate_));
        /*12*/ par.push_back(log(percentage_TNF_APCa_prod_rate_));
        /*13*/ par.push_back(log(percentage_TNF_APCbo_prod_rate_));


        /// 6) Proliferation rates
        /*14*/ par.push_back(log(APC_bound_proliferation_rate_));

        /// 7) Apoptosis rates
        /*15*/ par.push_back(log(APC0_apop_rate_));
        /*16*/ par.push_back(log(APCa_apop_rate_));
        /*17*/ par.push_back(log(APCbo_apop_rate_));
        /*18*/ par.push_back(log(APCbl_apop_rate_));
        /*19*/ par.push_back(log(APCexh_apop_rate_));

        /// 8) constant saturation of TNF for apoptosis
        /*20*/ par.push_back(log(Ks_APC_m_TNF_));

        /// 9) conversion rates
        /*21*/ par.push_back(log(APC_Ag_));
        /*22*/ par.push_back(log(APC_APC_));
        /*23*/ par.push_back(log(APC_NK_));
        /*24*/ par.push_back(log(APC_LT_1_));
        /*25*/ par.push_back(log(APC_LT_2_));
        /*26*/ par.push_back(log(APC_Ab_));
        /*27*/ par.push_back(log(APC_exh_));

        /// 10)Saturation constant of IFN and TNF for activation
        /*28*/ par.push_back(log(KsAPC_LT_));

        /// 11)Saturation constant of APC_LT interaction
        /*29*/ par.push_back(log(APC_Ksi_));
        /*30*/ par.push_back(log(APC_Kst_));

        /// 12) Percentages of cell expressing receptor
        /*31*/ par.push_back(log(APC0_expressing_receptor_));
        /*32*/ par.push_back(log(APCa_expressing_receptor_));
        /// 13) Apoptosis rate for TNF
        /*33*/ par.push_back(log(u_APC_TNF_));

   /// NK
        /// 1) Init ratio of cells
        /*1*/  par.push_back(log(init_ratio_NK_));

        /// 2) IFN Poductions rates of each type of NK
        /*2*/  par.push_back(log(IFN_NK0_prod_rate_));
        /*3*/  par.push_back(log(IFN_NKa_prod_rate_));
        /*4*/  par.push_back(log(IFN_NKbo_prod_rate_));

        /// 3) TNF Poductions rates of each type of NK
        /*5*/  par.push_back(log(TNF_NK0_prod_rate_));
        /*6*/  par.push_back(log(TNF_NKa_prod_rate_));
        /*7*/  par.push_back(log(TNF_NKbo_prod_rate_));

        /// 4) Percentages of IFN productions of each type of NK
        /*8*/  par.push_back(log(percentage_IFN_NK0_prod_rate_));
        /*9*/  par.push_back(log(percentage_IFN_AgNKa_prod_rate_));
        /*10*/  par.push_back(log(percentage_IFN_NKbo_prod_rate_));

        /// 5)Percentages of TNF productions of each type of NK
        /*11*/  par.push_back(log(percentage_TNF_NK0_prod_rate_));
        /*12*/  par.push_back(log(percentage_TNF_NKa_prod_rate_));
        /*13*/  par.push_back(log(percentage_TNF_NKbo_prod_rate_));

        /// 6) Proliferation rates
        /*13.5*/par.push_back (log(NK0_proliferation_rate_));
        /*14*/  par.push_back(log(NKa_proliferation_rate_));
        /*15*/  par.push_back(log(NKbo_proliferation_rate_));
        /*16*/  par.push_back(log(NKbl_proliferation_rate_));

        /// 7) Apoptosis rates
        /*17*/  par.push_back(log(NK0_apop_rate_));
        /*18*/  par.push_back(log(NKa_apop_rate_));
        /*19*/  par.push_back(log(NKbo_apop_rate_));
        /*20*/  par.push_back(log(NKbl_apop_rate_));
        /*21*/  par.push_back(log(NKexh_apop_rate_));

        /// 8) constant saturation of TNF for apoptosis
        /*22*/  par.push_back(log(Ks_NK_m_TNF_));

        /// 9) conversion rates
        /*23*/  par.push_back(log(KaNK_));
        /*24*/  par.push_back(log(NK_NK_));
        /*25*/  par.push_back(log(NK_Ab_));
        /*26*/  par.push_back(log(NK_exh_));

        /// 10)Saturation constant of NK interaction for activation
        /*27*/  par.push_back(log(KsAPC_NK_));

        /// 11)Saturation constant of NK_LT interaction
        /*28*/  par.push_back(log(NK_Ksi_));
        /*29*/  par.push_back(log(NK_Kst_));

        /// 12) Percentages of cell expressing receptor
        /*30*/  par.push_back(log(NK0_expressing_receptor_));
        /*31*/  par.push_back(log(NKa_expressing_receptor_));

        /// 13) Apoptosis rate for TNF
        /*32*/  par.push_back(log(u_NK_TNF_));

   /// LT
       /// 1) Init number of LT
       /*1*/  par.push_back(log(ratio_init_LTns_));
       /*2*/  par.push_back(log(ratio_initLTspecific_));

       /// 2) IFN Poductions rates of each type of LT
       /*3*/  par.push_back(log(IFN_LTns_prod_rate_));
       /*4*/  par.push_back(log(IFN_LTbo_prod_rate_));
       /*5*/  par.push_back(log(IFN_LTbl_prod_rate_));

       /// 3) TNF Poductions rates of each type of LT
       /*6*/  par.push_back(log(TNF_LTns_prod_rate_));
       /*7*/  par.push_back(log(TNF_LTbo_prod_rate_));
       /*8*/  par.push_back(log(TNF_LTbl_prod_rate_));

       /// 4) Percentages of IFN productions of each type of LT
       /*9*/  par.push_back(log(percentage_IFN_LTns_prod_rate_));
       /*10*/  par.push_back(log(percentage_IFN_LTbo_prod_rate_));
       /*11*/  par.push_back(log(percentage_IFN_LTbl_prod_rate_));

       /// 5)Percentages of TNF productions of each type of LT
       /*12*/  par.push_back(log(percentage_TNF_LTns_prod_rate_));
       /*13*/  par.push_back(log(percentage_TNF_LTbo_prod_rate_));
       /*14*/  par.push_back(log(percentage_TNF_LTbl_prod_rate_));

       /// 6) Proliferation rates
       /*15*/  par.push_back(log(LTns_proliferation_rate_));
       /*16*/  par.push_back(log(LTbo_proliferation_rate_));
       /*17*/  par.push_back(log(LTbl_proliferation_rate_));

       /// 7) Apoptosis rates
       /*18*/  par.push_back(log(LTns_apop_rate_));
       /*19*/  par.push_back(log(LTbo_apop_rate_));
       /*20*/  par.push_back(log(LTbl_apop_rate_));
       /*21*/  par.push_back(log(LTexh_apop_rate_));

       /// 8) constant saturation of TNF for apoptosis
       /*22*/  par.push_back(log(Ks_LT_m_TNF_));

       /// 9) Percentages of cell expressing receptor
       /*23*/  par.push_back(log(LTns_expressing_receptor_));

       /// 10) Apoptosis rate for TNF
       /*24*/  par.push_back(log(u_LT_TNF_));

       /// 11) LT exh rate
       /*25*/ par.push_back(log(LT_exh_rate_));

       /// 12) apoptosis related parameters
       /*26*/ par.push_back(log(t_apop_meas_));
       /*27*/ par.push_back(log(t_duration_apoptosis_));


        /// Media
        /*1*/ par.push_back(log(TNF_deg_));
        /*2*/ par.push_back(log(IFN_deg_));
        /*3*/ par.push_back(log(TymidineTriteate_));
        /*4*/ par.push_back(log(Prol_TymTr_));
    }

    else if (mode_=="PARTIAL")
    {   // multiplicar si estoy seguro del valor


   /// APC
        /// 1) Init ratio of cells
        /*1*/   par.push_back(log(init_ratio_APC_));
        /// 2) IFN Poductions rates of each type of APC
        /*2*/ par.push_back(log(IFN_APC0_prod_rate_));
        /*3*/ par.push_back(log(IFN_APCa_prod_rate_));
        /*4*/ par.push_back(log(IFN_APCbo_prod_rate_));


        /// 3) TNF Poductions rates of each type of APC
        /*5*/ par.push_back(log(TNF_APC0_prod_rate_));
        /*6*/ par.push_back(log(TNF_APCa_prod_rate_));
        /*7*/ par.push_back(log(TNF_APCbo_prod_rate_));


        /// 4) Percentages of IFN productions of each type of APC
        /*8*/ par.push_back(log(percentage_IFN_APC0_prod_rate_));
        /*9*/ par.push_back(log(percentage_IFN_APCa_prod_rate_));
        /*10*/ par.push_back(log(percentage_IFN_APCbo_prod_rate_));

        /// 5)Percentages of TNF productions of each type of APC
        /*11*/ par.push_back(log(percentage_TNF_APC0_prod_rate_));
        /*12*/ par.push_back(log(percentage_TNF_APCa_prod_rate_));
        /*13*/ par.push_back(log(percentage_TNF_APCbo_prod_rate_));


        /// 6) Proliferation rates
        /*14*/ par.push_back(log(APC_bound_proliferation_rate_));

        /// 7) Apoptosis rates
        /*15*/ par.push_back(log(APC0_apop_rate_));
        /*16*/ par.push_back(log(APCa_apop_rate_));
        /*17*/ par.push_back(log(APCbo_apop_rate_));
        /*18*/ par.push_back(log(APCbl_apop_rate_));
        /*19*/ par.push_back(log(APCexh_apop_rate_));

        /// 8) constant saturation of TNF for apoptosis
        /*20*/ par.push_back(log(Ks_APC_m_TNF_));

        /// 9) conversion rates
        /*21*/ par.push_back(log(APC_Ag_));
        /*22*/ par.push_back(log(APC_APC_));
        /*23*/ par.push_back(log(APC_NK_));
        /*24*/ par.push_back(log(APC_LT_1_));
        /*25*/ par.push_back(log(APC_LT_2_));
        /*26*/ par.push_back(log(APC_Ab_));
        /*27*/ par.push_back(log(APC_exh_));

        /// 10)Saturation constant of IFN and TNF for activation
        /*28*/ par.push_back(log(KsAPC_LT_));

        /// 11)Saturation constant of APC_LT interaction
        /*29*/ par.push_back(log(APC_Ksi_));
        /*30*/ par.push_back(log(APC_Kst_));

        /// 12) Percentages of cell expressing receptor
        /*31*/ par.push_back(log(APC0_expressing_receptor_));
        /*32*/ par.push_back(log(APCa_expressing_receptor_));
        /// 13) Apoptosis rate for TNF
        /*33*/ par.push_back(log(u_APC_TNF_));

   /// NK
       /// 1) Init ratio of cells
       /*1*/  par.push_back(log(init_ratio_NK_));

       /// 2) IFN Poductions rates of each type of NK
       /*2*/  par.push_back(log(IFN_NK0_prod_rate_));
       /*3*/  par.push_back(log(IFN_NKa_prod_rate_));
       /*4*/  par.push_back(log(IFN_NKbo_prod_rate_));

       /// 3) TNF Poductions rates of each type of NK
       /*5*/  par.push_back(log(TNF_NK0_prod_rate_));
       /*6*/  par.push_back(log(TNF_NKa_prod_rate_));
       /*7*/  par.push_back(log(TNF_NKbo_prod_rate_));

       /// 4) Percentages of IFN productions of each type of NK
       /*8*/  par.push_back(log(percentage_IFN_NK0_prod_rate_));
       /*9*/  par.push_back(log(percentage_IFN_AgNKa_prod_rate_));
       /*10*/  par.push_back(log(percentage_IFN_NKbo_prod_rate_));

       /// 5)Percentages of TNF productions of each type of NK
       /*11*/  par.push_back(log(percentage_TNF_NK0_prod_rate_));
       /*12*/  par.push_back(log(percentage_TNF_NKa_prod_rate_));
       /*13*/  par.push_back(log(percentage_TNF_NKbo_prod_rate_));

       /// 6) Proliferation rates
       /*13.5*/par.push_back (log(NK0_proliferation_rate_));
       /*14*/  par.push_back(log(NKa_proliferation_rate_));
       /*15*/  par.push_back(log(NKbo_proliferation_rate_));
       /*16*/  par.push_back(log(NKbl_proliferation_rate_));

       /// 7) Apoptosis rates
       /*17*/  par.push_back(log(NK0_apop_rate_));
       /*18*/  par.push_back(log(NKa_apop_rate_));
       /*19*/  par.push_back(log(NKbo_apop_rate_));
       /*20*/  par.push_back(log(NKbl_apop_rate_));
       /*21*/  par.push_back(log(NKexh_apop_rate_));

       /// 8) constant saturation of TNF for apoptosis
       /*22*/  par.push_back(log(Ks_NK_m_TNF_));

       /// 9) conversion rates
       /*23*/  par.push_back(log(KaNK_));
       /*24*/  par.push_back(log(NK_NK_));
       /*25*/  par.push_back(log(NK_Ab_));
       /*26*/  par.push_back(log(NK_exh_));

       /// 10)Saturation constant of NK interaction for activation
       /*27*/  par.push_back(log(KsAPC_NK_));

       /// 11)Saturation constant of NK_LT interaction
       /*28*/  par.push_back(log(NK_Ksi_));
       /*29*/  par.push_back(log(NK_Kst_));

       /// 12) Percentages of cell expressing receptor
       /*30*/  par.push_back(log(NK0_expressing_receptor_));
       /*31*/  par.push_back(log(NKa_expressing_receptor_));

       /// 13) Apoptosis rate for TNF
       /*32*/  par.push_back(log(u_NK_TNF_));

  /// LT
      /// 1) Init number of LT
      /*1*/  par.push_back(log(ratio_init_LTns_));
      /*2*/  par.push_back(log(ratio_initLTspecific_));

      /// 2) IFN Poductions rates of each type of LT
      /*3*/  par.push_back(log(IFN_LTns_prod_rate_));
      /*4*/  par.push_back(log(IFN_LTbo_prod_rate_));
      /*5*/  par.push_back(log(IFN_LTbl_prod_rate_));

      /// 3) TNF Poductions rates of each type of LT
      /*6*/  par.push_back(log(TNF_LTns_prod_rate_));
      /*7*/  par.push_back(log(TNF_LTbo_prod_rate_));
      /*8*/  par.push_back(log(TNF_LTbl_prod_rate_));

      /// 4) Percentages of IFN productions of each type of LT
      /*9*/  par.push_back(log(percentage_IFN_LTns_prod_rate_));
      /*10*/  par.push_back(log(percentage_IFN_LTbo_prod_rate_));
      /*11*/  par.push_back(log(percentage_IFN_LTbl_prod_rate_));

      /// 5)Percentages of TNF productions of each type of LT
      /*12*/  par.push_back(log(percentage_TNF_LTns_prod_rate_));
      /*13*/  par.push_back(log(percentage_TNF_LTbo_prod_rate_));
      /*14*/  par.push_back(log(percentage_TNF_LTbl_prod_rate_));

      /// 6) Proliferation rates
      /*15*/  par.push_back(log(LTns_proliferation_rate_));
      /*16*/  par.push_back(log(LTbo_proliferation_rate_));
      /*17*/  par.push_back(log(LTbl_proliferation_rate_));

      /// 7) Apoptosis rates
      /*18*/  par.push_back(log(LTns_apop_rate_));
      /*19*/  par.push_back(log(LTbo_apop_rate_));
      /*20*/  par.push_back(log(LTbl_apop_rate_));
      /*21*/  par.push_back(log(LTexh_apop_rate_));

      /// 8) constant saturation of TNF for apoptosis
      /*22*/  par.push_back(log(Ks_LT_m_TNF_));

      /// 9) Percentages of cell expressing receptor
      /*23*/  par.push_back(log(LTns_expressing_receptor_));

      /// 10) Apoptosis rate for TNF
      /*24*/  par.push_back(log(u_LT_TNF_));

      /// 11) LT exh rate
      /*25*/ par.push_back(log(LT_exh_rate_));

      /// 12) apoptosis related parameters
      /*26*/ par.push_back(log(t_apop_meas_));
      /*27*/ par.push_back(log(t_duration_apoptosis_));



        /// Media
        /*1*/ par.push_back(log(TNF_deg_));
        /*3*/ par.push_back(log(TymidineTriteate_));
        /*4*/ par.push_back(log(Prol_TymTr_));


    }

    return par;


}

SimParameters& SimParameters::applyParameters(const std::vector<double>& param)
{
    if (mode_=="FULL")
    {
        std::size_t i=0;
        // Dividir si estoy seguro el valor
/// APC
        /// 1) Init ratio of cells
        /*1*/  init_ratio_APC_ = exp(param[i++]);
        /// 2) IFN Poductions rates of each type of APC
        /*2*/  IFN_APC0_prod_rate_ = exp(param[i++]);
        /*3*/  IFN_APCa_prod_rate_ = exp(param[i++]);
        /*4*/  IFN_APCbo_prod_rate_ = exp(param[i++]);


        /// 3) TNF Poductions rates of each type of APC
        /*5*/  TNF_APC0_prod_rate_ = exp(param[i++]);
        /*6*/  TNF_APCa_prod_rate_ = exp(param[i++]);
        /*7*/  TNF_APCbo_prod_rate_ = exp(param[i++]);


        /// 4) Percentages of IFN productions of each type of APC
        /*8*/  percentage_IFN_APC0_prod_rate_ = exp(param[i++]);
        /*9*/  percentage_IFN_APCa_prod_rate_ = exp(param[i++]);
        /*10*/  percentage_IFN_APCbo_prod_rate_ = exp(param[i++]);

        /// 5)Percentages of TNF productions of each type of APC
        /*11*/  percentage_TNF_APC0_prod_rate_ = exp(param[i++]);
        /*12*/  percentage_TNF_APCa_prod_rate_ = exp(param[i++]);
        /*13*/  percentage_TNF_APCbo_prod_rate_ = exp(param[i++]);


        /// 6) Proliferation rates
        /*14*/  APC_bound_proliferation_rate_ = exp(param[i++]);

        /// 7) Apoptosis rates
        /*15*/  APC0_apop_rate_ = exp(param[i++]);
        /*16*/  APCa_apop_rate_ = exp(param[i++]);
        /*17*/  APCbo_apop_rate_ = exp(param[i++]);
        /*18*/  APCbl_apop_rate_ = exp(param[i++]);
        /*19*/  APCexh_apop_rate_ = exp(param[i++]);

        /// 8) constant saturation of TNF for apoptosis
        /*20*/  Ks_APC_m_TNF_ = exp(param[i++]);

        /// 9) conversion rates
        /*21*/  APC_Ag_ = exp(param[i++]);
        /*22*/  APC_APC_ = exp(param[i++]);
        /*23*/  APC_NK_ = exp(param[i++]);
        /*24*/  APC_LT_1_ = exp(param[i++]);
        /*25*/  APC_LT_2_ = exp(param[i++]);
        /*26*/  APC_Ab_ = exp(param[i++]);
        /*27*/  APC_exh_ = exp(param[i++]);

        /// 10)Saturation constant of IFN and TNF for activation
        /*28*/  KsAPC_LT_ = exp(param[i++]);

        /// 11)Saturation constant of APC_LT interaction
        /*29*/  APC_Ksi_ = exp(param[i++]);
        /*30*/  APC_Kst_ = exp(param[i++]);

        /// 12) Percentages of cell expressing receptor
        /*31*/  APC0_expressing_receptor_ = exp(param[i++]);
        /*32*/  APCa_expressing_receptor_ = exp(param[i++]);
        /// 13) Apoptosis rate for TNF
        /*33*/ u_APC_TNF_ = exp(param[i++]);

/// NK
        /// 1) Init ratio of cells
        /*1*/ init_ratio_NK_= exp(param[i++]);

        /// 2) IFN Poductions rates of each type of NK
        /*2*/ IFN_NK0_prod_rate_= exp(param[i++]);
        /*3*/ IFN_NKa_prod_rate_= exp(param[i++]);
        /*4*/ IFN_NKbo_prod_rate_= exp(param[i++]);

        /// 3) TNF Poductions rates of each type of NK
        /*5*/ TNF_NK0_prod_rate_= exp(param[i++]);
        /*6*/ TNF_NKa_prod_rate_= exp(param[i++]);
        /*7*/ TNF_NKbo_prod_rate_= exp(param[i++]);

        /// 4) Percentages of IFN productions of each type of NK
        /*8*/ percentage_IFN_NK0_prod_rate_= exp(param[i++]);
        /*9*/ percentage_IFN_AgNKa_prod_rate_= exp(param[i++]);
        /*10*/ percentage_IFN_NKbo_prod_rate_= exp(param[i++]);

        /// 5)Percentages of TNF productions of each type of NK
        /*11*/ percentage_TNF_NK0_prod_rate_= exp(param[i++]);
        /*12*/ percentage_TNF_NKa_prod_rate_= exp(param[i++]);
        /*13*/ percentage_TNF_NKbo_prod_rate_= exp(param[i++]);

        /// 6) Proliferation rates
        /*13.5*/NK0_proliferation_rate_= exp(param[i++]);
        /*14*/ NKa_proliferation_rate_= exp(param[i++]);
        /*15*/ NKbo_proliferation_rate_= exp(param[i++]);
        /*16*/ NKbl_proliferation_rate_= exp(param[i++]);

        /// 7) Apoptosis rates
        /*17*/ NK0_apop_rate_= exp(param[i++]);
        /*18*/ NKa_apop_rate_= exp(param[i++]);
        /*19*/ NKbo_apop_rate_= exp(param[i++]);
        /*20*/ NKbl_apop_rate_= exp(param[i++]);
        /*21*/ NKexh_apop_rate_= exp(param[i++]);

        /// 8) constant saturation of TNF for apoptosis
        /*22*/ Ks_NK_m_TNF_= exp(param[i++]);

        /// 9) conversion rates
        /*23*/ KaNK_= exp(param[i++]);
        /*24*/ NK_NK_= exp(param[i++]);
        /*25*/ NK_Ab_= exp(param[i++]);
        /*26*/ NK_exh_= exp(param[i++]);

        /// 10)Saturation constant of NK interaction for activation
        /*27*/ KsAPC_NK_= exp(param[i++]);

        /// 11)Saturation constant of NK_LT interaction
        /*28*/ NK_Ksi_= exp(param[i++]);
        /*29*/ NK_Kst_= exp(param[i++]);

        /// 12) Percentages of cell expressing receptor
        /*30*/ NK0_expressing_receptor_= exp(param[i++]);
        /*31*/ NKa_expressing_receptor_= exp(param[i++]);

        /// 13) Apoptosis rate for TNF
        /*32*/ u_NK_TNF_= exp(param[i++]);

/// LT
        /// 1) Init number of LT
        /*1*/  ratio_init_LTns_=exp(param[i++]);
        /*2*/  ratio_initLTspecific_=exp(param[i++]);

       /// 2) IFN Poductions rates of each type of LT
       /*3*/  IFN_LTns_prod_rate_=exp(param[i++]);
       /*4*/  IFN_LTbo_prod_rate_=exp(param[i++]);
       /*5*/  IFN_LTbl_prod_rate_=exp(param[i++]);

       /// 3) TNF Poductions rates of each type of LT
       /*6*/  TNF_LTns_prod_rate_=exp(param[i++]);
       /*7*/  TNF_LTbo_prod_rate_=exp(param[i++]);
       /*8*/  TNF_LTbl_prod_rate_=exp(param[i++]);

       /// 4) Percentages of IFN productions of each type of LT
       /*9*/  percentage_IFN_LTns_prod_rate_=exp(param[i++]);
       /*10*/  percentage_IFN_LTbo_prod_rate_=exp(param[i++]);
       /*11*/  percentage_IFN_LTbl_prod_rate_=exp(param[i++]);


       /// 5)Percentages of TNF productions of each type of LT
       /*12*/  percentage_TNF_LTns_prod_rate_=exp(param[i++]);
       /*13*/  percentage_TNF_LTbo_prod_rate_=exp(param[i++]);
       /*14*/  percentage_TNF_LTbl_prod_rate_=exp(param[i++]);

       /// 6) Proliferation rates
       /*15*/  LTns_proliferation_rate_=exp(param[i++]);
       /*16*/  LTbo_proliferation_rate_=exp(param[i++]);
       /*17*/  LTbl_proliferation_rate_=exp(param[i++]);

       /// 7) Apoptosis rates
       /*18*/  LTns_apop_rate_=exp(param[i++]);
       /*19*/  LTbo_apop_rate_=exp(param[i++]);
       /*20*/  LTbl_apop_rate_=exp(param[i++]);
       /*21*/  LTexh_apop_rate_=exp(param[i++]);

       /// 8) constant saturation of TNF for apoptosis
       /*22*/  Ks_LT_m_TNF_=exp(param[i++]);

       /// 9) Percentages of cell expressing receptor
       /*23*/  LTns_expressing_receptor_=exp(param[i++]);

       /// 10) Apoptosis rate for TNF
       /*24*/  u_LT_TNF_=exp(param[i++]);

       /// 11) LT exh rate
       /*25*/ LT_exh_rate_=exp(param[i++]);


/// Media
        /*1*/ TNF_deg_=exp(param[i++]);
        /*2*/ IFN_deg_=exp(param[i++]);
        /*3*/ TymidineTriteate_=exp(param[i++]);
        /*4*/ Prol_TymTr_=exp(param[i++]);


    }
    else if(mode_=="PARTIAL")
    {
        std::size_t i=0;
        // Dividir si estoy seguro el valor
        /// 1) Init ratio of cells
        /*1*/  init_ratio_APC_ = exp(param[i++]);
        /// 2) IFN Poductions rates of each type of APC
        /*2*/  IFN_APC0_prod_rate_ = exp(param[i++]);
        /*3*/  IFN_APCa_prod_rate_ = exp(param[i++]);
        /*4*/  IFN_APCbo_prod_rate_ = exp(param[i++]);


        /// 3) TNF Poductions rates of each type of APC
        /*5*/  TNF_APC0_prod_rate_ = exp(param[i++]);
        /*6*/  TNF_APCa_prod_rate_ = exp(param[i++]);
        /*7*/  TNF_APCbo_prod_rate_ = exp(param[i++]);


        /// 4) Percentages of IFN productions of each type of APC
        /*8*/  percentage_IFN_APC0_prod_rate_ = exp(param[i++]);
        /*9*/  percentage_IFN_APCa_prod_rate_ = exp(param[i++]);
        /*10*/  percentage_IFN_APCbo_prod_rate_ = exp(param[i++]);

        /// 5)Percentages of TNF productions of each type of APC
        /*11*/  percentage_TNF_APC0_prod_rate_ = exp(param[i++]);
        /*12*/  percentage_TNF_APCa_prod_rate_ = exp(param[i++]);
        /*13*/  percentage_TNF_APCbo_prod_rate_ = exp(param[i++]);


        /// 6) Proliferation rates
        /*14*/  APC_bound_proliferation_rate_ = exp(param[i++]);

        /// 7) Apoptosis rates
        /*15*/  APC0_apop_rate_ = exp(param[i++]);
        /*16*/  APCa_apop_rate_ = exp(param[i++]);
        /*17*/  APCbo_apop_rate_ = exp(param[i++]);
        /*18*/  APCbl_apop_rate_ = exp(param[i++]);
        /*19*/  APCexh_apop_rate_ = exp(param[i++]);

        /// 8) constant saturation of TNF for apoptosis
        /*20*/  Ks_APC_m_TNF_ = exp(param[i++]);

        /// 9) conversion rates
        /*21*/  APC_Ag_ = exp(param[i++]);
        /*22*/  APC_APC_ = exp(param[i++]);
        /*23*/  APC_NK_ = exp(param[i++]);
        /*24*/  APC_LT_1_ = exp(param[i++]);
        /*25*/  APC_LT_2_ = exp(param[i++]);
        /*26*/  APC_Ab_ = exp(param[i++]);
        /*27*/  APC_exh_ = exp(param[i++]);

        /// 10)Saturation constant of IFN and TNF for activation
        /*28*/  KsAPC_LT_ = exp(param[i++]);

        /// 11)Saturation constant of APC_LT interaction
        /*29*/  APC_Ksi_ = exp(param[i++]);
        /*30*/  APC_Kst_ = exp(param[i++]);

        /// 12) Percentages of cell expressing receptor
        /*31*/  APC0_expressing_receptor_ = exp(param[i++]);
        /*32*/  APCa_expressing_receptor_ = exp(param[i++]);
        /// 13) Apoptosis rate for TNF
        /*33*/ u_APC_TNF_ = exp(param[i++]);

/// NK
        /// 1) Init ratio of cells
        /*1*/ init_ratio_NK_= exp(param[i++]);

        /// 2) IFN Poductions rates of each type of NK
        /*2*/ IFN_NK0_prod_rate_= exp(param[i++]);
        /*3*/ IFN_NKa_prod_rate_= exp(param[i++]);
        /*4*/ IFN_NKbo_prod_rate_= exp(param[i++]);

        /// 3) TNF Poductions rates of each type of NK
        /*5*/ TNF_NK0_prod_rate_= exp(param[i++]);
        /*6*/ TNF_NKa_prod_rate_= exp(param[i++]);
        /*7*/ TNF_NKbo_prod_rate_= exp(param[i++]);

        /// 4) Percentages of IFN productions of each type of NK
        /*8*/ percentage_IFN_NK0_prod_rate_= exp(param[i++]);
        /*9*/ percentage_IFN_AgNKa_prod_rate_= exp(param[i++]);
        /*10*/ percentage_IFN_NKbo_prod_rate_= exp(param[i++]);

        /// 5)Percentages of TNF productions of each type of NK
        /*11*/ percentage_TNF_NK0_prod_rate_= exp(param[i++]);
        /*12*/ percentage_TNF_NKa_prod_rate_= exp(param[i++]);
        /*13*/ percentage_TNF_NKbo_prod_rate_= exp(param[i++]);

        /// 6) Proliferation rates
        /*13.5*/ NK0_proliferation_rate_= exp(param[i++]);
        /*14*/ NKa_proliferation_rate_= exp(param[i++]);
        /*15*/ NKbo_proliferation_rate_= exp(param[i++]);
        /*16*/ NKbl_proliferation_rate_= exp(param[i++]);

        /// 7) Apoptosis rates
        /*17*/ NK0_apop_rate_= exp(param[i++]);
        /*18*/ NKa_apop_rate_= exp(param[i++]);
        /*19*/ NKbo_apop_rate_= exp(param[i++]);
        /*20*/ NKbl_apop_rate_= exp(param[i++]);
        /*21*/ NKexh_apop_rate_= exp(param[i++]);

        /// 8) constant saturation of TNF for apoptosis
        /*22*/ Ks_NK_m_TNF_= exp(param[i++]);

        /// 9) conversion rates
        /*23*/ KaNK_= exp(param[i++]);
        /*24*/ NK_NK_= exp(param[i++]);
        /*25*/ NK_Ab_= exp(param[i++]);
        /*26*/ NK_exh_= exp(param[i++]);

        /// 10)Saturation constant of NK interaction for activation
        /*27*/ KsAPC_NK_= exp(param[i++]);

        /// 11)Saturation constant of NK_LT interaction
        /*28*/ NK_Ksi_= exp(param[i++]);
        /*29*/ NK_Kst_= exp(param[i++]);

        /// 12) Percentages of cell expressing receptor
        /*30*/ NK0_expressing_receptor_= exp(param[i++]);
        /*31*/ NKa_expressing_receptor_= exp(param[i++]);

        /// 13) Apoptosis rate for TNF
        /*32*/ u_NK_TNF_= exp(param[i++]);

 /// LT
        /// 1) Init number of LT
        /*1*/  ratio_init_LTns_=exp(param[i++]);
        /*2*/  ratio_initLTspecific_=exp(param[i++]);

        /// 2) IFN Poductions rates of each type of LT
        /*3*/  IFN_LTns_prod_rate_=exp(param[i++]);
        /*4*/  IFN_LTbo_prod_rate_=exp(param[i++]);
        /*5*/  IFN_LTbl_prod_rate_=exp(param[i++]);

        /// 3) TNF Poductions rates of each type of LT
        /*6*/  TNF_LTns_prod_rate_=exp(param[i++]);
        /*7*/  TNF_LTbo_prod_rate_=exp(param[i++]);
        /*8*/  TNF_LTbl_prod_rate_=exp(param[i++]);

        /// 4) Percentages of IFN productions of each type of LT
        /*9*/  percentage_IFN_LTns_prod_rate_=exp(param[i++]);
        /*10*/  percentage_IFN_LTbo_prod_rate_=exp(param[i++]);
        /*11*/  percentage_IFN_LTbl_prod_rate_=exp(param[i++]);


        /// 5)Percentages of TNF productions of each type of LT
        /*12*/  percentage_TNF_LTns_prod_rate_=exp(param[i++]);
        /*13*/  percentage_TNF_LTbo_prod_rate_=exp(param[i++]);
        /*14*/  percentage_TNF_LTbl_prod_rate_=exp(param[i++]);

        /// 6) Proliferation rates
        /*15*/  LTns_proliferation_rate_=exp(param[i++]);
        /*16*/  LTbo_proliferation_rate_=exp(param[i++]);
        /*17*/  LTbl_proliferation_rate_=exp(param[i++]);

        /// 7) Apoptosis rates
        /*18*/  LTns_apop_rate_=exp(param[i++]);
        /*19*/  LTbo_apop_rate_=exp(param[i++]);
        /*20*/  LTbl_apop_rate_=exp(param[i++]);
        /*21*/  LTexh_apop_rate_=exp(param[i++]);

        /// 8) constant saturation of TNF for apoptosis
        /*22*/  Ks_LT_m_TNF_=exp(param[i++]);

        /// 9) Percentages of cell expressing receptor
        /*23*/  LTns_expressing_receptor_=exp(param[i++]);

        /// 10) Apoptosis rate for TNF
        /*24*/  u_LT_TNF_=exp(param[i++]);

        /// 11) LT exh rate
        /*25*/ LT_exh_rate_=exp(param[i++]);


/// Media
       /*1*/ TNF_deg_=exp(param[i++]);
       /*2*/ IFN_deg_= TNF_deg_;
       /*3*/ TymidineTriteate_=exp(param[i++]);
       /*4*/ Prol_TymTr_=exp(param[i++]);

    }
    return *this;
}




SimParameters::SimParameters(const SimParameters& other):
    mode_(other.mode_),
    /// 1) Init ratio of cells
    /*1*/  init_ratio_APC_(other.init_ratio_APC_),

    /// 2) IFN Poductions rates of each type of APC
    /*2*/  IFN_APC0_prod_rate_ (other.IFN_APC0_prod_rate_),
    /*3*/  IFN_APCa_prod_rate_ (other.IFN_APCa_prod_rate_),
    /*4*/  IFN_APCbo_prod_rate_ (other.IFN_APCbo_prod_rate_),


    /// 3) TNF Poductions rates of each type of APC
    /*5*/  TNF_APC0_prod_rate_ (other.TNF_APC0_prod_rate_),
    /*6*/  TNF_APCa_prod_rate_ (other.TNF_APCa_prod_rate_),
    /*7*/  TNF_APCbo_prod_rate_ (other.TNF_APCbo_prod_rate_),

    /// 4) Percentages of IFN productions of each type of APC
    /*8*/  percentage_IFN_APC0_prod_rate_ (other.percentage_IFN_APC0_prod_rate_),
    /*9*/  percentage_IFN_APCa_prod_rate_ (other.percentage_IFN_APCa_prod_rate_),
    /*10*/  percentage_IFN_APCbo_prod_rate_ (other.percentage_IFN_APCbo_prod_rate_),

    /// 5)Percentages of TNF productions of each type of APC
    /*11*/  percentage_TNF_APC0_prod_rate_ (other.percentage_TNF_APC0_prod_rate_),
    /*12*/  percentage_TNF_APCa_prod_rate_ (other.percentage_TNF_APCa_prod_rate_),
    /*13*/  percentage_TNF_APCbo_prod_rate_ (other.percentage_TNF_APCbo_prod_rate_),

    /// 6) Proliferation rates
    /*14*/  APC_bound_proliferation_rate_ (other.APC_bound_proliferation_rate_),

    /// 7) Apoptosis rates
    /*15*/  APC0_apop_rate_ (other.APC0_apop_rate_),
    /*16*/  APCa_apop_rate_ (other.APCa_apop_rate_),
    /*17*/  APCbo_apop_rate_ (other.APCbo_apop_rate_),
    /*18*/  APCbl_apop_rate_ (other.APCbl_apop_rate_),
    /*19*/  APCexh_apop_rate_ (other.APCexh_apop_rate_),

    /// 8) constant saturation of TNF for apoptosis
    /*20*/  Ks_APC_m_TNF_ (other.Ks_APC_m_TNF_),

    /// 9) conversion rates
    /*21*/  APC_Ag_ (other.APC_Ag_),
    /*22*/  APC_APC_ (other.APC_APC_),
    /*23*/  APC_NK_ (other.APC_NK_),
    /*24*/  APC_LT_1_ (other.APC_LT_1_),
    /*25*/  APC_LT_2_ (other.APC_LT_2_),
    /*26*/  APC_Ab_ (other.APC_Ab_),
    /*27*/  APC_exh_ (other.APC_exh_),

    /// 10)Saturation constant of IFN and TNF for activation
    /*28*/  KsAPC_LT_ (other.KsAPC_LT_),

    /// 11)Saturation constant of APC_LT interaction
    /*29*/  APC_Ksi_ (other.APC_Ksi_),
    /*30*/  APC_Kst_ (other.APC_Kst_),

    /// 12) Percentages of cell expressing receptor
    /*31*/  APC0_expressing_receptor_ (other.APC0_expressing_receptor_),
    /*32*/  APCa_expressing_receptor_ (other.APCa_expressing_receptor_),
    /// 13) Apoptosis rate for TNF
    /*33*/ u_APC_TNF_ (other.u_APC_TNF_),

/// NK
    /// 1) Init ratio of cells
    /*1*/ init_ratio_NK_ (other.init_ratio_NK_),

    /// 2) IFN Poductions rates of each type of NK
    /*2*/ IFN_NK0_prod_rate_(other.IFN_NK0_prod_rate_),
    /*3*/ IFN_NKa_prod_rate_(other.IFN_NKa_prod_rate_),
    /*4*/ IFN_NKbo_prod_rate_(other.IFN_NKbo_prod_rate_),

    /// 3) TNF Poductions rates of each type of NK
    /*5*/ TNF_NK0_prod_rate_(other.TNF_NK0_prod_rate_),
    /*6*/ TNF_NKa_prod_rate_(other.TNF_NKa_prod_rate_),
    /*7*/ TNF_NKbo_prod_rate_(other.TNF_NKbo_prod_rate_),

    /// 4) Percentages of IFN productions of each type of NK
    /*8*/ percentage_IFN_NK0_prod_rate_(other.percentage_IFN_NK0_prod_rate_),
    /*9*/ percentage_IFN_AgNKa_prod_rate_(other.percentage_IFN_AgNKa_prod_rate_),
    /*10*/ percentage_IFN_NKbo_prod_rate_(other.percentage_IFN_NKbo_prod_rate_),

    /// 5)Percentages of TNF productions of each type of NK
    /*11*/ percentage_TNF_NK0_prod_rate_(other.percentage_TNF_NK0_prod_rate_),
    /*12*/ percentage_TNF_NKa_prod_rate_(other.percentage_TNF_NKa_prod_rate_),
    /*13*/ percentage_TNF_NKbo_prod_rate_(other.percentage_TNF_NKbo_prod_rate_),

    /// 6) Proliferation rates
    /*13.5*/ NK0_proliferation_rate_ (other. NK0_proliferation_rate_),
    /*14*/ NKa_proliferation_rate_(other.NKa_proliferation_rate_),
    /*15*/ NKbo_proliferation_rate_(other.NKbo_proliferation_rate_),
    /*16*/ NKbl_proliferation_rate_(other.NKbl_proliferation_rate_),

    /// 7) Apoptosis rates
    /*17*/ NK0_apop_rate_(other.NK0_apop_rate_),
    /*18*/ NKa_apop_rate_(other.NKa_apop_rate_),
    /*19*/ NKbo_apop_rate_(other.NKbo_apop_rate_),
    /*20*/ NKbl_apop_rate_(other.NKbl_apop_rate_),
    /*21*/ NKexh_apop_rate_(other.NKexh_apop_rate_),

    /// 8) constant saturation of TNF for apoptosis
    /*22*/ Ks_NK_m_TNF_(other.Ks_NK_m_TNF_),

    /// 9) conversion rates
    /*23*/ KaNK_(other.KaNK_),
    /*24*/ NK_NK_(other.NK_NK_),
    /*25*/ NK_Ab_(other.NK_Ab_),
    /*26*/ NK_exh_(other.NK_exh_),

    /// 10)Saturation constant of NK interaction for activation
    /*27*/ KsAPC_NK_(other.KsAPC_NK_),

    /// 11)Saturation constant of NK_LT interaction
    /*28*/ NK_Ksi_(other.NK_Ksi_),
    /*29*/ NK_Kst_(other.NK_Kst_),

    /// 12) Percentages of cell expressing receptor
    /*30*/ NK0_expressing_receptor_(other.NK0_expressing_receptor_),
    /*31*/ NKa_expressing_receptor_(other.NKa_expressing_receptor_),

    /// 13) Apoptosis rate for TNF
    /*32*/ u_NK_TNF_(other.u_NK_TNF_),


/// LT
    /// 1) Init number of LT
    /*1*/  ratio_init_LTns_(other.ratio_init_LTns_),
    /*2*/  ratio_initLTspecific_(other.ratio_initLTspecific_),

    /// 2) IFN Poductions rates of each type of LT
    /*3*/  IFN_LTns_prod_rate_(other.ratio_initLTspecific_),
    /*4*/  IFN_LTbo_prod_rate_(other.IFN_LTbo_prod_rate_),
    /*5*/  IFN_LTbl_prod_rate_(other.IFN_LTbl_prod_rate_),

    /// 3) TNF Poductions rates of each type of LT
    /*6*/  TNF_LTns_prod_rate_(other.TNF_LTns_prod_rate_),
    /*7*/  TNF_LTbo_prod_rate_(other.TNF_LTbo_prod_rate_),
    /*8*/  TNF_LTbl_prod_rate_(other.TNF_LTbl_prod_rate_),

    /// 4) Percentages of IFN productions of each type of LT
    /*9*/  percentage_IFN_LTns_prod_rate_(other.percentage_IFN_LTns_prod_rate_),
    /*10*/  percentage_IFN_LTbo_prod_rate_(other.percentage_IFN_LTbo_prod_rate_),
    /*11*/  percentage_IFN_LTbl_prod_rate_(other.percentage_IFN_LTbl_prod_rate_),


    /// 5)Percentages of TNF productions of each type of LT
    /*12*/  percentage_TNF_LTns_prod_rate_(other.percentage_TNF_LTns_prod_rate_),
    /*13*/  percentage_TNF_LTbo_prod_rate_(other.percentage_TNF_LTbo_prod_rate_),
    /*14*/  percentage_TNF_LTbl_prod_rate_(other.percentage_TNF_LTbl_prod_rate_),

    /// 6) Proliferation rates
    /*15*/  LTns_proliferation_rate_(other.LTns_proliferation_rate_),
    /*16*/  LTbo_proliferation_rate_(other.LTbo_proliferation_rate_),
    /*17*/  LTbl_proliferation_rate_(other.LTbl_proliferation_rate_),

    /// 7) Apoptosis rates
    /*18*/  LTns_apop_rate_(other.LTns_apop_rate_),
    /*19*/  LTbo_apop_rate_(other.LTbo_apop_rate_),
    /*20*/  LTbl_apop_rate_(other.LTbl_apop_rate_),
    /*21*/  LTexh_apop_rate_(other.LTexh_apop_rate_),

    /// 8) constant saturation of TNF for apoptosis
    /*22*/  Ks_LT_m_TNF_(other.LTexh_apop_rate_),

    /// 9) Percentages of cell expressing receptor
    /*23*/  LTns_expressing_receptor_(other.LTexh_apop_rate_),

    /// 10) Apoptosis rate for TNF
    /*24*/  u_LT_TNF_(other.LTexh_apop_rate_),

    /// 11) LT exh rate
    /*25*/ LT_exh_rate_(other.LT_exh_rate_),

    /// 12) apoptosis related parameters
    /*26*/ t_apop_meas_(other.t_apop_meas_),
    /*27*/ t_duration_apoptosis_(other.t_duration_apoptosis_),

/// Media
    /*1*/ TNF_deg_(other.TNF_deg_),
    /*2*/ IFN_deg_(other.IFN_deg_),
    /*3*/ TymidineTriteate_(other.TymidineTriteate_),
    /*4*/ Prol_TymTr_(other.Prol_TymTr_)

{}

SimParameters&
SimParameters::operator=(const SimParameters& other)
{
    if (this!=&other)
    {
        SimParameters tmp(other);
        swap(*this,tmp);
    }
    return *this;
}

void swap(SimParameters& one, SimParameters& other)
{
    std::swap(one.mode_,other.mode_);
    /// 1) Init ratio of cells
    /*1*/  std::swap(one.init_ratio_APC_,other.init_ratio_APC_);
    /// 2) IFN Poductions rates of each type of APC
    /*2*/  std::swap(one.IFN_APC0_prod_rate_ ,other.IFN_APC0_prod_rate_);
    /*3*/  std::swap(one.IFN_APCa_prod_rate_ ,other.IFN_APCa_prod_rate_);
    /*4*/  std::swap(one.IFN_APCbo_prod_rate_ ,other.IFN_APCbo_prod_rate_);


    /// 3) TNF Poductions rates of each type of APC
    /*5*/  std::swap(one.TNF_APC0_prod_rate_ ,other.TNF_APC0_prod_rate_);
    /*6*/  std::swap(one.TNF_APCa_prod_rate_ ,other.TNF_APCa_prod_rate_);
    /*7*/  std::swap(one.TNF_APCbo_prod_rate_ ,other.TNF_APCbo_prod_rate_);


    /// 4) Percentages of IFN productions of each type of APC
    /*8*/  std::swap(one.percentage_IFN_APC0_prod_rate_ ,other.percentage_IFN_APC0_prod_rate_);
    /*9*/  std::swap(one.percentage_IFN_APCa_prod_rate_ ,other.percentage_IFN_APCa_prod_rate_);
    /*10*/  std::swap(one.percentage_IFN_APCbo_prod_rate_ ,other.percentage_IFN_APCbo_prod_rate_);

    /// 5)Percentages of TNF productions of each type of APC
    /*11*/  std::swap(one.percentage_TNF_APC0_prod_rate_ ,other.percentage_TNF_APC0_prod_rate_);
    /*12*/  std::swap(one.percentage_TNF_APCa_prod_rate_ ,other.percentage_TNF_APCa_prod_rate_);
    /*13*/  std::swap(one.percentage_TNF_APCbo_prod_rate_ ,other.percentage_TNF_APCbo_prod_rate_);


    /// 6) Proliferation rates
    /*14*/  std::swap(one.APC_bound_proliferation_rate_ ,other.APC_bound_proliferation_rate_);

    /// 7) Apoptosis rates
    /*15*/  std::swap(one.APC0_apop_rate_ ,other.APC0_apop_rate_);
    /*16*/  std::swap(one.APCa_apop_rate_ ,other.APCa_apop_rate_);
    /*17*/  std::swap(one.APCbo_apop_rate_ ,other.APCbo_apop_rate_);
    /*18*/  std::swap(one.APCbl_apop_rate_ ,other.APCbl_apop_rate_);
    /*19*/  std::swap(one.APCexh_apop_rate_ ,other.APCexh_apop_rate_);

    /// 8) constant saturation of TNF for apoptosis
    /*20*/  std::swap(one.Ks_APC_m_TNF_ ,other.Ks_APC_m_TNF_);

    /// 9) conversion rates
    /*21*/  std::swap(one.APC_Ag_ ,other.APC_Ag_);
    /*22*/  std::swap(one.APC_APC_ ,other.APC_APC_);
    /*23*/  std::swap(one.APC_NK_ ,other.APC_NK_);
    /*24*/  std::swap(one.APC_LT_1_ ,other.APC_LT_1_);
    /*25*/  std::swap(one.APC_LT_2_ ,other.APC_LT_2_);
    /*26*/  std::swap(one.APC_Ab_ ,other.APC_Ab_);
    /*27*/  std::swap(one.APC_exh_ ,other.APC_exh_);

    /// 10)Saturation constant of IFN and TNF for activation
    /*28*/  std::swap(one.KsAPC_LT_ ,other.KsAPC_LT_);

    /// 11)Saturation constant of APC_LT interaction
    /*29*/  std::swap(one.APC_Ksi_ ,other.APC_Ksi_);
    /*30*/  std::swap(one.APC_Kst_ ,other.APC_Kst_);

    /// 12) Percentages of cell expressing receptor
    /*31*/  std::swap(one.APC0_expressing_receptor_ ,other.APC0_expressing_receptor_);
    /*32*/  std::swap(one.APCa_expressing_receptor_ ,other.APCa_expressing_receptor_);
    /// 13) Apoptosis rate for TNF
    /*33*/ std::swap(one.u_APC_TNF_ ,other.u_APC_TNF_);

/// NK
   /// 1) Init ratio of cells
   /*1*/ std::swap(one.init_ratio_NK_ ,other.init_ratio_NK_);

   /// 2) IFN Poductions rates of each type of NK
   /*2*/ std::swap(one.IFN_NK0_prod_rate_,other.IFN_NK0_prod_rate_);
   /*3*/ std::swap(one.IFN_NKa_prod_rate_,other.IFN_NKa_prod_rate_);
   /*4*/ std::swap(one.IFN_NKbo_prod_rate_,other.IFN_NKbo_prod_rate_);

   /// 3) TNF Poductions rates of each type of NK
   /*5*/ std::swap(one.TNF_NK0_prod_rate_,other.TNF_NK0_prod_rate_);
   /*6*/ std::swap(one.TNF_NKa_prod_rate_,other.TNF_NKa_prod_rate_);
   /*7*/ std::swap(one.TNF_NKbo_prod_rate_,other.TNF_NKbo_prod_rate_);

   /// 4) Percentages of IFN productions of each type of NK
   /*8*/ std::swap(one.percentage_IFN_NK0_prod_rate_,other.percentage_IFN_NK0_prod_rate_);
   /*9*/ std::swap(one.percentage_IFN_AgNKa_prod_rate_,other.percentage_IFN_AgNKa_prod_rate_);
   /*10*/ std::swap(one.percentage_IFN_NKbo_prod_rate_,other.percentage_IFN_NKbo_prod_rate_);

   /// 5)Percentages of TNF productions of each type of NK
   /*11*/ std::swap(one.percentage_TNF_NK0_prod_rate_,other.percentage_TNF_NK0_prod_rate_);
   /*12*/ std::swap(one.percentage_TNF_NKa_prod_rate_,other.percentage_TNF_NKa_prod_rate_);
   /*13*/ std::swap(one.percentage_TNF_NKbo_prod_rate_,other.percentage_TNF_NKbo_prod_rate_);

    /// 6) Proliferation rates
    /*13.5*/ std::swap(one.NK0_proliferation_rate_,other.NK0_proliferation_rate_);
    /*14*/ std::swap(one.NKa_proliferation_rate_,other.NKa_proliferation_rate_);
    /*15*/ std::swap(one.NKbo_proliferation_rate_,other.NKbo_proliferation_rate_);
    /*16*/ std::swap(one.NKbl_proliferation_rate_,other.NKbl_proliferation_rate_);

    /// 7) Apoptosis rates
    /*17*/ std::swap(one.NK0_apop_rate_,other.NK0_apop_rate_);
    /*18*/ std::swap(one.NKa_apop_rate_,other.NKa_apop_rate_);
    /*19*/ std::swap(one.NKbo_apop_rate_,other.NKbo_apop_rate_);
    /*20*/ std::swap(one.NKbl_apop_rate_,other.NKbl_apop_rate_);
    /*21*/ std::swap(one.NKexh_apop_rate_,other.NKexh_apop_rate_);

    /// 8) constant saturation of TNF for apoptosis
    /*22*/ std::swap(one.Ks_NK_m_TNF_,other.Ks_NK_m_TNF_);

    /// 9) conversion rates
    /*23*/ std::swap(one.KaNK_,other.KaNK_);
    /*24*/ std::swap(one.NK_NK_,other.NK_NK_);
    /*25*/ std::swap(one.NK_Ab_,other.NK_Ab_);
    /*26*/ std::swap(one.NK_exh_,other.NK_exh_);

    /// 10)Saturation constant of NK interaction for activation
    /*27*/ std::swap(one.KsAPC_NK_,other.KsAPC_NK_);

    /// 11)Saturation constant of NK_LT interaction
    /*28*/ std::swap(one.NK_Ksi_,other.NK_Ksi_);
    /*29*/ std::swap(one.NK_Kst_,other.NK_Kst_);

    /// 12) Percentages of cell expressing receptor
    /*30*/ std::swap(one.NK0_expressing_receptor_,other.NK0_expressing_receptor_);
    /*31*/ std::swap(one.NKa_expressing_receptor_,other.NKa_expressing_receptor_);

    /// 13) Apoptosis rate for TNF
    /*32*/ std::swap(one.u_NK_TNF_,other.u_NK_TNF_);

/// LT
    /// 1) Init number of LT
    /*1*/   std::swap(one.ratio_init_LTns_,other.ratio_init_LTns_);
    /*2*/   std::swap(one.ratio_initLTspecific_,other.ratio_initLTspecific_);

    /// 2) IFN Poductions rates of each type of LT
    /*3*/   std::swap(one.IFN_LTns_prod_rate_,other.ratio_initLTspecific_);
    /*4*/   std::swap(one.IFN_LTbo_prod_rate_,other.IFN_LTbo_prod_rate_);
    /*5*/   std::swap(one.IFN_LTbl_prod_rate_,other.IFN_LTbl_prod_rate_);

    /// 3) TNF Poductions rates of each type of LT
    /*6*/  std::swap(one.TNF_LTns_prod_rate_,other.TNF_LTns_prod_rate_);
    /*7*/   std::swap(one.TNF_LTbo_prod_rate_,other.TNF_LTbo_prod_rate_);
    /*8*/   std::swap(one.TNF_LTbl_prod_rate_,other.TNF_LTbl_prod_rate_);

    /// 4) Percentages of IFN productions of each type of LT
    /*9*/   std::swap(one.percentage_IFN_LTns_prod_rate_,other.percentage_IFN_LTns_prod_rate_);
    /*10*/  std::swap(one.percentage_IFN_LTbo_prod_rate_,other.percentage_IFN_LTbo_prod_rate_);
    /*11*/   std::swap(one.percentage_IFN_LTbl_prod_rate_,other.percentage_IFN_LTbl_prod_rate_);


    /// 5)Percentages of TNF productions of each type of LT
    /*12*/   std::swap(one.percentage_TNF_LTns_prod_rate_,other.percentage_TNF_LTns_prod_rate_);
    /*13*/   std::swap(one.percentage_TNF_LTbo_prod_rate_,other.percentage_TNF_LTbo_prod_rate_);
    /*14*/   std::swap(one.percentage_TNF_LTbl_prod_rate_,other.percentage_TNF_LTbl_prod_rate_);

    /// 6) Proliferation rates
    /*15*/   std::swap(one.LTns_proliferation_rate_,other.LTns_proliferation_rate_);
    /*16*/   std::swap(one.LTbo_proliferation_rate_,other.LTbo_proliferation_rate_);
    /*17*/   std::swap(one.LTbl_proliferation_rate_,other.LTbl_proliferation_rate_);

    /// 7) Apoptosis rates
    /*18*/   std::swap(one.LTns_apop_rate_,other.LTns_apop_rate_);
    /*19*/   std::swap(one.LTbo_apop_rate_,other.LTbo_apop_rate_);
    /*20*/   std::swap(one.LTbl_apop_rate_,other.LTbl_apop_rate_);
    /*21*/   std::swap(one.LTexh_apop_rate_,other.LTexh_apop_rate_);

    /// 8) constant saturation of TNF for apoptosis
    /*22*/   std::swap(one.Ks_LT_m_TNF_,other.LTexh_apop_rate_);

    /// 9) Percentages of cell expressing receptor
    /*23*/   std::swap(one.LTns_expressing_receptor_,other.LTexh_apop_rate_);

    /// 10) Apoptosis rate for TNF
    /*24*/  std::swap(one.u_LT_TNF_,other.LTexh_apop_rate_);

    /// 11) LT exh rate
    /*25*/ std::swap(one.LT_exh_rate_,other.LT_exh_rate_);

    /// 12) apoptosis related parameters
    /*26*/ std::swap(one.t_apop_meas_,other.t_apop_meas_);
    /*27*/ std::swap(one.t_duration_apoptosis_,other.t_duration_apoptosis_);

/// Media
    /*1*/  std::swap(one.TNF_deg_,other.TNF_deg_);
    /*2*/  std::swap(one.IFN_deg_,other.IFN_deg_);
    /*3*/  std::swap(one.TymidineTriteate_,other.TymidineTriteate_);
    /*4*/  std::swap(one.Prol_TymTr_,other.Prol_TymTr_);
}

std::ostream& operator<<(std::ostream& s,SimParameters p)
{   /// APC
    /*1*/   s<<"\n init_ratio_APC_ \t " <<p.init_ratio_APC_;
    /*2*/   s<<"\n IFN_APC0_prod_rate_ \t "<<p.IFN_APC0_prod_rate_;
    /*3*/   s<<"\n IFN_APCa_prod_rate_ \t"<<p.IFN_APCa_prod_rate_;
    /*4*/   s<<"\n IFN_APCbo_prod_rate_ \t"<<p.IFN_APCbo_prod_rate_;
    /*5*/   s<<"\n TNF_APC0_prod_rate_ \t"<<p.TNF_APC0_prod_rate_;
    /*6*/   s<<"\n TNF_APCa_prod_rate_ \t"<<p.TNF_APCa_prod_rate_;
    /*7*/   s<<"\n TNF_APCbo_prod_rate_ \t"<<p.TNF_APCbo_prod_rate_;
    /*8*/   s<<"\n percentage_IFN_APC0_prod_rate_ \t"<<p.percentage_IFN_APC0_prod_rate_;
    /*9*/   s<<"\n percentage_IFN_APCa_prod_rate_ \t"<<p.percentage_IFN_APCa_prod_rate_;
    /*10*/  s<<"\n percentage_IFN_APCbo_prod_rate_\t"<<p.percentage_IFN_APCbo_prod_rate_;
    /*11*/  s<<"\n percentage_TNF_APC0_prod_rate_\t"<<p.percentage_TNF_APC0_prod_rate_;
    /*12*/  s<<"\n percentage_TNF_APCa_prod_rate_\t"<<p.percentage_TNF_APCa_prod_rate_;
    /*13*/  s<<"\n percentage_TNF_APCbo_prod_rate_\t"<<p.percentage_TNF_APCbo_prod_rate_;
    /*14*/  s<<"\n APC_bound_proliferation_rate_\t"<<p.APC_bound_proliferation_rate_;
    /*15*/  s<<"\n APC0_apop_rate_\t"<<p.APC0_apop_rate_;
    /*16*/  s<<"\n APCa_apop_rate_\t"<<p.APCa_apop_rate_;
    /*17*/  s<<"\n APCbo_apop_rate_\t"<<p.APCbo_apop_rate_;
    /*18*/  s<<"\n APCbl_apop_rate_\t"<<p.APCbl_apop_rate_;
    /*19*/  s<<"\n APCexh_apop_rate_\t"<<p.APCexh_apop_rate_;
    /*20*/  s<<"\n Ks_APC_m_TNF_d\t"<<p.Ks_APC_m_TNF_;
    /*21*/  s<<"\n APC_Ag_\t"<<p.APC_Ag_;
    /*22*/  s<<"\n APC_APC_\t"<<p.APC_APC_;
    /*23*/  s<<"\n APC_NK_\t"<<p.APC_NK_;
    /*24*/  s<<"\n APC_LT_1_\t"<<p.APC_LT_1_;
    /*25*/  s<<"\n APC_LT_2_\t"<<p.APC_LT_2_;
    /*26*/  s<<"\n APC_Ab_\t"<<p.APC_Ab_;
    /*27*/  s<<"\n  APC_exh_\t"<<p.APC_exh_;
    /*28*/  s<<"\n KsAPC_LT_\t"<<p.KsAPC_LT_;
    /*29*/  s<<"\n Ksi_\t"<<p.APC_Ksi_;
    /*30*/  s<<"\n Kst_\t"<<p.APC_Kst_;
    /*31*/  s<<"\n APC0_expressing_receptor_\t"<<p.APC0_expressing_receptor_;
    /*32*/  s<<"\n APCa_expressing_receptor_\t"<<p.APCa_expressing_receptor_;
    /*33*/  s<<"\n u_APC_TNF_\t"<<p.u_APC_TNF_;
    /// NK
    /*1*/   s<<"\n init_ratio_NK_\t"<<p.init_ratio_NK_;
    /*2*/   s<<"\n IFN_NK0_prod_rate_t"<<p.IFN_NK0_prod_rate_;
    /*3*/   s<<"\n IFN_NKa_prod_rate_\t"<<p.IFN_NKa_prod_rate_;
    /*4*/   s<<"\n IFN_NKbo_prod_rate_\t"<<p.IFN_NKbo_prod_rate_;
    /*5*/   s<<"\n TNF_NK0_prod_rate_\t"<<p.TNF_NK0_prod_rate_;
    /*6*/   s<<"\n TNF_NKa_prod_rate_\t"<<p.TNF_NKa_prod_rate_;
    /*7*/   s<<"\n TNF_NKbo_prod_rate_\t"<<p.TNF_NKbo_prod_rate_;
    /*8*/   s<<"\n percentage_IFN_NK0_prod_rate_\t"<<p.percentage_IFN_NK0_prod_rate_;
    /*9*/   s<<"\n percentage_IFN_AgNKa_prod_rate_\t"<<p.percentage_IFN_AgNKa_prod_rate_;
    /*10*/  s<<"\n percentage_IFN_NKbo_prod_rate_\t"<<p.percentage_IFN_NKbo_prod_rate_;
    /*11*/  s<<"\n percentage_TNF_NK0_prod_rate_\t"<<p.percentage_TNF_NK0_prod_rate_;
    /*12*/  s<<"\n percentage_TNF_NKa_prod_rate_\t"<<p.percentage_TNF_NKa_prod_rate_;
    /*13*/  s<<"\n percentage_TNF_NKbo_prod_rate_\t"<<p.percentage_TNF_NKbo_prod_rate_;
    /*13.5*/  s<<"\n NK0_proliferation_rate_\t"<<p.NK0_proliferation_rate_;
    /*14*/  s<<"\n NKa_proliferation_rate_\t"<<p.NKa_proliferation_rate_;
    /*15*/  s<<"\n NKbo_proliferation_rate_\t"<<p.NKbo_proliferation_rate_;
    /*16*/  s<<"\n NKbl_proliferation_rate_\t"<<p.NKbl_proliferation_rate_;
    /*17*/  s<<"\n NK0_apop_rate_\t"<<p.NK0_apop_rate_;
    /*18*/  s<<"\n NKa_apop_rate_\t"<<p.NKa_apop_rate_;
    /*19*/  s<<"\n NKbo_apop_rate_\t"<<p.NKbo_apop_rate_;
    /*20*/  s<<"\n NKbl_apop_rate_t"<<p.NKbl_apop_rate_;
    /*21*/  s<<"\n NKexh_apop_rate_\t"<<p.NKexh_apop_rate_;
    /*22*/  s<<"\n Ks_NK_m_TNF_\t"<<p.Ks_NK_m_TNF_;
    /*23*/  s<<"\n KaNK_\t"<<p.KaNK_;
    /*24*/  s<<"\n NK_NK_\t"<<p.NK_NK_;
    /*25*/  s<<"\n NK_Ab\t"<<p.NK_Ab_;
    /*26*/  s<<"\n NK_exh_\t"<<p.NK_exh_;
    /*27*/  s<<"\n KsAPC_NK_\t"<<p.KsAPC_NK_;
    /*28*/  s<<"\n Ksi_\t"<<p.NK_Ksi_;
    /*29*/  s<<"\n Kst_\t"<<p.NK_Kst_;
    /*30*/  s<<"\n NK0_expressing_receptor_\t"<<p.NK0_expressing_receptor_;
    /*31*/  s<<"\n NKa_expressing_receptor_\t"<<p.NKa_expressing_receptor_;
    /*32*/  s<<"\n u_NK_TNF_\t " <<p.u_NK_TNF_;
    /// LT
    /*1*/   s<<"\n ratio_init_LTns_\t"<<p. ratio_init_LTns_;
    /*2*/   s<<"\n ratio_initLTspecific_\t"<<p. ratio_initLTspecific_;
    /*3*/   s<<"\n IFN_LTns_prod_rate_\t"<<p. ratio_initLTspecific_;
    /*4*/   s<<"\n IFN_LTbo_prod_rate_\t"<<p. IFN_LTbo_prod_rate_;
    /*5*/   s<<"\n IFN_LTbl_prod_rate_\t"<<p. IFN_LTbl_prod_rate_;
    /*6*/   s<<"\n TNF_LTns_prod_rate_\t"<<p. TNF_LTns_prod_rate_;
    /*7*/   s<<"\n TNF_LTbo_prod_rate_\t"<<p. TNF_LTbo_prod_rate_;
    /*8*/   s<<"\n TNF_LTbl_prod_rate_\t"<<p. TNF_LTbl_prod_rate_;
    /*9*/   s<<"\n percentage_IFN_LTns_prod_rate_\t"<<p. percentage_IFN_LTns_prod_rate_;
    /*10*/  s<<"\n percentage_IFN_LTbo_prod_rate_\t"<<p. percentage_IFN_LTbo_prod_rate_;
    /*11*/  s<<"\n percentage_IFN_LTbl_prod_rate_\t"<<p. percentage_IFN_LTbl_prod_rate_;
    /*12*/  s<<"\n percentage_TNF_LTns_prod_rate_\t"<<p. percentage_TNF_LTns_prod_rate_;
    /*13*/  s<<"\n percentage_TNF_LTbo_prod_rate_\t"<<p. percentage_TNF_LTbo_prod_rate_;
    /*14*/  s<<"\n percentage_TNF_LTbl_prod_rate_\t"<<p. percentage_TNF_LTbl_prod_rate_;
    /*15*/  s<<"\n LTns_proliferation_rate_\t"<<p. LTns_proliferation_rate_;
    /*16*/  s<<"\n LTbo_proliferation_rate_\t"<<p. LTbo_proliferation_rate_;
    /*17*/  s<<"\n LTbl_proliferation_rate_\t"<<p. LTbl_proliferation_rate_;
    /*18*/  s<<"\n LTns_apop_rate_\t"<<p. LTns_apop_rate_;
    /*19*/  s<<"\n LTbo_apop_rate_\t"<<p. LTbo_apop_rate_;
    /*20*/  s<<"\n LTbl_apop_rate_\t"<<p. LTbl_apop_rate_;
    /*21*/  s<<"\n LTexh_apop_rate_\t"<<p. LTexh_apop_rate_;
    /*22*/  s<<"\n Ks_LT_m_TNF_\t"<<p. LTexh_apop_rate_;
    /*23*/  s<<"\n LTns_expressing_receptor_\t"<<p. LTexh_apop_rate_;
    /*24*/  s<<"\n u_LT_TNF_\t"<<p.LTexh_apop_rate_;
    /*25*/  s<<"\n LT_exh_rate_\t"<<p.LT_exh_rate_;
    /*26*/  s<<"\n t_apop_meas_\t"<<p.t_apop_meas_;
    /*27*/  s<<"\n t_duration_apoptosis_\t"<<p.t_duration_apoptosis_;
    /// Media
    /*1*/  s<<"\n TNF_deg_\t"<<p.TNF_deg_;
    /*2*/  s<<"\n IFN_deg_\t"<<p.IFN_deg_;
    /*3*/  s<<"\n TymidineTriteate_\t"<<p.TymidineTriteate_;
    /*4*/  s<<"\n Prol_TymTr_\t"<<p.Prol_TymTr_;

    return s;

}
std::vector<double> SimParameters::getRandomParameters(double range)const
{
    std::vector<double> result=getParameters();
    srand ( time(NULL) );

    for (std::size_t i=0; i<result.size();++i)
    {
        result[i]+=2*range*((1.0*rand())/RAND_MAX-0.5);
    }
    return result;
}