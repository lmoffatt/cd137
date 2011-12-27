#include <fstream>
#include "Includes/experiment1.h"
#include "Includes/Cell_simulator.h"
#include "Includes/Treatment.h"
#include "Includes/Experiment.h"
#include "Includes/OptimizationResults.h"

void Bayes()
{   Treatment  med;
    med.Ag=0.0;
    med.Ab=0.0;
    med.sim_duration_d=120;
    med.time_step_d=1.0/120;
    med.init_cells=1e6;

    Treatment Mtb;
    Mtb.Ag=10.0;
    Mtb.Ab=0.0;
    Mtb.sim_duration_d=120;
    Mtb.time_step_d=1.0/120;
    Mtb.init_cells=1e6;

    Treatment block;
    block.Ag=10.0;
    block.Ab=10.0;
    block.sim_duration_d=120;
    block.time_step_d=1.0/12;
    block.init_cells=1e6;

    Results MediaRes ("media");
    Results MtbRes("mtb");
    Results blockRes ("block");
    Experiment E;
    E.push_back(block,blockRes);
    E.push_back(Mtb,MtbRes);
    E.push_back(med,MediaRes);


    SimParameters sp;
    sp.mode_="FULL";
    /// APC
    /// 1) Init ratio of cells
    /*1*/ sp.init_ratio_APC_=1e5/1e6;

    /// 2) IFN Poductions rates of each type of APC
    /*2*/ sp.IFN_APC0_prod_rate_=1e-10;
    /*3*/ sp.IFN_APCa_prod_rate_=1e-10;
    /*4*/ sp.IFN_APCbo_prod_rate_=1e-8;


    /// 3) TNF Poductions rates of each type of APC
    /*5*/ sp.TNF_APC0_prod_rate_=1e-4;
    /*6*/ sp.TNF_APCa_prod_rate_=1e-4;
    /*7*/ sp.TNF_APCbo_prod_rate_=1e-2;


    /// 4) Percentages of IFN productions of each type of APC
    /*8*/ sp.percentage_IFN_APC0_prod_rate_=0.001;
    /*9*/ sp.percentage_IFN_APCa_prod_rate_=0.01;
    /*10*/ sp.percentage_IFN_APCbo_prod_rate_=0.035;

    /// 5)Percentages of TNF productions of each type of APC
    /*11*/ sp.percentage_TNF_APC0_prod_rate_=0.001;
    /*12*/ sp.percentage_TNF_APCa_prod_rate_=0.01;
    /*13*/ sp.percentage_TNF_APCbo_prod_rate_=0.035;


    /// 6) Proliferation rates
    /*14*/ sp.APC_bound_proliferation_rate_=1e-6;

    /// 7) Apoptosis rates
    /*15*/ sp.APC0_apop_rate_=1e-16;
    /*16*/ sp.APCa_apop_rate_=1e-16;
    /*17*/ sp.APCbo_apop_rate_=1e-16;
    /*18*/ sp.APCbl_apop_rate_=1e-16;
    /*19*/ sp.APCexh_apop_rate_=1e-16;

    /// 8) constant saturation of TNF for apoptosis
    /*20*/ sp.Ks_APC_m_TNF_=10.0e10;

    /// 9) conversion rates
    /*21*/ sp.APC_Ag_=1e-6;
    /*22*/ sp.APC_APC_=1e-6;
    /*23*/ sp.APC_NK_=1e-6;
    /*24*/ sp.APC_LT_1_=1e-6;
    /*25*/ sp.APC_LT_2_=1e-6;
    /*26*/ sp.APC_Ab_=1e-6;
    /*27*/ sp.APC_exh_=1e-6;

    /// 10)Saturation constant of IFN and TNF for activation
    /*28*/ sp.KsAPC_LT_=1e-6;

    /// 11)Saturation constant of APC_LT interaction
    /*29*/ sp.APC_Ksi_=1e-6;
    /*30*/ sp.APC_Kst_=1e-6;

    /// 12) Percentages of cell expressing receptor
    /*31*/ sp.APC0_expressing_receptor_=0.01;
    /*32*/ sp.APCa_expressing_receptor_=0.01;
    /// 13) Apoptosis rate for TNF
    /*33*/ sp.u_APC_TNF_ =1e-12;

    /// NK
    /// 1) Init ratio of cells
    /*1*/  sp.init_ratio_NK_=1e5/1e6;

    /// 2) IFN Poductions rates of each type of NK
    /*2*/  sp.IFN_NK0_prod_rate_=1e-10;
    /*3*/  sp.IFN_NKa_prod_rate_=1e-2;
    /*4*/  sp.IFN_NKbo_prod_rate_=1e-2;

    /// 3) TNF Poductions rates of each type of NK
    /*5*/  sp.TNF_NK0_prod_rate_=1e-10;
    /*6*/  sp.TNF_NKa_prod_rate_=1e-10;
    /*7*/  sp.TNF_NKbo_prod_rate_=1e-10;

    /// 4) Percentages of IFN productions of each type of NK
    /*8*/  sp.percentage_IFN_NK0_prod_rate_=0.01;
    /*9*/  sp.percentage_IFN_AgNKa_prod_rate_=0.01;
    /*10*/  sp.percentage_IFN_NKbo_prod_rate_=0.01;

    /// 5)Percentages of TNF productions of each type of NK
    /*11*/  sp.percentage_TNF_NK0_prod_rate_=0.01;
    /*12*/  sp.percentage_TNF_NKa_prod_rate_=0.01;
    /*13*/  sp.percentage_TNF_NKbo_prod_rate_=0.01;

    /// 6) Proliferation rates
    /*13.5*/  sp.NK0_proliferation_rate_=1e-6;
    /*14*/  sp.NKa_proliferation_rate_=1e-6;
    /*15*/  sp.NKbo_proliferation_rate_=1e-6;
    /*16*/  sp.NKbl_proliferation_rate_=1e-6;

    /// 7) Apoptosis rates
    /*17*/  sp.NK0_apop_rate_=1e-16;
    /*18*/  sp.NKa_apop_rate_=1e-16;
    /*19*/  sp.NKbo_apop_rate_=1e-16;
    /*20*/  sp.NKbl_apop_rate_=1e-16;
    /*21*/  sp.NKexh_apop_rate_=1e-16;

    /// 8) constant saturation of TNF for apoptosis
    /*22*/  sp.Ks_NK_m_TNF_=1e6;

    /// 9) conversion rates
    /*23*/  sp.KaNK_=1e-6;
    /*24*/  sp.NK_NK_=1e-6;
    /*25*/  sp.NK_Ab_=1e-6;
    /*26*/  sp.NK_exh_=1e-6;

    /// 10)Saturation constant of APC NK interaction for activation
    /*27*/  sp.KsAPC_NK_=1e-6;

    /// 11)Saturation constant of NK_LT interaction
    /*28*/  sp.NK_Ksi_=1e-6;
    /*29*/  sp.NK_Kst_=1e-6;

    /// 12) Percentages of cell expressing receptor
    /*30*/  sp.NK0_expressing_receptor_=0.01;
    /*31*/  sp.NKa_expressing_receptor_=0.01;

    /// 13) Apoptosis rate for TNF
    /*32*/  sp.u_NK_TNF_=1e-16;

    /// LT
    /// 1) Init number of LT
       /*1*/  sp.ratio_init_LTns_=0.9e6/1e6;
       /*2*/  sp.ratio_initLTspecific_=1e2/1e6;

    /// 2) IFN Poductions rates of each type of LT
       /*3*/  sp.IFN_LTns_prod_rate_=1e-2;
       /*4*/  sp.IFN_LTbo_prod_rate_=1e-2;
       /*5*/  sp.IFN_LTbl_prod_rate_=1e-2;

   /// 3) TNF Poductions rates of each type of LT
       /*6*/  sp.TNF_LTns_prod_rate_=1e-10;
       /*7*/  sp.TNF_LTbo_prod_rate_=1e-10;
       /*8*/  sp.TNF_LTbl_prod_rate_=1e-10;

   /// 4) Percentages of IFN productions of each type of LT
       /*9*/  sp.percentage_IFN_LTns_prod_rate_=0.05;
       /*10*/  sp.percentage_IFN_LTbo_prod_rate_=0.8;
       /*11*/  sp.percentage_IFN_LTbl_prod_rate_=0.8;


   /// 5)Percentages of TNF productions of each type of LT
       /*12*/  sp.percentage_TNF_LTns_prod_rate_=0.05;
       /*13*/  sp.percentage_TNF_LTbo_prod_rate_=0.8;
       /*14*/  sp.percentage_TNF_LTbl_prod_rate_=0.8;

   /// 6) Proliferation rates
       /*15*/  sp.LTns_proliferation_rate_=1e-2;
       /*16*/  sp.LTbo_proliferation_rate_=1e-2;
       /*17*/  sp.LTbl_proliferation_rate_=1e-2;

   /// 7) Apoptosis rates
       /*18*/  sp.LTns_apop_rate_=1e-16;
       /*19*/  sp.LTbo_apop_rate_=1e-16;
       /*20*/  sp.LTbl_apop_rate_=1e-16;
       /*21*/  sp.LTexh_apop_rate_=1e-16;

   /// 8) constant saturation of TNF for apoptosis
       /*22*/  sp.Ks_LT_m_TNF_=1e-10;

   /// 9) Percentages of cell expressing receptor
       /*23*/  sp.LTns_expressing_receptor_=1e-10;

   /// 10) Apoptosis rate for TNF
       /*24*/  sp.u_LT_TNF_=1e-16;

   /// 11) LT exh rate
       /*25*/ sp.LT_exh_rate_=1e-10;

   /// 12) apoptosis related parameters
       /*26*/ sp.t_apop_meas_=120;
       /*27*/ sp.t_duration_apoptosis_=2;

    /// Media
    /*1*/ sp.TNF_deg_=0.5;
    /*2*/ sp.IFN_deg_=0.5;
    /*4*/ sp.Prol_TymTr_=0.5;




    Cell_simulator cell(sp, E);
    SimParameters perturbedPar=sp.applyParameters(sp.getRandomParameters(0));
    Experiment simulExp=cell.Simulate(perturbedPar ,E);

    std::cout<<E.Result_i(0);

     //Modificar num iteracines
    //simulExp: simulado  E:experimental
    //OptimizationResults O=cell.Optimize(sp,sp,simulExp,1,500);

     OptimizationResults O=cell.Optimize(sp,sp,E,1,10);




    /*std::cout<<"Lower to 1\n";
    O=cell.Optimize(O.OptimalParameters(),E,1,200);
    std::cout<<"Lower to 0.5\n";

    O=cell.Optimize(O.OptimalParameters(),E,0.5,200);
*/


    std::ofstream f;
    f.open("resultssim.txt");
    f<<O;

    cell.applyParameters(O.OptimalParameters(),E.Treatment_i(0));
    cell.run();
    cell.applyParameters(O.OptimalParameters(),E.Treatment_i(1));
    cell.run();
    cell.applyParameters(O.OptimalParameters(),E.Treatment_i(2));
    cell.run();
    f.close();
}


void TestParametersIO()
{
    Parameters sp=Cell_simulator::getStandardParameters();
    std::string filename="parameters.txt";
    std::ofstream f;
    f.open(filename.c_str(),std::ios_base::app);
    f<<sp;
    f.close();


    std::ifstream fi;
    fi.open(filename.c_str(),std::ios_base::in);

    Parameters sp2;
    fi>>sp2;
    fi.close();
    std::cout<<sp2;



}



void BayesParameters()
{   Treatment  med;
    med.Ag=0.0;
    med.Ab=0.0;
    med.sim_duration_d=120;
    med.time_step_d=1.0/1200;
    med.init_cells=1.0e6;
    med.t_apop_meas_d=119;

    Treatment Mtb;
    Mtb.Ag=10.0;
    Mtb.Ab=0.0;
    Mtb.sim_duration_d=120;
    Mtb.time_step_d=1.0/1200;
    Mtb.init_cells=1.0e6;
    Mtb.t_apop_meas_d=119;

    Treatment block;
    block.Ag=10.0;
    block.Ab=10.0;
    block.sim_duration_d=120;
    block.time_step_d=1.0/1200;
    block.init_cells=1.0e6;
    block.t_apop_meas_d=119;

    Results MediaRes ("media");
    Results MtbRes("mtb");
    Results blockRes ("block");
    Experiment E;
    E.push_back(block,blockRes);
    E.push_back(Mtb,MtbRes);
    E.push_back(med,MediaRes);

    std::cout<<E;
    char ch;
  //  std::cin>>ch;

    Parameters sp=Cell_simulator::getStandardParameters();


    Cell_simulator cell(sp, E);

     //Modificar num iteracines
    //simulExp: simulado  E:experimental
    //OptimizationResults O=cell.Optimize(sp,sp,simulExp,1,500);

    cell.Optimize(sp,E,"optimParam.txt");




    /*std::cout<<"Lower to 1\n";
    O=cell.Optimize(O.OptimalParameters(),E,1,200);
    std::cout<<"Lower to 0.5\n";

    O=cell.Optimize(O.OptimalParameters(),E,0.5,200);
*/



}


void compareTimeStep(
        )
{   Treatment  med;
    med.Ag=0.0;
    med.Ab=0.0;
    med.sim_duration_d=120;
    med.time_step_d=1.0/120*4.0;
    med.init_cells=1.0e6;
    med.t_apop_meas_d=119;

    Treatment Mtb;
    Mtb.Ag=10.0;
    Mtb.Ab=0.0;
    Mtb.sim_duration_d=120;
    Mtb.time_step_d=1.0/120*4.0;
    Mtb.init_cells=1.0e6;
    Mtb.t_apop_meas_d=119;

    Treatment block;
    block.Ag=10.0;
    block.Ab=10.0;
    block.sim_duration_d=120;
    block.time_step_d=1.0/120*4.0;
    block.init_cells=1.0e6;
    block.t_apop_meas_d=119;

    Results MediaRes ("media");
    Results MtbRes("mtb");
    Results blockRes ("block");
    Experiment E;
    E.push_back(block,blockRes);
    E.push_back(Mtb,MtbRes);
    E.push_back(med,MediaRes);

    Parameters sp=Cell_simulator::getStandardParameters();


    Cell_simulator cell(sp, E);

    std::string filename="comparetimestep.txt";
    std::ofstream f;
    f.open(filename.c_str(),std::ios_base::app);

    cell.run(f,sp);
    f.close();



     //Modificar num iteracines
    //simulExp: simulado  E:experimental
    //OptimizationResults O=cell.Optimize(sp,sp,simulExp,1,500);

    cell.Optimize(sp,SimE,"testOptimization.txt");






}




































































