#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "Includes/experiment1.h"
#include "Includes/Cell_simulator.h"
#include "Includes/Treatment.h"
#include "Includes/Experiment.h"
#include "Includes/OptimizationResults.h"

void loadModel()

{
    srand ( time(NULL) );

    Treatment  med;
    med.Ag=0.0;
    med.Ab=0.0;
    med.sim_duration_d=120;
    med.time_step_d=1.0/120;
    med.init_cells=1.0e6;
    med.t_apop_meas_d=119;

    Treatment Mtb;
    Mtb.Ag=10.0;
    Mtb.Ab=0.0;
    Mtb.sim_duration_d=120;
    Mtb.time_step_d=1.0/120;
    Mtb.init_cells=1.0e6;
    Mtb.t_apop_meas_d=119;

    Treatment block;
    block.Ag=10.0;
    block.Ab=1.0;
    block.sim_duration_d=120;
    block.time_step_d=1.0/120;
    block.init_cells=1.0e6;
    block.t_apop_meas_d=119;

    Results MediaRes ("media");
    Results MtbRes("mtb");
    Results blockRes ("block");
    Experiment E;
    E.push_back(med,MediaRes);
    E.push_back(block,blockRes);
    E.push_back(Mtb,MtbRes);



    std::string filenamePrior="priorMODEL.txt";

    std::ifstream f;
    f.open(filenamePrior.c_str());

    Parameters prior;
    f>>prior;

    std::cout<<prior;
    f.close();

 /*   std::string filenameCurrent="curentMODEL.txt";

    f.open(filenameCurrent.c_str());

    Parameters current;
    f>>current;

    std::cout<<current;
    f.close();

*/

    Cell_simulator cell(prior,prior, E);



     //Modificar num iteracines
    //simulExp: simulado  E:experimental
    //OptimizationResults O=cell.Optimize(sp,sp,simulExp,1,500);

    if (1)
    {
        if (true)
            cell.Optimize(prior,E,"MODELOptimization.txt");
        else if (false)
        {
            BayesIteration b(&cell,prior,&E,"ModeloOptimizationCont.txt");

            std::string filenameStartingParameter="resultMODEL.txt";

            std::ifstream f;
            f.open(filenameStartingParameter.c_str());

            Parameters seedPar;
            f>>seedPar;

            std::cout<<seedPar;
            f.close();

            b.getPosterior(seedPar);
        }
        else
        {
            BayesIteration b(&cell,prior,&E,"ModeloOptimizationContRand.txt");

            std::string filenameStartingParameter="resultMODEL.txt";

            std::ifstream f;
            f.open(filenameStartingParameter.c_str());

            Parameters seedPar;
            f>>seedPar;

            std::cout<<seedPar;
            f.close();
            double factor=0.1;
            std::size_t numseeds=10;
            double probParameterChange=1;
            b.getPosterior(seedPar,factor,numseeds,probParameterChange);

        }


    }



    else
    {
        std::string filename="MODEL_Run.txt";
        std::ofstream f;
        f.open(filename.c_str());


        cell.run(f,prior);
        f.close();

    }




}









//void Bayes()
//{   Treatment  med;
//    med.Ag=0.0;
//    med.Ab=0.0;
//    med.sim_duration_d=120;
//    med.time_step_d=1.0/120;
//    med.init_cells=1e6;

//    Treatment Mtb;
//    Mtb.Ag=10.0;
//    Mtb.Ab=0.0;
//    Mtb.sim_duration_d=120;
//    Mtb.time_step_d=1.0/120;
//    Mtb.init_cells=1e6;

//    Treatment block;
//    block.Ag=10.0;
//    block.Ab=10.0;
//    block.sim_duration_d=120;
//    block.time_step_d=1.0/12;
//    block.init_cells=1e6;

//    Results MediaRes ("media");
//    Results MtbRes("mtb");
//    Results blockRes ("block");
//    Experiment E;
//    E.push_back(block,blockRes);
//    E.push_back(Mtb,MtbRes);
//    E.push_back(med,MediaRes);


//    SimParameters sp;
//    sp.mode_="FULL";
//    /// APC
//    /// 1) Init ratio of cells
//    /*1*/ sp.init_ratio_APC_=1e5/1e6;

//    /// 2) IFN Poductions rates of each type of APC
//    /*2*/ sp.IFN_APC0_prod_rate_=1e-10;
//    /*3*/ sp.IFN_APCa_prod_rate_=1e-10;
//    /*4*/ sp.IFN_APCbo_prod_rate_=1e-8;


//    /// 3) TNF Poductions rates of each type of APC
//    /*5*/ sp.TNF_APC0_prod_rate_=1e-4;
//    /*6*/ sp.TNF_APCa_prod_rate_=1e-4;
//    /*7*/ sp.APC_TNF_Induction_CD137=1e-2;


//    /// 4) Percentages of IFN productions of each type of APC
//    /*8*/ sp.percentage_IFN_APC0_prod_rate_=0.001;
//    /*9*/ sp.percentage_IFN_APCa_prod_rate_=0.01;
//    /*10*/ sp.percentage_IFN_APCbo_prod_rate_=0.035;

//    /// 5)Percentages of TNF productions of each type of APC
//    /*11*/ sp.percentage_TNF_APC0_prod_rate_=0.001;
//    /*12*/ sp.percentage_TNF_APCa_prod_rate_=0.01;
//    /*13*/ sp.percentage_APC_TNF_Induction_CD137=0.035;


//    /// 6) Proliferation rates
//    /*14*/ sp.APC_bound_proliferation_rate_=1e-6;

//    /// 7) Apoptosis rates
//    /*15*/ sp.APC0_apop_rate_=1e-16;
//    /*16*/ sp.APCa_apop_rate_=1e-16;
//    /*17*/ sp.APCbo_apop_rate_=1e-16;
//    /*18*/ sp.APCbl_apop_rate_=1e-16;
////    /*19*/ sp.APCexh_apop_rate_=1e-16;

//    /// 8) constant saturation of TNF for apoptosis
//    /*20*/ sp.Ks_APC_m_TNF_=10.0e10;

//    /// 9) conversion rates
//    /*21*/ sp.APC_Ag_=1e-6;
//    /*22*/ sp.APC_APC_=1e-6;
//    /*23*/ sp.APC_NK_=1e-6;
//    /*24*/ sp.APC_LT_1_=1e-6;
//    /*25*/ sp.APC_LT_2_=1e-6;
//    /*26*/ sp.APC_Ab_=1e-6;
////    /*27*/ sp.APC_exh_=1e-6;

//    /// 10)Saturation constant of IFN and TNF for activation
//    /*28*/ sp.KsAPC_LT_=1e-6;

//    /// 11)Saturation constant of APC_LT interaction
//    /*29*/ sp.APC_Ksi_=1e-6;
//    /*30*/ sp.APC_Kst_=1e-6;

//    /// 12) Percentages of cell expressing receptor
//    /*31*/ sp.APC0_expressing_receptor_=0.01;
//    /*32*/ sp.APCa_expressing_receptor_=0.01;
//    /// 13) Apoptosis rate for TNF
//    /*33*/ sp.u_APC_TNF_ =1e-12;

//    /// NK
//    /// 1) Init ratio of cells
//    /*1*/  sp.init_ratio_NK_=1e5/1e6;

//    /// 2) IFN Poductions rates of each type of NK
//    /*2*/  sp.IFN_NK0_prod_rate_=1e-10;
//    /*3*/  sp.IFN_NKa_prod_rate_=1e-2;
//    /*4*/  sp.IFN_NKbo_prod_rate_=1e-2;

//    /// 3) TNF Poductions rates of each type of NK
//    /*5*/  sp.TNF_NK0_prod_rate_=1e-10;
//    /*6*/  sp.TNF_NKa_prod_rate_=1e-10;
//    /*7*/  sp.TNF_NKbo_prod_rate_=1e-10;

//    /// 4) Percentages of IFN productions of each type of NK
//    /*8*/  sp.percentage_IFN_NK0_prod_rate_=0.01;
//    /*9*/  sp.percentage_IFN_AgNKa_prod_rate_=0.01;
//    /*10*/  sp.percentage_IFN_NKbo_prod_rate_=0.01;

//    /// 5)Percentages of TNF productions of each type of NK
//    /*11*/  sp.percentage_TNF_NK0_prod_rate_=0.01;
//    /*12*/  sp.percentage_TNF_NKa_prod_rate_=0.01;
//    /*13*/  sp.percentage_TNF_NKbo_prod_rate_=0.01;

//    /// 6) Proliferation rates
//    /*13.5*/  sp.NK0_proliferation_rate_=1e-6;
//    /*14*/  sp.NKa_proliferation_rate_=1e-6;
//    /*15*/  sp.NKbo_proliferation_rate_=1e-6;
//    /*16*/  sp.NKbl_proliferation_rate_=1e-6;

//    /// 7) Apoptosis rates
//    /*17*/  sp.NK0_apop_rate_=1e-16;
//    /*18*/  sp.NKa_apop_rate_=1e-16;
//    /*19*/  sp.NKbo_apop_rate_=1e-16;
//    /*20*/  sp.NKbl_apop_rate_=1e-16;
////    /*21*/  sp.NKexh_apop_rate_=1e-16;

//    /// 8) constant saturation of TNF for apoptosis
//    /*22*/  sp.Ks_NK_m_TNF_=1e6;

//    /// 9) conversion rates
//    /*23*/  sp.KaNK_=1e-6;
//    /*24*/  sp.NK_NK_=1e-6;
//    /*25*/  sp.NK_Ab_=1e-6;
////    /*26*/  sp.NK_exh_=1e-6;

//    /// 10)Saturation constant of APC NK interaction for activation
//    /*27*/  sp.KsAPC_NK_=1e-6;

//    /// 11)Saturation constant of NK_LT interaction
//    /*28*/  sp.NK_Ksi_=1e-6;
//    /*29*/  sp.NK_Kst_=1e-6;

//    /// 12) Percentages of cell expressing receptor
//    /*30*/  sp.NK0_expressing_receptor_=0.01;
//    /*31*/  sp.NKa_expressing_receptor_=0.01;

//    /// 13) Apoptosis rate for TNF
//    /*32*/  sp.u_NK_TNF_=1e-16;

//    /// LT
//    /// 1) Init number of LT
//       /*1*/  sp.ratio_init_LTns_=0.9e6/1e6;
//       /*2*/  sp.ratio_initLTspecific_=1e2/1e6;

//    /// 2) IFN Poductions rates of each type of LT
//       /*3*/  sp.IFN_LTns_prod_rate_=1e-2;
//       /*4*/  sp.IFN_LTbo_prod_rate_=1e-2;
//       /*5*/  sp.IFN_LTbl_prod_rate_=1e-2;

//   /// 3) TNF Poductions rates of each type of LT
//       /*6*/  sp.TNF_LTns_prod_rate_=1e-10;
//       /*7*/  sp.TNF_LTbo_prod_rate_=1e-10;
//       /*8*/  sp.TNF_LTbl_prod_rate_=1e-10;

//   /// 4) Percentages of IFN productions of each type of LT
//       /*9*/  sp.percentage_IFN_LTns_prod_rate_=0.05;
//       /*10*/  sp.percentage_IFN_LTbo_prod_rate_=0.8;
//       /*11*/  sp.percentage_IFN_LTbl_prod_rate_=0.8;


//   /// 5)Percentages of TNF productions of each type of LT
//       /*12*/  sp.percentage_TNF_LTns_prod_rate_=0.05;
//       /*13*/  sp.percentage_TNF_LTbo_prod_rate_=0.8;
//       /*14*/  sp.percentage_TNF_LTbl_prod_rate_=0.8;

//   /// 6) Proliferation rates
//       /*15*/  sp.LTns_proliferation_rate_=1e-2;
//       /*16*/  sp.LTbo_proliferation_rate_=1e-2;
//       /*17*/  sp.LTbl_proliferation_rate_=1e-2;

//   /// 7) Apoptosis rates
//       /*18*/  sp.LTns_apop_rate_=1e-16;
//       /*19*/  sp.LTbo_apop_rate_=1e-16;
//       /*20*/  sp.LTbl_apop_rate_=1e-16;
////       /*21*/  sp.LTexh_apop_rate_=1e-16;

//   /// 8) constant saturation of TNF for apoptosis
//       /*22*/  sp.Ks_LT_m_TNF_=1e-10;

//   /// 9) Percentages of cell expressing receptor
//       /*23*/  sp.LTns_expressing_receptor_=1e-10;

//   /// 10) Apoptosis rate for TNF
//       /*24*/  sp.u_LT_TNF_=1e-16;

////   /// 11) LT exh rate
////       /*25*/ sp.LT_exh_rate_=1e-10;

//   /// 12) apoptosis related parameters
//       /*26*/ sp.t_apop_meas_=120;
//       /*27*/ sp.t_duration_apoptosis_=2;
//              sp.LT_Ab_=0.1;

//    /// Media
//    /*1*/ sp.TNF_deg_=0.5;
//    /*2*/ sp.IFN_deg_=0.5;
//          sp.Ag_deg_=0.5;
//    /*4*/ sp.Prol_TymTr_=0.5;




//    Cell_simulator cell(sp, E);
//    SimParameters perturbedPar=sp.applyParameters(sp.getRandomParameters(0));
//    Experiment simulExp=cell.Simulate(perturbedPar ,E);

//    std::cout<<E.Result_i(0);

//     //Modificar num iteracines
//    //simulExp: simulado  E:experimental
//    //OptimizationResults O=cell.Optimize(sp,sp,simulExp,1,500);

//     OptimizationResults O=cell.Optimize(sp,sp,E,1,10);




//    /*std::cout<<"Lower to 1\n";
//    O=cell.Optimize(O.OptimalParameters(),E,1,200);
//    std::cout<<"Lower to 0.5\n";

//    O=cell.Optimize(O.OptimalParameters(),E,0.5,200);
//*/


//    std::ofstream f;
//    f.open("resultssim.txt");
//    f<<O;

//    cell.applyParameters(O.OptimalParameters(),E.Treatment_i(0));
//    cell.run();
//    cell.applyParameters(O.OptimalParameters(),E.Treatment_i(1));
//    cell.run();
//    cell.applyParameters(O.OptimalParameters(),E.Treatment_i(2));
//    cell.run();
//    f.close();
//}


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


void compareTimeStep()
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


    Cell_simulator cell(sp,sp, E);

    std::string filename="comparetimestep.txt";
    std::ofstream f;
    f.open(filename.c_str(),std::ios_base::app);

    cell.run(f,sp);
    f.close();



     //Modificar num iteracines
    //simulExp: simulado  E:experimental
    //OptimizationResults O=cell.Optimize(sp,sp,simulExp,1,500);

    cell.Optimize(sp,E,"testOptimization.txt");






}








void BayesParameters()
{   Treatment  med;
    med.Ag=0.0;
    med.Ab=0.0;
    med.sim_duration_d=120;
    med.time_step_d=1.0/120;
    med.init_cells=1.0e6;
    med.t_apop_meas_d=119;

    Treatment Mtb;
    Mtb.Ag=10.0;
    Mtb.Ab=0.0;
    Mtb.sim_duration_d=120;
    Mtb.time_step_d=1.0/120;
    Mtb.init_cells=1.0e6;
    Mtb.t_apop_meas_d=119;

    Treatment block;
    block.Ag=10.0;
    block.Ab=10.0;
    block.sim_duration_d=120;
    block.time_step_d=1.0/120;
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


    Cell_simulator cell(sp,sp, E);

     //Modificar num iteracines
    //simulExp: simulado  E:experimental
    //OptimizationResults O=cell.Optimize(sp,sp,simulExp,1,500);

    cell.Optimize(sp,E,"optimParam2.txt");




    /*std::cout<<"Lower to 1\n";
    O=cell.Optimize(O.OptimalParameters(),E,1,200);
    std::cout<<"Lower to 0.5\n";

    O=cell.Optimize(O.OptimalParameters(),E,0.5,200);
*/



}




void BayesParametersTest()
{   Treatment  med;
    med.Ag=0.0;
    med.Ab=0.0;
    med.sim_duration_d=120;
    med.time_step_d=1.0/120;
    med.init_cells=1.0e6;
    med.t_apop_meas_d=119;

    Treatment Mtb;
    Mtb.Ag=10.0;
    Mtb.Ab=0.0;
    Mtb.sim_duration_d=120;
    Mtb.time_step_d=1.0/120;
    Mtb.init_cells=1.0e6;
    Mtb.t_apop_meas_d=119;

    Treatment block;
    block.Ag=10.0;
    block.Ab=10.0;
    block.sim_duration_d=120;
    block.time_step_d=1.0/120;
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


    Cell_simulator cell(sp, sp,E);

    Experiment SimE=cell.Simulate(sp,E);


     //Modificar num iteracines
    //simulExp: simulado  E:experimental
    //OptimizationResults O=cell.Optimize(sp,sp,simulExp,1,500);

    cell.Optimize(sp,SimE,"testOptimization.txt");




    /*std::cout<<"Lower to 1\n";
    O=cell.Optimize(O.OptimalParameters(),E,1,200);
    std::cout<<"Lower to 0.5\n";

    O=cell.Optimize(O.OptimalParameters(),E,0.5,200);
*/



}





































































///*void experiment1()
//{
//    Cell_simulator simulation;
//    simulation.ask_parameters();
//    simulation.run();
//  }*/



//void experiment2()
//// Me dio SS=110
//{   Treatment  med;
//    med.Ag=0.0;
//    med.Ab=0.0;
//    med.sim_duration_d=120;
//    med.time_step_d=1.0/12;
//    med.init_cells=1e6;

//    Treatment Mtb;
//    Mtb.Ag=10.0;
//    Mtb.Ab=0.0;
//    Mtb.sim_duration_d=120;
//    Mtb.time_step_d=1.0/12;
//    Mtb.init_cells=1e6;

//    Treatment block;
//    block.Ag=10.0;
//    block.Ab=10.0;
//    block.sim_duration_d=120;
//    block.time_step_d=1.0/12;
//    block.init_cells=1e6;

//    Results MediaRes ("media");
//    Results MtbRes("mtb");
//    Results blockRes ("block");
//    Experiment E;
//    E.push_back(Mtb,MtbRes);
//    E.push_back(block,blockRes);
//    E.push_back(med,MediaRes);


//    SimParameters sp;
//    sp.mode_="FULL";
//    sp.max_num_cells_=2e6;
//    sp.init_ratio_APC_cells_=1e5/1e6;
//    sp.init_ratio_NK_cells_=1e5/1e6;
//    sp.init_ratio_LT_cells_=7.9e5/1e6;
//    sp.LT_ratio_specific_ =1000/1e6;

//    sp.APC_free_proliferation_rate_=1.0/120;
//    sp.APC_Ag_proliferation_rate_=1.0/120;
//    sp.APC_bound_proliferation_rate_=1.0/120;
//    sp.APC_blocked_proliferation_rate_=1.0/120;
//    sp.APC_exh_proliferation_rate_=1.0/120;
//    sp.NK_max_proliferation_rate_=1.0/120;
//    sp.LT_max_no_receptor_prol_rate_=1.0/480;
//    sp.LT_max_free_prol_rate_=0.00000000000000000001/1e11;
//    sp.LT_max_bound_prol_rate_=1.0/6;
//    sp.LT_max_blocked_prol_rate_=1.0/8;

//    sp.APC_free_apoptosis_rate_=1.0/240;
//    sp.APC_Ag_apoptosis_rate_=1.0/240;
//    sp.APC_bound_apoptosis_rate_=1.0/240;
//    sp.APC_blocked_apoptosis_rate_=1.0/240;
//    sp.APC_exh_apoptosis_rate_=1.0/240;

//    sp.APC_no_to_free_rate_per_Ag_=1.0/120;
//    sp.APC_free_to_bound_rate_per_LT_=1.0/36000;
//    sp.APC_Ab_binding_rate_=1.0/120;
//    sp.APC_exh_rate=1.0/6;
//    sp.NK_no_to_free_rate_per_Ag_=1.0/120;
//    sp.NK_free_to_bound_rate_per_LT_=1.0/36000;
//    sp.NK_Ab_binding_rate=1.0/120;
//    sp.NK_exh_rate=1.0/6;
//    sp.LT_no_to_free_rate_per_APC_=1.2/10000;
//    sp.LT_free_to_bound_rate_per_APC_=1.0/10;
//    sp.LT_mAb_binding_rate_=1.0/30;
//    sp.LT_exh_rate_=1.0/1e18;


//    sp.LT_IFN_no_rec_prod_rate_=0.001/1e5;
//    sp.LT_IFN_free_prod_rate_=0.000000000000000000001/10e11;
//    sp.LT_IFN_bound_prod_rate_=160.0/1e5;
//    sp.LT_IFN_blocked_prod_rate_=100.0/1e5;
//    sp.LT_TNF_no_rec_prod_rate_=0.0001/1e5;
//    sp.LT_TNF_free_prod_rate_=0.000000000000000000001/1e11;
//    sp.LT_TNF_bound_prod_rate_=6.0/1e7;
//    sp.LT_TNF_blocked_prod_rate_=4.0/1e7;
//    sp.APC_IFN_free_prod_rate_=0.00005/1e5;
//    sp.APC_IFN_Ag_prod_rate_=0.19/1e5;
//    sp.APC_IFN_bound_prod_rate_=0.39/1e5;
//    sp.APC_IFN_blocked_prod_rate_=0.19/1e5;
//    sp.APC_IFN_exh_prod_rate_=1.0/1e120;
//    sp.APC_TNF_free_prod_rate_=0.005/1.0e5;
//    sp.APC_TNF_Ag_prod_rate_=8.0/1e5;
//    sp.APC_TNF_bound_prod_rate_=6.0/1e5;
//    sp.APC_TNF_blocked_prod_rate_=8.0/1e5;
//    sp.APC_TNF_exh_prod_rate_=1.0/1e120;
//    sp.NK_IFN_free_prod_rate_=0.001/1e5;
//    sp.NK_IFN_Ag_prod_rate_=1.9/1e5;
//    sp.NK_IFN_blocked_prod_rate_=3.9/1e5;
//    sp.NK_IFN_bound_prod_rate_=1.9/1e5;
//    sp.NK_TNF_free_prod_rate_=0.05/1.0e5;
//    sp.NK_TNF_Ag_prod_rate_=0.008/1e5;
//    sp.NK_TNF_bound_prod_rate_=0.006/1e5;
//    sp.NK_TNF_blocked_prod_rate_=0.008/1e5;

//    sp.percentage_APC_IFN_free_prod_rate_=0.01;
//    sp.percentage_APC_IFN_Ag_prod_rate_=0.07;
//    sp.percentage_APC_IFN_bound_prod_rate_=0.11;
//    sp.percentage_APC_IFN_blocked_prod_rate_=0.07;
//    sp.percentage_APC_IFN_exh_prod_rate_=1.0/1e120;
//    sp.percentage_APC_TNF_free_prod_rate_=0.02;
//    sp.percentage_APC_TNF_Ag_prod_rate_=0.16;
//    sp.percentage_APC_TNF_bound_prod_rate_=0.41;
//    sp.percentage_APC_TNF_blocked_prod_rate_=0.16;
//    sp.percentage_APC_TNF_exh_prod_rate_=1.0/1e120;

//    //  sp.Ag_internalization_rate=0.5,
//    sp.TNF_deg=1.0/12;
//    sp.IFN_deg=1.0/12;




//    Cell_simulator cell(sp, E);
//    SimParameters perturbedPar=sp.applyParameters(sp.getRandomParameters(1));
//    Experiment simulExp=cell.Simulate(perturbedPar ,E);

//    std::cout<<E.Result_i(0);

//     //Modificar num iteracines
//    //simulExp: simulado  E:experimental
//    //OptimizationResults O=cell.Optimize(sp,sp,simulExp,1,500);

//     OptimizationResults O=cell.Optimize(sp,sp,E,1,10);




//    /*std::cout<<"Lower to 1\n";
//    O=cell.Optimize(O.OptimalParameters(),E,1,200);
//    std::cout<<"Lower to 0.5\n";

//    O=cell.Optimize(O.OptimalParameters(),E,0.5,200);
//*/


//    std::ofstream f;
//    f.open("resultssim.txt");
//    f<<O;

//    cell.applyParameters(O.OptimalParameters(),E.Treatment_i(0));
//    cell.run();
//    cell.applyParameters(O.OptimalParameters(),E.Treatment_i(1));
//    cell.run();
//    cell.applyParameters(O.OptimalParameters(),E.Treatment_i(2));
//    cell.run();
//    f.close();
//}


//void Fitted_Parameters_SS_110()
//{   Treatment  med;
//    med.Ag=0.0;
//    med.Ab=0.0;
//    med.sim_duration_d=120;
//    med.time_step_d=1.0/1200;
//    med.init_cells=1e6;

//    Treatment Mtb;
//    Mtb.Ag=10.0;
//    Mtb.Ab=0.0;
//    Mtb.sim_duration_d=120;
//    Mtb.time_step_d=1.0/1200;
//    Mtb.init_cells=1e6;

//    Treatment block;
//    block.Ag=10.0;
//    block.Ab=10.0;
//    block.sim_duration_d=120;
//    block.time_step_d=1.0/1200;
//    block.init_cells=1e6;

//    Results MediaRes ("media");
//    Results MtbRes("mtb");
//    Results blockRes ("block");
//    Experiment E;
//    E.push_back(Mtb,MtbRes);
//    E.push_back(block,blockRes);
//    E.push_back(med,MediaRes);


//    SimParameters sp;
//    sp.mode_="FULL";
//    sp.max_num_cells_=2.10511e+006;
//    sp.init_ratio_APC_cells_=0.0922583;
//    sp.init_ratio_NK_cells_=0.0970037;
//    sp.init_ratio_LT_cells_=0.759624;
//    sp.LT_ratio_specific_ =0.00567458;


//    sp.NK_max_proliferation_rate_=0.105814;
//    sp.LT_max_no_receptor_prol_rate_=0.003574690;
//    sp.LT_max_free_prol_rate_=6.46014e-010;
//    sp.LT_max_bound_prol_rate_=0.223669;
//    sp.LT_max_blocked_prol_rate_=0.117619;

//    sp.APC_no_to_free_rate_per_Ag_=0.0353429;
//    sp.APC_free_to_bound_rate_per_LT_=1.99072e-005;
//    sp.APC_Ab_binding_rate_=0.00556831;
//    sp.APC_exh_rate=0.180911;
//    sp.NK_no_to_free_rate_per_Ag_=0.0146091;
//    sp.NK_free_to_bound_rate_per_LT_=4.45254e-005;
//    sp.NK_Ab_binding_rate=0.0289378;
//    sp.NK_exh_rate=0.159299;
//    sp.LT_no_to_free_rate_per_APC_=6.44566e-005;
//    sp.LT_free_to_bound_rate_per_APC_=0.517016;
//    sp.LT_mAb_binding_rate_=0.0698594;
//    sp.LT_exh_rate_=2.3708e-010;



//    sp.APC_IFN_free_prod_rate_=2.3708e-010;
//    sp.APC_IFN_Ag_prod_rate_=3.68678e-006;
//    sp.APC_IFN_bound_prod_rate_=7.34169e-006;
//    sp.APC_IFN_blocked_prod_rate_=4.16568e-006;
//    sp.NK_IFN_free_prod_rate_=1.2142e-008;
//    sp.NK_IFN_Ag_prod_rate_=8.36737e-006;
//    sp.NK_IFN_bound_prod_rate_=0.000146465;
//    sp.NK_IFN_blocked_prod_rate_=9.93579e-005;

//    sp.LT_IFN_no_rec_prod_rate_=1.87039e-008;
//    sp.LT_IFN_free_prod_rate_=4.96099e-034;
//    sp.LT_IFN_bound_prod_rate_=0.000251217;
//    sp.LT_IFN_blocked_prod_rate_=8.79267e-005;


//    sp.APC_TNF_free_prod_rate_=3.14539e-008;
//    sp.APC_TNF_Ag_prod_rate_=0.000254004;
//    sp.APC_TNF_bound_prod_rate_=3.3537e-005;
//    sp.APC_TNF_blocked_prod_rate_=0.000285036;
//    sp.NK_TNF_free_prod_rate_=2.70394e-008;
//    sp.NK_TNF_Ag_prod_rate_=5.5439e-008;
//    sp.NK_TNF_bound_prod_rate_=4.99242e-008;
//    sp.NK_TNF_blocked_prod_rate_=6.80975e-008;
//    sp.LT_TNF_no_rec_prod_rate_=1.78249e-009;
//    sp.LT_TNF_free_prod_rate_=2.47759e-032;
//    sp.LT_TNF_bound_prod_rate_=3.42161e-006;
//    sp.LT_TNF_blocked_prod_rate_=0.055471e-006;

//    //  sp.Ag_internalization_rate=0.5,
//    sp.TNF_deg=0.053476;
//    sp.IFN_deg=0.0013241;




//    Cell_simulator cell(sp, E);
//    SimParameters perturbedPar=sp.applyParameters(sp.getRandomParameters(1));
//    Experiment simulExp=cell.Simulate(perturbedPar ,E);

//    std::cout<<E.Result_i(0);

//     //Modificar num iteracines
//    //simulExp: simulado  E:experimental
//    //OptimizationResults O=cell.Optimize(sp,sp,simulExp,1,500);

//     OptimizationResults O=cell.Optimize(sp,sp,E,1,150);




//    /*std::cout<<"Lower to 1\n";
//    O=cell.Optimize(O.OptimalParameters(),E,1,200);
//    std::cout<<"Lower to 0.5\n";

//    O=cell.Optimize(O.OptimalParameters(),E,0.5,200);
//*/


//    std::ofstream f;
//    f.open("resultssim.txt");
//    f<<O;

//    cell.applyParameters(O.OptimalParameters(),E.Treatment_i(0));
//    cell.run();
//    cell.applyParameters(O.OptimalParameters(),E.Treatment_i(1));
//    cell.run();
//    cell.applyParameters(O.OptimalParameters(),E.Treatment_i(2));
//    cell.run();
//    f.close();
//}

//void experiment3()
//{   Treatment  med;
//    med.Ag=0.0;
//    med.Ab=0.0;
//    med.sim_duration_d=120;
//    med.time_step_d=1.0/120;
//    med.init_cells=1e6;

//    Treatment Mtb;
//    Mtb.Ag=10.0;
//    Mtb.Ab=0.0;
//    Mtb.sim_duration_d=120;
//    Mtb.time_step_d=1.0/120;
//    Mtb.init_cells=1e6;

//    Treatment block;
//    block.Ag=10.0;
//    block.Ab=10.0;
//    block.sim_duration_d=120;
//    block.time_step_d=1.0/120;
//    block.init_cells=1e6;

//    Results MediaRes ("media");
//    Results MtbRes("mtb");
//    Results blockRes ("block");
//    Experiment E;
//    E.push_back(Mtb,MtbRes);
//    E.push_back(block,blockRes);
//    E.push_back(med,MediaRes);


//    SimParameters sp;
//    sp.mode_="PARTIAL";
//    sp.max_num_cells_=2e6;
//    sp.init_ratio_APC_cells_=1e5/1e6;
//    sp.init_ratio_NK_cells_=1e5/1e6;
//    sp.init_ratio_LT_cells_=7.9e5/1e6;
//    sp.LT_ratio_specific_ =1000/1e6;

//    sp.NK_max_proliferation_rate_=1.0/120;
//    sp.LT_max_no_receptor_prol_rate_=1.0/480;
//    sp.LT_max_free_prol_rate_=0.00000000000000000001/1e-11;
//    sp.LT_max_bound_prol_rate_=1.0/6;
//    sp.LT_max_blocked_prol_rate_=1.0/8;
//    sp.APC_no_to_free_rate_per_Ag_=1.0/120;
//    sp.APC_free_to_bound_rate_per_LT_=1.0/36000;
//    sp.APC_Ab_binding_rate_=1.0/120;
//    sp.NK_no_to_free_rate_per_Ag_=1.0/120;
//    sp.NK_free_to_bound_rate_per_LT_=1.0/36000;
//    sp.NK_Ab_binding_rate=1.0/120;
//    sp.LT_no_to_free_rate_per_APC_=1.2/10000;
//    sp.LT_free_to_bound_rate_per_APC_=1.0/1e2;
//    sp.LT_mAb_binding_rate_=1.0/30;
//    sp.LT_exh_rate_=1.0/1e4;
//    sp.APC_exh_rate=1.0/6;
//    sp.NK_exh_rate=1.0/6;
//    sp.LT_IFN_no_rec_prod_rate_=0.001/1e5;
//    sp.LT_IFN_free_prod_rate_=0.000000000000000000001/10e11;
//    sp.LT_IFN_bound_prod_rate_=160.0/1e5;
//    sp.LT_IFN_blocked_prod_rate_=100.0/1e5;
//    sp.LT_TNF_no_rec_prod_rate_=0.0001/1e5;
//    sp.LT_TNF_free_prod_rate_=0.000000000000000000001/1e11;
//    sp.LT_TNF_bound_prod_rate_=6.0/1e7;
//    sp.LT_TNF_blocked_prod_rate_=4.0/1e7;
//    sp.APC_IFN_free_prod_rate_=0.00005/1e5;
//    sp.APC_IFN_Ag_prod_rate_=0.19/1e5;
//    sp.APC_IFN_bound_prod_rate_=0.39/1e5;
//    //sp.APC_IFN_blocked_prod_rate_=0.19/1e5;
//    sp.APC_TNF_free_prod_rate_=0.005/1.0e5;
//    sp.APC_TNF_Ag_prod_rate_=8.0/1e5;
//    sp.APC_TNF_bound_prod_rate_=6.0/1e5;
//    //sp.APC_TNF_blocked_prod_rate_=8.0/1e5;
//    sp.NK_IFN_free_prod_rate_=0.001/1e5;
//    sp.NK_IFN_Ag_prod_rate_=1.9/1e5;
//    //sp.NK_IFN_blocked_prod_rate_=3.9/1e5;
//    sp.NK_IFN_bound_prod_rate_=1.9/1e5;
//    sp.NK_TNF_free_prod_rate_=0.05/1.0e5;
//    sp.NK_TNF_Ag_prod_rate_=0.008/1e5;
//    sp.NK_TNF_bound_prod_rate_=0.006/1e5;
//    //sp.NK_TNF_blocked_prod_rate_=0.008/1e5;
//    sp.max_num_cells_=2e6;
//    //  sp.Ag_internalization_rate=0.5,
//    sp.TNF_deg=1.0/12;
//    sp.IFN_deg=1.0/12;




//    Cell_simulator cell(sp, E);
//    SimParameters perturbedPar=sp.applyParameters(sp.getRandomParameters(1));
//    Experiment simulExp=cell.Simulate(perturbedPar ,E);

//    std::cout<<E.Result_i(0);

//     //Modificar num iteracines
//    //simulExp: simulado  E:experimental
//    //OptimizationResults O=cell.Optimize(sp,sp,simulExp,1,500);

//     OptimizationResults O=cell.Optimize(sp,sp,E,1,10);




//    /*std::cout<<"Lower to 1\n";
//    O=cell.Optimize(O.OptimalParameters(),E,1,200);
//    std::cout<<"Lower to 0.5\n";

//    O=cell.Optimize(O.OptimalParameters(),E,0.5,200);
//*/


//    std::ofstream f;
//    f.open("resultssim.txt");
//    f<<O;

//    cell.applyParameters(O.OptimalParameters(),E.Treatment_i(0));
//    cell.run();
//    cell.applyParameters(O.OptimalParameters(),E.Treatment_i(1));
//    cell.run();
//    cell.applyParameters(O.OptimalParameters(),E.Treatment_i(2));
//    cell.run();
//    f.close();
//}
