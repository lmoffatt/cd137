#include <fstream>
#include "Includes/experiment1.h"
#include "Includes/Cell_simulator.h"
#include "Includes/Treatment.h"
#include "Includes/Experiment.h"
#include "Includes/OptimizationResults.h"

void con_APC()
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
    E.push_back(Mtb,MtbRes);
    E.push_back(block,blockRes);
    E.push_back(med,MediaRes);


    SimParameters sp;
    sp.mode_="PARTIAL";
    sp.max_num_cells_=1.99819e6;
    sp.init_ratio_APC_cells_=0.0867002;
    sp.init_ratio_NK_cells_=0.103059;
    sp.init_ratio_LT_cells_=0.819948;
    sp.LT_ratio_specific_ =0.00402011;

    sp.APC_free_proliferation_rate_=0.0134578;
    //sp.APC_Ag_proliferation_rate_=0.0134578;
    //sp.APC_bound_proliferation_rate_=0.0134578;
    //sp.APC_blocked_proliferation_rate_=0.0134578;
    //sp.APC_exh_proliferation_rate_=0.0134578;
    sp.NK_max_proliferation_rate_=0.106358;
    sp.LT_max_no_receptor_prol_rate_=0.106358;
    sp.LT_max_free_prol_rate_=2.68021e-010;
    sp.LT_max_bound_prol_rate_=0.290289;
    sp.LT_max_blocked_prol_rate_=0.124652;

   /* sp.APC_free_apoptosis_rate_=1.0/1e24;
    sp.APC_Ag_apoptosis_rate_=1.0/1e24;
    sp.APC_bound_apoptosis_rate_=1.0/1e24;
    sp.APC_blocked_apoptosis_rate_=1.0/1e24;
    sp.APC_exh_apoptosis_rate_=1.0/1e24;*/

    sp.APC_no_to_free_rate_per_Ag_=0.0242049;
    sp.APC_free_to_bound_rate_per_LT_=4.51511e-005;
    sp.APC_Ab_binding_rate_=0.0057502;
    sp.APC_exh_rate=0.187432;
    sp.NK_no_to_free_rate_per_Ag_=0.00763233;
    sp.NK_free_to_bound_rate_per_LT_=4.41789e-005;
    sp.NK_Ab_binding_rate=0.085456;
    sp.NK_exh_rate=0.225187;
    sp.LT_no_to_free_rate_per_APC_=9.96647e-005;
    sp.LT_free_to_bound_rate_per_APC_=0.00964449;
    sp.LT_mAb_binding_rate_=0.027425;
    sp.LT_exh_rate_=5.44769e-010;


    sp.LT_IFN_no_rec_prod_rate_=3.61199e-008;
    sp.LT_IFN_free_prod_rate_=3.91117e-034;
    sp.LT_IFN_bound_prod_rate_=0.000291428;
    sp.LT_IFN_blocked_prod_rate_=0.000291428;
    sp.LT_TNF_no_rec_prod_rate_=1.54684e-009;
    sp.LT_TNF_free_prod_rate_=3.63419e-032;
    sp.LT_TNF_bound_prod_rate_=3.37704e-006;
    sp.LT_TNF_blocked_prod_rate_=1.26275e-005;
    sp.APC_IFN_free_prod_rate_=1.81668e-010;
    sp.APC_IFN_Ag_prod_rate_=	9.05073e-006;
    sp.APC_IFN_bound_prod_rate_=3.49917e-006;
    sp.APC_IFN_blocked_prod_rate_=2.03066e-006;
    //sp.APC_IFN_exh_prod_rate_=1.0/1e120;
    sp.APC_TNF_free_prod_rate_=2.03171e-008;
    sp.APC_TNF_Ag_prod_rate_=	0.000495342;
    sp.APC_TNF_bound_prod_rate_=2.22827e-005;
    sp.APC_TNF_blocked_prod_rate_=0.00015285;
    //sp.APC_TNF_exh_prod_rate_=1.0/1e120;
    sp.NK_IFN_free_prod_rate_=3.18401e-008;
    sp.NK_IFN_Ag_prod_rate_=1.26801e-005;
    sp.NK_IFN_blocked_prod_rate_=0.000294432;
    sp.NK_IFN_bound_prod_rate_=0.000207201;
    sp.NK_TNF_free_prod_rate_=1.24834e-008;
    sp.NK_TNF_Ag_prod_rate_=3.52325e-008;
    sp.NK_TNF_bound_prod_rate_=3.52325e-008;
    sp.NK_TNF_blocked_prod_rate_=6.74271e-008;

    /*sp.percentage_APC_IFN_free_prod_rate_=0.01;
    sp.percentage_APC_IFN_Ag_prod_rate_=0.07;
    sp.percentage_APC_IFN_bound_prod_rate_=0.11;
    sp.percentage_APC_IFN_blocked_prod_rate_=0.07;
    sp.percentage_APC_IFN_exh_prod_rate_=1.0/1e24;
    sp.percentage_APC_TNF_free_prod_rate_=0.02;
    sp.percentage_APC_TNF_Ag_prod_rate_=0.16;
    sp.percentage_APC_TNF_bound_prod_rate_=0.41;
    sp.percentage_APC_TNF_blocked_prod_rate_=0.16;
    sp.percentage_APC_TNF_exh_prod_rate_=1.0/1e24;*/

    //  sp.Ag_internalization_rate=0.5,
    sp.TNF_deg=0.0113988;
    sp.IFN_deg=0.9067;




    Cell_simulator cell(sp, E);
    SimParameters perturbedPar=sp.applyParameters(sp.getRandomParameters(1));
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















































































/*void experiment1()
{
    Cell_simulator simulation;
    simulation.ask_parameters();
    simulation.run();
  }*/



void experiment2()
// Me dio SS=110
{   Treatment  med;
    med.Ag=0.0;
    med.Ab=0.0;
    med.sim_duration_d=120;
    med.time_step_d=1.0/12;
    med.init_cells=1e6;

    Treatment Mtb;
    Mtb.Ag=10.0;
    Mtb.Ab=0.0;
    Mtb.sim_duration_d=120;
    Mtb.time_step_d=1.0/12;
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
    E.push_back(Mtb,MtbRes);
    E.push_back(block,blockRes);
    E.push_back(med,MediaRes);


    SimParameters sp;
    sp.mode_="FULL";
    sp.max_num_cells_=2e6;
    sp.init_ratio_APC_cells_=1e5/1e6;
    sp.init_ratio_NK_cells_=1e5/1e6;
    sp.init_ratio_LT_cells_=7.9e5/1e6;
    sp.LT_ratio_specific_ =1000/1e6;

    sp.APC_free_proliferation_rate_=1.0/120;
    sp.APC_Ag_proliferation_rate_=1.0/120;
    sp.APC_bound_proliferation_rate_=1.0/120;
    sp.APC_blocked_proliferation_rate_=1.0/120;
    sp.APC_exh_proliferation_rate_=1.0/120;
    sp.NK_max_proliferation_rate_=1.0/120;
    sp.LT_max_no_receptor_prol_rate_=1.0/480;
    sp.LT_max_free_prol_rate_=0.00000000000000000001/1e11;
    sp.LT_max_bound_prol_rate_=1.0/6;
    sp.LT_max_blocked_prol_rate_=1.0/8;

    sp.APC_free_apoptosis_rate_=1.0/240;
    sp.APC_Ag_apoptosis_rate_=1.0/240;
    sp.APC_bound_apoptosis_rate_=1.0/240;
    sp.APC_blocked_apoptosis_rate_=1.0/240;
    sp.APC_exh_apoptosis_rate_=1.0/240;

    sp.APC_no_to_free_rate_per_Ag_=1.0/120;
    sp.APC_free_to_bound_rate_per_LT_=1.0/36000;
    sp.APC_Ab_binding_rate_=1.0/120;
    sp.APC_exh_rate=1.0/6;
    sp.NK_no_to_free_rate_per_Ag_=1.0/120;
    sp.NK_free_to_bound_rate_per_LT_=1.0/36000;
    sp.NK_Ab_binding_rate=1.0/120;
    sp.NK_exh_rate=1.0/6;
    sp.LT_no_to_free_rate_per_APC_=1.2/10000;
    sp.LT_free_to_bound_rate_per_APC_=1.0/10;
    sp.LT_mAb_binding_rate_=1.0/30;
    sp.LT_exh_rate_=1.0/1e18;


    sp.LT_IFN_no_rec_prod_rate_=0.001/1e5;
    sp.LT_IFN_free_prod_rate_=0.000000000000000000001/10e11;
    sp.LT_IFN_bound_prod_rate_=160.0/1e5;
    sp.LT_IFN_blocked_prod_rate_=100.0/1e5;
    sp.LT_TNF_no_rec_prod_rate_=0.0001/1e5;
    sp.LT_TNF_free_prod_rate_=0.000000000000000000001/1e11;
    sp.LT_TNF_bound_prod_rate_=6.0/1e7;
    sp.LT_TNF_blocked_prod_rate_=4.0/1e7;
    sp.APC_IFN_free_prod_rate_=0.00005/1e5;
    sp.APC_IFN_Ag_prod_rate_=0.19/1e5;
    sp.APC_IFN_bound_prod_rate_=0.39/1e5;
    sp.APC_IFN_blocked_prod_rate_=0.19/1e5;
    sp.APC_IFN_exh_prod_rate_=1.0/1e120;
    sp.APC_TNF_free_prod_rate_=0.005/1.0e5;
    sp.APC_TNF_Ag_prod_rate_=8.0/1e5;
    sp.APC_TNF_bound_prod_rate_=6.0/1e5;
    sp.APC_TNF_blocked_prod_rate_=8.0/1e5;
    sp.APC_TNF_exh_prod_rate_=1.0/1e120;
    sp.NK_IFN_free_prod_rate_=0.001/1e5;
    sp.NK_IFN_Ag_prod_rate_=1.9/1e5;
    sp.NK_IFN_blocked_prod_rate_=3.9/1e5;
    sp.NK_IFN_bound_prod_rate_=1.9/1e5;
    sp.NK_TNF_free_prod_rate_=0.05/1.0e5;
    sp.NK_TNF_Ag_prod_rate_=0.008/1e5;
    sp.NK_TNF_bound_prod_rate_=0.006/1e5;
    sp.NK_TNF_blocked_prod_rate_=0.008/1e5;

    sp.percentage_APC_IFN_free_prod_rate_=0.01;
    sp.percentage_APC_IFN_Ag_prod_rate_=0.07;
    sp.percentage_APC_IFN_bound_prod_rate_=0.11;
    sp.percentage_APC_IFN_blocked_prod_rate_=0.07;
    sp.percentage_APC_IFN_exh_prod_rate_=1.0/1e120;
    sp.percentage_APC_TNF_free_prod_rate_=0.02;
    sp.percentage_APC_TNF_Ag_prod_rate_=0.16;
    sp.percentage_APC_TNF_bound_prod_rate_=0.41;
    sp.percentage_APC_TNF_blocked_prod_rate_=0.16;
    sp.percentage_APC_TNF_exh_prod_rate_=1.0/1e120;

    //  sp.Ag_internalization_rate=0.5,
    sp.TNF_deg=1.0/12;
    sp.IFN_deg=1.0/12;




    Cell_simulator cell(sp, E);
    SimParameters perturbedPar=sp.applyParameters(sp.getRandomParameters(1));
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


void Fitted_Parameters_SS_110()
{   Treatment  med;
    med.Ag=0.0;
    med.Ab=0.0;
    med.sim_duration_d=120;
    med.time_step_d=1.0/1200;
    med.init_cells=1e6;

    Treatment Mtb;
    Mtb.Ag=10.0;
    Mtb.Ab=0.0;
    Mtb.sim_duration_d=120;
    Mtb.time_step_d=1.0/1200;
    Mtb.init_cells=1e6;

    Treatment block;
    block.Ag=10.0;
    block.Ab=10.0;
    block.sim_duration_d=120;
    block.time_step_d=1.0/1200;
    block.init_cells=1e6;

    Results MediaRes ("media");
    Results MtbRes("mtb");
    Results blockRes ("block");
    Experiment E;
    E.push_back(Mtb,MtbRes);
    E.push_back(block,blockRes);
    E.push_back(med,MediaRes);


    SimParameters sp;
    sp.mode_="FULL";
    sp.max_num_cells_=2.10511e+006;
    sp.init_ratio_APC_cells_=0.0922583;
    sp.init_ratio_NK_cells_=0.0970037;
    sp.init_ratio_LT_cells_=0.759624;
    sp.LT_ratio_specific_ =0.00567458;


    sp.NK_max_proliferation_rate_=0.105814;
    sp.LT_max_no_receptor_prol_rate_=0.003574690;
    sp.LT_max_free_prol_rate_=6.46014e-010;
    sp.LT_max_bound_prol_rate_=0.223669;
    sp.LT_max_blocked_prol_rate_=0.117619;

    sp.APC_no_to_free_rate_per_Ag_=0.0353429;
    sp.APC_free_to_bound_rate_per_LT_=1.99072e-005;
    sp.APC_Ab_binding_rate_=0.00556831;
    sp.APC_exh_rate=0.180911;
    sp.NK_no_to_free_rate_per_Ag_=0.0146091;
    sp.NK_free_to_bound_rate_per_LT_=4.45254e-005;
    sp.NK_Ab_binding_rate=0.0289378;
    sp.NK_exh_rate=0.159299;
    sp.LT_no_to_free_rate_per_APC_=6.44566e-005;
    sp.LT_free_to_bound_rate_per_APC_=0.517016;
    sp.LT_mAb_binding_rate_=0.0698594;
    sp.LT_exh_rate_=2.3708e-010;



    sp.APC_IFN_free_prod_rate_=2.3708e-010;
    sp.APC_IFN_Ag_prod_rate_=3.68678e-006;
    sp.APC_IFN_bound_prod_rate_=7.34169e-006;
    sp.APC_IFN_blocked_prod_rate_=4.16568e-006;
    sp.NK_IFN_free_prod_rate_=1.2142e-008;
    sp.NK_IFN_Ag_prod_rate_=8.36737e-006;
    sp.NK_IFN_bound_prod_rate_=0.000146465;
    sp.NK_IFN_blocked_prod_rate_=9.93579e-005;

    sp.LT_IFN_no_rec_prod_rate_=1.87039e-008;
    sp.LT_IFN_free_prod_rate_=4.96099e-034;
    sp.LT_IFN_bound_prod_rate_=0.000251217;
    sp.LT_IFN_blocked_prod_rate_=8.79267e-005;


    sp.APC_TNF_free_prod_rate_=3.14539e-008;
    sp.APC_TNF_Ag_prod_rate_=0.000254004;
    sp.APC_TNF_bound_prod_rate_=3.3537e-005;
    sp.APC_TNF_blocked_prod_rate_=0.000285036;
    sp.NK_TNF_free_prod_rate_=2.70394e-008;
    sp.NK_TNF_Ag_prod_rate_=5.5439e-008;
    sp.NK_TNF_bound_prod_rate_=4.99242e-008;
    sp.NK_TNF_blocked_prod_rate_=6.80975e-008;
    sp.LT_TNF_no_rec_prod_rate_=1.78249e-009;
    sp.LT_TNF_free_prod_rate_=2.47759e-032;
    sp.LT_TNF_bound_prod_rate_=3.42161e-006;
    sp.LT_TNF_blocked_prod_rate_=0.055471e-006;

    //  sp.Ag_internalization_rate=0.5,
    sp.TNF_deg=0.053476;
    sp.IFN_deg=0.0013241;




    Cell_simulator cell(sp, E);
    SimParameters perturbedPar=sp.applyParameters(sp.getRandomParameters(1));
    Experiment simulExp=cell.Simulate(perturbedPar ,E);

    std::cout<<E.Result_i(0);

     //Modificar num iteracines
    //simulExp: simulado  E:experimental
    //OptimizationResults O=cell.Optimize(sp,sp,simulExp,1,500);

     OptimizationResults O=cell.Optimize(sp,sp,E,1,150);




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

void experiment3()
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
    block.time_step_d=1.0/120;
    block.init_cells=1e6;

    Results MediaRes ("media");
    Results MtbRes("mtb");
    Results blockRes ("block");
    Experiment E;
    E.push_back(Mtb,MtbRes);
    E.push_back(block,blockRes);
    E.push_back(med,MediaRes);


    SimParameters sp;
    sp.mode_="PARTIAL";
    sp.max_num_cells_=2e6;
    sp.init_ratio_APC_cells_=1e5/1e6;
    sp.init_ratio_NK_cells_=1e5/1e6;
    sp.init_ratio_LT_cells_=7.9e5/1e6;
    sp.LT_ratio_specific_ =1000/1e6;

    sp.NK_max_proliferation_rate_=1.0/120;
    sp.LT_max_no_receptor_prol_rate_=1.0/480;
    sp.LT_max_free_prol_rate_=0.00000000000000000001/1e-11;
    sp.LT_max_bound_prol_rate_=1.0/6;
    sp.LT_max_blocked_prol_rate_=1.0/8;
    sp.APC_no_to_free_rate_per_Ag_=1.0/120;
    sp.APC_free_to_bound_rate_per_LT_=1.0/36000;
    sp.APC_Ab_binding_rate_=1.0/120;
    sp.NK_no_to_free_rate_per_Ag_=1.0/120;
    sp.NK_free_to_bound_rate_per_LT_=1.0/36000;
    sp.NK_Ab_binding_rate=1.0/120;
    sp.LT_no_to_free_rate_per_APC_=1.2/10000;
    sp.LT_free_to_bound_rate_per_APC_=1.0/1e2;
    sp.LT_mAb_binding_rate_=1.0/30;
    sp.LT_exh_rate_=1.0/1e4;
    sp.APC_exh_rate=1.0/6;
    sp.NK_exh_rate=1.0/6;
    sp.LT_IFN_no_rec_prod_rate_=0.001/1e5;
    sp.LT_IFN_free_prod_rate_=0.000000000000000000001/10e11;
    sp.LT_IFN_bound_prod_rate_=160.0/1e5;
    sp.LT_IFN_blocked_prod_rate_=100.0/1e5;
    sp.LT_TNF_no_rec_prod_rate_=0.0001/1e5;
    sp.LT_TNF_free_prod_rate_=0.000000000000000000001/1e11;
    sp.LT_TNF_bound_prod_rate_=6.0/1e7;
    sp.LT_TNF_blocked_prod_rate_=4.0/1e7;
    sp.APC_IFN_free_prod_rate_=0.00005/1e5;
    sp.APC_IFN_Ag_prod_rate_=0.19/1e5;
    sp.APC_IFN_bound_prod_rate_=0.39/1e5;
    //sp.APC_IFN_blocked_prod_rate_=0.19/1e5;
    sp.APC_TNF_free_prod_rate_=0.005/1.0e5;
    sp.APC_TNF_Ag_prod_rate_=8.0/1e5;
    sp.APC_TNF_bound_prod_rate_=6.0/1e5;
    //sp.APC_TNF_blocked_prod_rate_=8.0/1e5;
    sp.NK_IFN_free_prod_rate_=0.001/1e5;
    sp.NK_IFN_Ag_prod_rate_=1.9/1e5;
    //sp.NK_IFN_blocked_prod_rate_=3.9/1e5;
    sp.NK_IFN_bound_prod_rate_=1.9/1e5;
    sp.NK_TNF_free_prod_rate_=0.05/1.0e5;
    sp.NK_TNF_Ag_prod_rate_=0.008/1e5;
    sp.NK_TNF_bound_prod_rate_=0.006/1e5;
    //sp.NK_TNF_blocked_prod_rate_=0.008/1e5;
    sp.max_num_cells_=2e6;
    //  sp.Ag_internalization_rate=0.5,
    sp.TNF_deg=1.0/12;
    sp.IFN_deg=1.0/12;




    Cell_simulator cell(sp, E);
    SimParameters perturbedPar=sp.applyParameters(sp.getRandomParameters(1));
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
