#include <fstream>
#include "Includes/experiment1.h"
#include "Includes/Cell_simulator.h"
#include "Includes/Treatment.h"
#include "Includes/Experiment.h"
#include "Includes/OptimizationResults.h"

void experiment1()
{
    Cell_simulator simulation;
    simulation.ask_parameters();
    simulation.run();
  }



void experiment2()
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
    sp.mode_="FULL";
    sp.max_num_cells_=2e6;
    sp.init_ratio_APC_cells_=1e5/1e6;
    sp.init_ratio_NK_cells_=1e5/1e6;
    sp.init_ratio_LT_cells_=7.9e5/1e6;
    sp.LT_ratio_specific_ =1000/1e6;
    sp.APC_max_proliferation_rate_=1.0/120;
    sp.NK_max_proliferation_rate_=1.0/120;
    sp.LT_max_no_receptor_prol_rate_=1.0/240;
    sp.LT_max_free_prol_rate_=1.0/8;
    sp.LT_max_bound_prol_rate_=1.0/6;
    sp.LT_max_blocked_prol_rate_=1.0/8;
    sp.APC_no_to_free_rate_per_Ag_=1.0/120;
    sp.APC_free_to_bound_rate_per_LT_=1.0/36000;
    sp.APC_Ab_binding_rate_=1.0/120;
    sp.NK_no_to_free_rate_per_Ag_=1.0/120;
    sp.NK_free_to_bound_rate_per_LT_=1.0/36000;
    sp.NK_Ab_binding_rate=1.0/120;
    sp.LT_no_to_free_rate_per_APC_=1.2/1.0e4;
    sp.LT_free_to_bound_rate_per_APC_=1.0/1e6;
    sp.LT_mAb_binding_rate_=1.0/30;
    sp.APC_exh_rate=1.0/6;
    sp.NK_exh_rate=1.0/6;
    sp.LT_IFN_no_rec_prod_rate_=0.0001/240e6;
    sp.LT_IFN_free_prod_rate_=1000.0/1e5;
    sp.LT_IFN_bound_prod_rate_=1600.0/1e5;
    sp.LT_IFN_blocked_prod_rate_=1000.0/1e5;
    sp.LT_TNF_no_rec_prod_rate_=0.0000001/240e6;
    sp.LT_TNF_free_prod_rate_=0.00001/240e6;
    sp.LT_TNF_bound_prod_rate_=0.00002/240e6;
    sp.LT_TNF_blocked_prod_rate_=0.00001/240e6;
    sp.APC_IFN_free_prod_rate_=0.00005/120e5;
    sp.APC_IFN_Ag_prod_rate_=0.000005/120e5;
    sp.APC_IFN_bound_prod_rate_=0.000001/120e5;
    sp.APC_IFN_blocked_prod_rate_=0.000005/120e5;
    sp.APC_TNF_free_prod_rate_=5.0/1.0e7;
    sp.APC_TNF_Ag_prod_rate_=8.0/100000;
    sp.APC_TNF_bound_prod_rate_=6.0/100000;
    sp.APC_TNF_blocked_prod_rate_=8.0/100000;
    sp.NK_IFN_free_prod_rate_=0.001/1e5;
    sp.NK_IFN_Ag_prod_rate_=1.9/1e5;
    sp.NK_IFN_blocked_prod_rate_=3.9/1e5;
    sp.NK_IFN_bound_prod_rate_=1.9/1e5;
    sp.NK_TNF_free_prod_rate_=0.00001/120e5;
    sp.NK_TNF_Ag_prod_rate_=0.00002/120e5;
    sp.NK_TNF_bound_prod_rate_=0.000015/120e5;
    sp.NK_TNF_blocked_prod_rate_=0.00002/120e5;
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


void fitpar()
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
    sp.APC_max_proliferation_rate_=0.00855905;
    sp.NK_max_proliferation_rate_=0.0413318;
    sp.LT_max_no_receptor_prol_rate_=0.00693542;
    sp.LT_max_free_prol_rate_=1.06285;
    sp.LT_max_bound_prol_rate_=0.678239;
    sp.LT_max_blocked_prol_rate_=1.06285;
    sp.APC_no_to_free_rate_per_Ag_=1.49128;
    sp.APC_free_to_bound_rate_per_LT_=2.55635e-5;
    sp.APC_Ab_binding_rate_=0.056328;
    sp.APC_exh_rate=0.156994;
    sp.NK_no_to_free_rate_per_Ag_=0.00585877;
    sp.NK_free_to_bound_rate_per_LT_=2.56117e-5;
    sp.NK_Ab_binding_rate=0.126604;
    sp.NK_exh_rate=0.820873;
    sp.LT_no_to_free_rate_per_APC_=0.000160848;
    sp.LT_free_to_bound_rate_per_APC_=1.21087e-6;
    sp.LT_mAb_binding_rate_=0.207377;
    sp.APC_IFN_free_prod_rate_=3.6471e-12;
    sp.APC_IFN_Ag_prod_rate_=	1.68612e-12;
    sp.APC_IFN_bound_prod_rate_=9.55028e-14;
    sp.APC_IFN_blocked_prod_rate_=1.68612e-12;
    sp.NK_IFN_free_prod_rate_=7.02324e-5;
    sp.NK_IFN_Ag_prod_rate_=2.15214e-8;
    sp.NK_IFN_bound_prod_rate_=4.11724e-10;
    sp.NK_IFN_blocked_prod_rate_=2.15214e-8;
    sp.LT_IFN_no_rec_prod_rate_=4.56079e-12;
    sp.LT_IFN_free_prod_rate_=2.43116e-9;
    sp.LT_IFN_bound_prod_rate_=1.17006e-7;
    sp.LT_IFN_blocked_prod_rate_=2.43116e-9;
    sp.APC_TNF_free_prod_rate_=1.31233e-6;
    sp.APC_TNF_Ag_prod_rate_=0.000157061;
    sp.APC_TNF_bound_prod_rate_=8.51921e-5;
    sp.APC_TNF_blocked_prod_rate_=0.000157061;
    sp.NK_TNF_free_prod_rate_=3.45682e-11;
    sp.NK_TNF_Ag_prod_rate_=7.84051e-12;
    sp.NK_TNF_bound_prod_rate_=2.41617e-11;
    sp.NK_TNF_blocked_prod_rate_=7.84051e-12;
    sp.LT_TNF_no_rec_prod_rate_=4.327e-15;
    sp.LT_TNF_free_prod_rate_=1.35635e-014;
    sp.LT_TNF_bound_prod_rate_=9.8869e-014;
    sp.LT_TNF_blocked_prod_rate_=1.35635e-014;





    sp.max_num_cells_=2e6;
    //  sp.Ag_internalization_rate=0.5,
    sp.TNF_deg=0.00763902;
    sp.IFN_deg=0.00763902;




    Cell_simulator cell(sp, E);
    SimParameters perturbedPar=sp.applyParameters(sp.getRandomParameters(1));
    Experiment simulExp=cell.Simulate(perturbedPar ,E);

    std::cout<<simulExp.Result_i(0);

     //Modificar num iteracines
    //simulExp: simulado  E:experimental
    //OptimizationResults O=cell.Optimize(sp,sp,simulExp,1,500);

     OptimizationResults O=cell.Optimize(sp,sp,E,0.01,5);
     cell.applyParameters(O.OptimalParameters(),E.Treatment_i(0));
     cell.run();

    /*std::cout<<"Lower to 1\n";
    O=cell.Optimize(O.OptimalParameters(),E,1,200);
    std::cout<<"Lower to 0.5\n";

    O=cell.Optimize(O.OptimalParameters(),E,0.5,200);
*/


    std::ofstream f;
    f.open("resultssim.txt");
    f<<O;
}
/*void experiment3()
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
    sp.APC_max_proliferation_rate_=1.0/120;
    sp.NK_max_proliferation_rate_=1.0/120;
    sp.LT_max_no_receptor_prol_rate_=1.0/1200;
    sp.LT_max_free_prol_rate_=1.0/4;
    sp.LT_max_bound_prol_rate_=1.0/2;
    sp.LT_max_blocked_prol_rate_=1.0/4;
    sp.APC_no_to_free_rate_per_Ag_=1.0/400;
    sp.APC_free_to_bound_rate_per_LT_=0.25/1e5;
    sp.APC_Ab_binding_rate_=1.0/30;
    sp.NK_no_to_free_rate_per_Ag_=1.0/600;
    sp.NK_free_to_bound_rate_per_LT_=0.25/1e5;
    sp.NK_Ab_binding_rate=1.0/30;
    sp.LT_no_to_free_rate_per_APC_=1.0/6e5;
    sp.LT_free_to_bound_rate_per_APC_=1.0/1e6;
    sp.LT_mAb_binding_rate_=1.0/30;
    sp.APC_exh_rate=1.0/25;
    sp.NK_exh_rate=1.0/25;
    sp.LT_IFN_no_rec_prod_rate_=0.0001/240e6;
    sp.LT_IFN_free_prod_rate_=3.0/240e6;
    sp.LT_IFN_bound_prod_rate_=9.0/240e6;
    sp.LT_IFN_blocked_prod_rate_=3.0/240e6;
    sp.LT_TNF_no_rec_prod_rate_=0.0000001/240e6;
    sp.LT_TNF_free_prod_rate_=0.00001/240e6;
    sp.LT_TNF_bound_prod_rate_=0.00002/240e6;
    sp.LT_TNF_blocked_prod_rate_=0.00001/240e6;
    sp.APC_IFN_free_prod_rate_=0.00005/120e5;
    sp.APC_IFN_Ag_prod_rate_=0.000005/120e5;
    sp.APC_IFN_bound_prod_rate_=0.000001/120e5;
    sp.APC_IFN_blocked_prod_rate_=0.000005/120e5;
    sp.APC_TNF_free_prod_rate_=5.0/120e5;
    sp.APC_TNF_Ag_prod_rate_=8000.0/120e5;
    sp.APC_TNF_bound_prod_rate_=6000.0/120e5;
    sp.APC_TNF_blocked_prod_rate_=8000.0/120e5;
    sp.NK_IFN_free_prod_rate_=0.0001/240e6;
    sp.NK_IFN_Ag_prod_rate_=0.3/240e6;
    sp.NK_IFN_blocked_prod_rate_=0.9/240e6;
    sp.NK_IFN_bound_prod_rate_=0.3/240e6;
    sp.NK_TNF_free_prod_rate_=0.00001/120e5;
    sp.NK_TNF_Ag_prod_rate_=0.00002/120e5;
    sp.NK_TNF_bound_prod_rate_=0.000015/120e5;
    sp.NK_TNF_blocked_prod_rate_=0.00002/120e5;
    sp.max_num_cells_=2e6;
    sp.Ag_internalization_rate=0.5,
    sp.TNF_deg=1.0/120;
    sp.IFN_deg=1.0/120;

    sp.APC_TNF_free_prod_rate_=5.0/1.0e7;
    sp.APC_TNF_Ag_prod_rate_=8.0/100000;
    sp.APC_TNF_bound_prod_rate_=6.0/100000;
    sp.APC_TNF_blocked_prod_rate_=8.0/100000;


    Cell_simulator cell(sp, E);
    SimParameters perturbedPar=sp.applyParameters(sp.getRandomParameters(0.01));
    Experiment simulExp=cell.Simulate(perturbedPar
                                      ,E);

    std::cout<<simulExp.Result_i(0);


    Modificar num iteracines
    OptimizationResults O=cell.Optimize(sp,sp,simulExp,0.01,500);

    OptimizationResults O=cell.Optimize(sp,E,4,500);

    std::cout<<"Lower to 1\n";
    O=cell.Optimize(O.OptimalParameters(),E,1,200);
    std::cout<<"Lower to 0.5\n";

    O=cell.Optimize(O.OptimalParameters(),E,0.5,200);



    std::ofstream f;
    f.open("resultssim.txt");
    f<<O;
}*/
/*void con_bloqueo ()
{   SimParameters sp;
    sp.Ag=10,
    sp.Ab=5,
    sp.max_num_cells=2e6,
    sp.init_ratio_LT_cells=9e5,
    sp.init_ratio_NK_cells=1e5,
    sp.LT_ratio_specific=1000,
    sp.init_ratio_APC_cells=1e5,
    sp.sim_duration_d=120,
    sp.time_step_d=1.0/2e5,
    sp.APC_max_proliferation_rate_=1.0/60,
    sp.NK_max_proliferation_rate_=1.0/60,
    sp.LT_max_no_receptor_prol_rate_=1.0/240,
    sp.LT_max_free_prol_rate_=1.0/20,
    sp.LT_max_bound_prol_rate_=1.0/10,
    sp.LT_max_blocked_prol_rate_=1.0/20,
    sp.APC_no_to_free_rate_per_Ag_=1.0/300,
    sp.APC_free_to_bound_rate_per_LT_=1.0/6,
    sp.APC_Ab_binding_rate_=1.0/30;
    sp.NK_no_to_free_rate_per_Ag_=1.0/400,
    sp.NK_free_to_bound_rate_per_LT_=1.0/6,
    sp.NK_Ab_binding_rate=1.0/6,
    sp.LT_no_to_free_rate_per_APC_=1.0/100,
    sp.LT_free_to_bound_rate_per_APC_=1.0/100,
    sp.LT_mAb_binding_rate_=1.0/10,
    sp.APC_exh_rate=1.0/10,
    sp.NK_exh_rate=1.0/12,
    sp.LT_IFN_no_rec_prod_rate_=0.0001/240e6,
    sp.LT_IFN_free_prod_rate_=300.0/240e7,
    sp.LT_IFN_bound_prod_rate_=900.0/240e7,
    sp.LT_IFN_blocked_prod_rate_=300.0/240e7,
    sp.LT_TNF_no_rec_prod_rate_=0.0000001/240e6,
    sp.LT_TNF_free_prod_rate_=0.00001/240e6,
    sp.LT_TNF_bound_prod_rate_=0.00002/240e6,
    sp.LT_TNF_blocked_prod_rate_=0.00001/240e6,
    sp.APC_IFN_free_prod_rate_=0.00005/120e5,
    sp.APC_IFN_Ag_prod_rate_=0.000005/120e5,
    sp.APC_IFN_bound_prod_rate_=0.000001/120e5,
    sp.APC_IFN_blocked_prod_rate_=0.000005/120e5,
    sp.APC_TNF_free_prod_rate_=5.0/120e5,
    sp.APC_TNF_Ag_prod_rate_=8000.0/120e5,
    sp.APC_TNF_bound_prod_rate_=6000.0/120e5,
    sp.APC_TNF_blocked_prod_rate_=8000.0/120e5,
    sp.NK_IFN_free_prod_rate_=0.0001/240e6,
    sp.NK_IFN_Ag_prod_rate_=0.9/1e4,
    sp.NK_IFN_blocked_prod_rate_=0.9/1e4,
    sp.NK_IFN_bound_prod_rate_=0.3/1e4,
    sp.NK_TNF_free_prod_rate_=0.00001/120e5,
    sp.NK_TNF_Ag_prod_rate_=0.00002/120e5,
    sp.NK_TNF_bound_prod_rate_=0.000015/120e5,
    sp.NK_TNF_blocked_prod_rate_=0.00002/120e5,
    sp.max_num_cells_=2e6,
  //  sp.Ag_internalization_rate=0.5,
    sp.TNF_deg=1.0/6;
    sp.IFN_deg=1.0/7;
    Cell_simulator simulation(sp);
    simulation.run();
    Results myRes("Fig_1_2");
    Results mySim=simulation.Simulate(sp,myRes);
    SumSquareTXT (myRes,mySim);

};
*/
