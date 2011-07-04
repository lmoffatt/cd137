#include "Includes/experiment1.h"
#include "Includes/Cell_simulator.h"


void experiment1()
{
    Cell_simulator simulation;
    simulation.ask_parameters();
    simulation.run();
  };

void experiment2()
{
    SimParameters sp;
    sp.Ag=10;
    sp.max_num_cells=2e6;
    sp.init_num_LT_cells=9e5;
    sp.LT_num_specific=1000;
    sp.init_num_APC_cells=1e5;
    sp.sim_duration_d=120;
    sp.time_step_d=0.01;
    sp.APC_max_proliferation_rate_=1.0/240;
    sp.LT_max_no_receptor_prol_rate_=1.0/96;
    sp.LT_max_free_prol_rate_=1.0/2;
    sp.LT_max_bound_prol_rate_=1.0/1.2;
    sp.APC_no_to_free_rate_per_Ag_=1.0/30;
    sp.APC_free_to_bound_rate_per_LT_=1.0/1e5;
    sp.APC_Ab_binding_rate_=9/1e5;
    sp.LT_no_to_free_rate_per_APC_=1.0/6e5;
    sp.LT_free_to_bound_rate_per_APC_=1.0/1e5;
    sp.APC_exh_rate=1/1e5;
    sp.LT_IFN_no_rec_prod_rate_=0.001/1e5;
    sp.LT_IFN_free_prod_rate_=101.0/1e5;
    sp.LT_IFN_bound_prod_rate_=200.0/1e5;
    sp.LT_TNF_no_rec_prod_rate_=0.001/1e5;
    sp.LT_TNF_free_prod_rate_=10.0/1e5;
    sp.LT_TNF_bound_prod_rate_=20.0/1e5;
    sp.APC_IFN_free_prod_rate_=0.5/1e5;
    sp.APC_IFN_Ag_prod_rate_=5.0/1e5;
    sp.APC_IFN_bound_prod_rate_=10.0/1e5;
    sp.APC_TNF_free_prod_rate_=5/1e5;
    sp.APC_TNF_Ag_prod_rate_=570/1e5;
    sp.APC_TNF_bound_prod_rate_=1110/1e5;
    Cell_simulator simulation(sp);
    simulation.run();
};

void experiment3()
{
    SimParameters sp;
    sp.Ag=10;
    sp.max_num_cells=2e6;
    sp.init_num_LT_cells=9e5;
    sp.LT_num_specific=1000;
    sp.init_num_APC_cells=1e5;
    sp.sim_duration_d=120;
    sp.time_step_d=0.01;
    sp.APC_max_proliferation_rate_=0;
    sp.LT_max_no_receptor_prol_rate_=1.0/96;
    sp.LT_max_free_prol_rate_=1.0/2;
    sp.LT_max_bound_prol_rate_=1.0/1.2;
    sp.APC_no_to_free_rate_per_Ag_=1.0/30;
    sp.APC_free_to_bound_rate_per_LT_=1.0/1e5;
    sp.APC_Ab_binding_rate_=9/1e5;
    sp.LT_no_to_free_rate_per_APC_=1.0/6e5;
    sp.LT_free_to_bound_rate_per_APC_=1.0/1e5;
    sp.APC_exh_rate=1/1e5;
    sp.LT_IFN_no_rec_prod_rate_=0.001/1e5;
    sp.LT_IFN_free_prod_rate_=101.0/1e5;
    sp.LT_IFN_bound_prod_rate_=200.0/1e5;
    sp.LT_TNF_no_rec_prod_rate_=0.001/1e5;
    sp.LT_TNF_free_prod_rate_=10.0/1e5;
    sp.LT_TNF_bound_prod_rate_=20.0/1e5;
    sp.APC_IFN_free_prod_rate_=0.5/1e5;
    sp.APC_IFN_Ag_prod_rate_=5.0/1e5;
    sp.APC_IFN_bound_prod_rate_=10.0/1e5;
    sp.APC_TNF_free_prod_rate_=5/1e5;
    sp.APC_TNF_Ag_prod_rate_=570/1e5;
    sp.APC_TNF_bound_prod_rate_=1110/1e5;
    //sp.Ag_internalization_rate=0.000005;
    Cell_simulator simulation(sp);
    simulation.run();
};


void experiment4()
{
    SimParameters sp;
    sp.Ag=10;
    sp.max_num_cells=2e6;
    sp.init_num_LT_cells=9e5;
    sp.LT_num_specific=1000;
    sp.init_num_APC_cells=1e5;
    sp.sim_duration_d=120;
    sp.time_step_d=0.01;
    sp.APC_max_proliferation_rate_=0;
    sp.LT_max_no_receptor_prol_rate_=1.0/96;
    sp.LT_max_free_prol_rate_=1.0/2;
    sp.LT_max_bound_prol_rate_=1.0/1.2;
    sp.APC_no_to_free_rate_per_Ag_=1.0/30;
    sp.APC_free_to_bound_rate_per_LT_=1.0/1e5;
    sp.APC_Ab_binding_rate_=9/1e5;
    sp.LT_no_to_free_rate_per_APC_=1.0/6e5;
    sp.LT_free_to_bound_rate_per_APC_=1.0/1e5;
    sp.APC_exh_rate=1/1e5;
    sp.LT_IFN_no_rec_prod_rate_=0.001/1e5;
    sp.LT_IFN_free_prod_rate_=25.0/1e5;
    sp.LT_IFN_bound_prod_rate_=100.0/1e5;
    sp.LT_TNF_no_rec_prod_rate_=0.001/1e5;
    sp.LT_TNF_free_prod_rate_=10.0/1e5;
    sp.LT_TNF_bound_prod_rate_=20.0/1e5;
    sp.APC_IFN_free_prod_rate_=0.5/1e5;
    sp.APC_IFN_Ag_prod_rate_=5.0/1e5;
    sp.APC_IFN_bound_prod_rate_=10.0/1e5;
    sp.APC_TNF_free_prod_rate_=5/1e5;
    sp.APC_TNF_Ag_prod_rate_=570/1e5;
    sp.APC_TNF_bound_prod_rate_=1110/1e5;
    //sp.Ag_internalization_rate=0.000005;
    sp.TNF_deg=0.5;
    Cell_simulator simulation(sp);
    simulation.run();
};

void experiment5()
{   SimParameters sp;
    sp.Ag=10,
    sp.max_num_cells=2e6,
    sp.init_num_LT_cells=9e5,
    sp.init_num_NK_cells=1e5,
    sp.LT_num_specific=1000,
    sp.init_num_APC_cells=1e5,
    sp.sim_duration_d=120,
    sp.time_step_d=0.01,
    sp.APC_max_proliferation_rate_=1.0/240,
    sp.NK_max_proliferation_rate_=1.0/240,
    sp.LT_max_no_receptor_prol_rate_=1.0/96,
    sp.LT_max_free_prol_rate_=1.0/2,
    sp.LT_max_bound_prol_rate_=1.0/1.2,
    sp.APC_no_to_free_rate_per_Ag_=1.0/30,
    sp.APC_free_to_bound_rate_per_LT_=1.0/1e5,
    sp.NK_no_to_free_rate_per_Ag_=1.0/30,
    sp.NK_free_to_bound_rate_per_LT_=1.0/1e5,
    sp.LT_no_to_free_rate_per_APC_=1.0/6e5,
    sp.LT_free_to_bound_rate_per_APC_=1.0/1e5,
    sp.APC_Ab_binding_rate_=9/1e5;
    sp.APC_exh_rate=1/1e5,
    sp.NK_exh_rate=1/1e5,
    sp.LT_IFN_no_rec_prod_rate_=0.001/1e5,
    sp.LT_IFN_free_prod_rate_=101.0/1e5,
    sp.LT_IFN_bound_prod_rate_=200.0/1e5,
    sp.LT_TNF_no_rec_prod_rate_=0.001/1e5,
    sp.LT_TNF_free_prod_rate_=10.0/1e5,
    sp.LT_TNF_bound_prod_rate_=20.0/1e5,
    sp.APC_IFN_free_prod_rate_=0.5/1e5,
    sp.APC_IFN_Ag_prod_rate_=5.0/1e5,
    sp.APC_IFN_bound_prod_rate_=10.0/1e5,
    sp.APC_TNF_free_prod_rate_=5/1e5,
    sp.APC_TNF_Ag_prod_rate_=570/1e5,
    sp.APC_TNF_bound_prod_rate_=1110/1e5,
    sp.NK_IFN_free_prod_rate_=0.5/1e5,
    sp.NK_IFN_Ag_prod_rate_=5.0/1e5,
    sp.NK_IFN_bound_prod_rate_=10.0/1e5,
    sp.NK_TNF_free_prod_rate_=5/1e5,
    sp.NK_TNF_Ag_prod_rate_=570/1e5,
    sp.NK_TNF_bound_prod_rate_=1110/1e5,
    sp.max_num_cells_=2e6,
    //sp.Ag_internalization_rate=0.5,
    sp.TNF_deg=10;
    Cell_simulator simulation(sp);
    simulation.run();
};

void experiment6()
{   SimParameters sp;
    sp.Ag=10,
    sp.Ab=5,
    sp.max_num_cells=2e6,
    sp.init_num_LT_cells=9e5,
    sp.init_num_NK_cells=1e5,
    sp.LT_num_specific=1000,
    sp.init_num_APC_cells=1e5,
    sp.sim_duration_d=120,
    sp.time_step_d=0.01,
    sp.APC_max_proliferation_rate_=1.0/240,
    sp.NK_max_proliferation_rate_=1.0/240,
    sp.LT_max_no_receptor_prol_rate_=1.0/96,
    sp.LT_max_free_prol_rate_=1.0/2,
    sp.LT_max_bound_prol_rate_=1.0/1.2,
    sp.APC_no_to_free_rate_per_Ag_=1.0/30,
    sp.APC_free_to_bound_rate_per_LT_=1.0/1e5,
    sp.APC_Ab_binding_rate_=9/1e5;
    sp.NK_no_to_free_rate_per_Ag_=1.0/30,
    sp.NK_free_to_bound_rate_per_LT_=1.0/1e5,
    sp.NK_Ab_binding_rate=1/5,
    sp.LT_no_to_free_rate_per_APC_=1.0/6e5,
    sp.LT_free_to_bound_rate_per_APC_=1.0/1e5,
    sp.APC_exh_rate=4200/1e5,
    sp.NK_exh_rate=4200/1e5,
    sp.LT_IFN_no_rec_prod_rate_=0.001/1e5,
    sp.LT_IFN_free_prod_rate_=30.0/1e5,
    sp.LT_IFN_bound_prod_rate_=70.0/1e5,
    sp.LT_TNF_no_rec_prod_rate_=0.001/1e5,
    sp.LT_TNF_free_prod_rate_=10.0/1e5,
    sp.LT_TNF_bound_prod_rate_=20.0/1e5,
    sp.APC_IFN_free_prod_rate_=0.5/1e5,
    sp.APC_IFN_Ag_prod_rate_=5.0/1e5,
    sp.APC_IFN_bound_prod_rate_=10.0/1e5,
    sp.APC_TNF_free_prod_rate_=5/1e5,
    sp.APC_TNF_Ag_prod_rate_=570/1e5,
    sp.APC_TNF_bound_prod_rate_=1110/1e5,
    sp.NK_IFN_free_prod_rate_=0.5/1e5,
    sp.NK_IFN_Ag_prod_rate_=15.0/1e5,
    sp.NK_IFN_blocked_prod_rate_=0.5/1e5,
    sp.NK_IFN_bound_prod_rate_=30.0/1e5,
    sp.NK_TNF_free_prod_rate_=5/1e5,
    sp.NK_TNF_Ag_prod_rate_=570/1e5,
    sp.NK_TNF_bound_prod_rate_=1110/1e5,
    sp.NK_TNF_blocked_prod_rate_=5/1e5,
    sp.max_num_cells_=2e6,
  //  sp.Ag_internalization_rate=0.5,
    sp.TNF_deg=0.015;
    Cell_simulator simulation(sp);
    simulation.run();
};
