#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <cstddef>
#include <vector>
#include "Results.h"
#include "Includes/Cell_simulator.h"
#include "Includes/SimParameters.h"
#include "Includes/Media.h"
#include "Includes/APC.h"
#include "Includes/NK.h"
#include "Includes/LT.h"
#include "Includes/LevenbergMarquardt.h"
#include "Includes/OptimizationResults.h"


/*void Cell_simulator::ask_parameters()
{
    double max_num_cells_;
    double init_num_APC_cells;
    double init_num_NK_cells;
    double init_num_LT_cells;
    double LT_num_specific;
    double Ag;
    double Ab;
    //  double Ag_internalization_rate;
    double APC_free_proliferation_rate_;
    double APC_Ag_proliferation_rate_;
    double APC_bound_proliferation_rate_;
    double APC_blocked_proliferation_rate_;
    double APC_exh_proliferation_rate_;
    double APC_free_apop_rate_;
    double APC_Ag_apop_rate_;
    double APC_bound_apop_rate_;
    double APC_blocked_apop_rate_;
    double APC_exh_apop_rate_;
    double NK_max_proliferation_rate_;
    double LT_max_no_receptor_prol_rate_;
    double LT_max_free_prol_rate_;
    double LT_max_bound_prol_rate_;
    double LT_max_blocked_prol_rate_;
    double APC_no_to_free_rate_per_Ag_;
    double APC_free_to_bound_rate_per_LT_;
    double APC_Ab_binding_rate_;
    double APC_exh_rate;
    double NK_no_to_free_rate_per_Ag_;
    double NK_free_to_bound_rate_per_LT_;
    double NK_Ab_binding_rate;
    double NK_exh_rate;
    double LT_no_to_free_rate_per_APC_;
    double LT_free_to_bound_rate_per_APC_;
    double LT_mAb_binding_rate_;
    double LT_IFN_no_rec_prod_rate_;
    double LT_IFN_free_prod_rate_;
    double LT_IFN_bound_prod_rate_;
    double LT_IFN_blocked_prod_rate_;
    double LT_TNF_no_rec_prod_rate_;
    double LT_TNF_free_prod_rate_;
    double LT_TNF_bound_prod_rate_;
    double LT_TNF_blocked_prod_rate_;
    double APC_IFN_free_prod_rate_;
    double APC_IFN_Ag_prod_rate_;
    double APC_IFN_bound_prod_rate_;
    double APC_IFN_blocked_prod_rate_;
    double APC_IFN_exh_prod_rate_;
    double percentage_APC_IFN_free_prod_rate_;
    double percentage_APC_IFN_Ag_prod_rate_;
    double percentage_APC_IFN_bound_prod_rate_;
    double percentage_APC_IFN_blocked_prod_rate_;
    double percentage_APC_IFN_exh_prod_rate_;

    double APC_TNF_free_prod_rate_;
    double APC_TNF_Ag_prod_rate_;
    double APC_TNF_bound_prod_rate_;
    double APC_TNF_blocked_prod_rate_;
    double percentage_APC_TNF_free_prod_rate_;
    double percentage_APC_TNF_Ag_prod_rate_;
    double percentage_APC_TNF_bound_prod_rate_;
    double percentage_APC_TNF_blocked_prod_rate_;
    double percentage_APC_TNF_exh_prod_rate_;
    double NK_IFN_free_prod_rate_;
    double NK_IFN_Ag_prod_rate_;
    double NK_IFN_bound_prod_rate_;
    double NK_IFN_blocked_prod_rate_;
    double NK_TNF_free_prod_rate_;
    double NK_TNF_Ag_prod_rate_;
    double NK_TNF_bound_prod_rate_;
    double NK_TNF_blocked_prod_rate_;
    double TNF_deg;
    double IFN_deg;
    double LT_exh_rate_;

    ///Asking the user the parameters that define the system
    std::cout<<"Enter the filename where you want to store the values\n";
    std::cin>>this->filename;
    ///Media Parameters
    std::cout<<"enter the number of cells that the media supports\n";
    std::cout<<"default is 2e6, enter -1 to keep this value\n";
    std::cin>>max_num_cells_;
    if (max_num_cells_==-1)
        max_num_cells_=2e6;

    std::cout<<"enter the applied concentration of Antigen\n";
    std::cout<<"default is 10 ug/ml, enter -1 to keep this value\n";
    std::cin>>Ag;
    if (Ag==-1)
        Ag=5;

    std::cout<<"enter the applied concentration of blocking Ab\n";
    std::cout<<"default is 10 ug/ml, enter -1 to keep this value\n";
    std::cin>>Ab;
    if (Ab==-1)
        Ab=5;

      std::cout<<"Enter de Antigen internalization rate";
    std::cout<<"default is, 0.000005, enter -1 to keep this value";
    std::cin>>this->Ag_internalization_rate;
    if (Ag_internalization_rate==-1)
        Ag_internalization_rate==0.000005;

    std::cout<<"Enter the degradation rate of TNF";
    std::cin>>TNF_deg;

    /// Initial number of cells
    std::cout<<"enter the intial number of APC cells\n";
    std::cout<<"default is 1e5, enter -1 to keep this value\n";
    std::cin>>init_num_APC_cells;
    if (init_num_APC_cells==-1)
        init_num_APC_cells=1e5;

    std::cout<<"enter the intial number of NK cells\n";
    std::cout<<"default is 1e5, enter -1 to keep this value\n";
    std::cin>>init_num_NK_cells;
    if (init_num_NK_cells==-1)
        init_num_NK_cells=1e5;

    std::cout<<"enter the intial number of LT cells\n";
    std::cout<<"default is 9e5, enter -1 to keep this value\n";
    std::cin>>init_num_LT_cells;
    if (init_num_LT_cells==-1)
        init_num_LT_cells=9e5;

    std::cout<<"enter the intial number of LT cells specific for the antigen \n";
    std::cout<<"default is 1000, enter -1  to keep this value\n";
    std::cin>>LT_num_specific;
    if (LT_num_specific==-1)
        LT_num_specific=1000;


    ///Parameters of time
    std::cout<<"enter the duration of the simulation in hours\n";
    std::cout<<"default is 120, enter -1 to keep this value\n";
    std::cin>>this->sim_duration_d;
    if (sim_duration_d==-1)
        sim_duration_d=120;

    std::cout<<"enter the duration of the time step of the simulation in hours\n";
    std::cout<<"default is 0.01, enter -1 to keep this value\n";
    std::cin>>this->time_step_d;
    if (time_step_d==-1)
        time_step_d=0.01;

    ///Proliferation rates
    std::cout<<"enter the following proliferation rates\n";
    std::cout<<"of APC cells default value =1/240 h, enter -1 to keep this value\n";
    std::cin>>APC_free_proliferation_rate_;
    if (APC_free_proliferation_rate_==-1)
        APC_free_proliferation_rate_=1.0/240;

    std::cout<<"of NK cells default value =1/240 h, enter -1 to keep this value\n";
    std::cin>>NK_max_proliferation_rate_;
    if (NK_max_proliferation_rate_==-1)
        NK_max_proliferation_rate_=1.0/240;

    std::cout<<"of LT cells that do not express the receptor, default value  1/96 1/h, enter -1 to keep this value \n";
    std::cin>>LT_max_no_receptor_prol_rate_;
    if (LT_max_no_receptor_prol_rate_==-1)
        LT_max_no_receptor_prol_rate_=1.0/96;

    std::cout<<"of LT cells that do express the receptor and it it is free, default value 1/2, enter -1 to keep this value \n";
    std::cin>>LT_max_free_prol_rate_;
    if (LT_max_free_prol_rate_==-1)
        LT_max_free_prol_rate_=1.0/2;

    std::cout<<"of LT cells where the receptor have been bound, default value 1/1.2, enter -1 to keep this value \n";
    std::cin>>LT_max_bound_prol_rate_;
    if (LT_max_bound_prol_rate_==-1)
        LT_max_bound_prol_rate_=1.0/1.2;

    ///Conversion rates
    ///APC
    std::cout<<"\n\n enter the following conversion rates in APC cells\n";
    std::cout<<"of Antigen internalization default it needs 30 h for 1 ug/ml of antigen, enter -1 to keep this value \n";
    std::cin>>APC_no_to_free_rate_per_Ag_;
    if (APC_no_to_free_rate_per_Ag_==-1)
        APC_no_to_free_rate_per_Ag_=1.0/30;
    std::cout<<"of lingand receptor binding default value it needs 1 h for a population of 100e3, enter -1 to keep this value \n";
    std::cin>>APC_free_to_bound_rate_per_LT_;
    if (APC_free_to_bound_rate_per_LT_==-1)
        APC_free_to_bound_rate_per_LT_=1.0/1e5;
    std::cout<<"of Ab binding";
    std::cout<<"defaut is 1/1e5, enter -1 to keep this value\n";
    std::cin>>APC_Ab_binding_rate_;
    if (APC_Ab_binding_rate_==-1)
        APC_Ab_binding_rate_=1/1e5;
    std::cout<<"of exhausted APC";
    std::cout<<"defaut is 1/1e5, enter -1 to keep this value\n";
    std::cin>>APC_exh_rate;
    if (APC_exh_rate==-1)
        APC_exh_rate=1/1e5;

    ///NK
    std::cout<<"\n\n enter the following conversion rates in NK cells\n";
    std::cout<<"of Antigen internalization default it needs 30 h for 1 ug/ml of antigen, enter -1 to keep this value \n";
    std::cin>>NK_no_to_free_rate_per_Ag_;
    if (NK_no_to_free_rate_per_Ag_==-1)
        NK_no_to_free_rate_per_Ag_=1.0/30;
    std::cout<<"of lingand receptor binding default value it needs 1 h for a population of 100e3, enter -1 to keep this value \n";
    std::cin>>NK_free_to_bound_rate_per_LT_;
    if (NK_free_to_bound_rate_per_LT_==-1)
        NK_free_to_bound_rate_per_LT_=1.0/1e5;
    std::cout<<"of Ab binding";
    std::cout<<"defaut is 1/1e5, enter -1 to keep this value\n";
    std::cin>>NK_Ab_binding_rate;
    if (NK_Ab_binding_rate==-1)
        NK_Ab_binding_rate=1/1e5;
    std::cout<<"of exhausted NK";
    std::cout<<"defaut is 1/1e5, enter -1 to keep this value\n";
    std::cin>>NK_exh_rate;
    if (NK_exh_rate==-1)
        NK_exh_rate=1/1e5;

    ///LT
    std::cout<<"\n\n enter the following conversion rates in LT cells\n";
    std::cout<<"of receptor expression per APC cell default value is 6 hours for 1e5 cells, enter -1 to keep this value\n";
    std::cin>>LT_no_to_free_rate_per_APC_;
    if (LT_no_to_free_rate_per_APC_==-1)
        LT_no_to_free_rate_per_APC_=1.0/6e5;
    std::cout<<"of lingand receptor binding per APC cell  default value is 1h for a population of 1e5 cells, enter -1 to keep this value\n";
    std::cin>>LT_free_to_bound_rate_per_APC_;
    if (LT_free_to_bound_rate_per_APC_==-1)
        LT_free_to_bound_rate_per_APC_=1.0/1e5;

    ///IFNg production rates
    ///APC
    std::cout<<"\n\n enter the following IFN production rates in APC cells\n";
    std::cout<<"of cells without receptor or free receptor default value is 0.5 pg per hour per 1e5 cells, enter -1 to keep this value \n";
    std::cin>>APC_IFN_free_prod_rate_;
    if (APC_IFN_free_prod_rate_==-1)
        APC_IFN_free_prod_rate_=0.5/1e5;
    std::cout<<"of cells with the internalized antigen  default value is 5 pg per hour per 1e5 cells, enter -1 to keep this value  \n";
    std::cin>>APC_IFN_Ag_prod_rate_;
    if (APC_IFN_Ag_prod_rate_==-1)
        APC_IFN_Ag_prod_rate_=5.0/1e5;
    std::cout<<"of cells with bound receptor  default value is 10 pg per hour per 1e5 cells, enter -1 to keep this value  \n";
    std::cin>>APC_IFN_bound_prod_rate_;
    if (APC_IFN_bound_prod_rate_==-1)
        APC_IFN_bound_prod_rate_=10.0/1e5;

    ///NK
    std::cout<<"\n\n enter the following IFN production rates in APC cells\n";
    std::cout<<"of cells without receptor or free receptor default value is 0.5 pg per hour per 1e5 cells, enter -1 to keep this value \n";
    std::cin>>NK_IFN_free_prod_rate_;
    if (NK_IFN_free_prod_rate_==-1)
        NK_IFN_free_prod_rate_=0.5/1e5;
    std::cout<<"of cells with the internalized antigen  default value is 5 pg per hour per 1e5 cells, enter -1 to keep this value  \n";
    std::cin>>NK_IFN_Ag_prod_rate_;
    if (NK_IFN_Ag_prod_rate_==-1)
        NK_IFN_Ag_prod_rate_=5.0/1e5;
    std::cout<<"of cells with bound receptor  default value is 10 pg per hour per 1e5 cells, enter -1 to keep this value  \n";
    std::cin>>NK_IFN_bound_prod_rate_;
    if (NK_IFN_bound_prod_rate_==-1)
        NK_IFN_bound_prod_rate_=10.0/1e5;

    ///LT
    std::cout<<"\n\n enter the following IFN production rates in LT cells\n";
    std::cout<<"of cells without receptor default value is 0.001 pg per hour per 1e5 cells, enter -1 to keep this value \n";
    std::cin>>LT_IFN_no_rec_prod_rate_;
    if (LT_IFN_no_rec_prod_rate_==-1)
        LT_IFN_no_rec_prod_rate_=0.001/1e5;
    std::cout<<"of cells with free receptor  default value is 500 pg per hour per 1e5 cells, enter -1 to keep this value\n";
    std::cin>>LT_IFN_free_prod_rate_;
    if (LT_IFN_free_prod_rate_==-1)
        LT_IFN_free_prod_rate_=101.0/1e5;
    std::cout<<"of cells with bound receptor,  default value is 1000 pg per hour per 1e5 cells, enter -1 to keep this value \n";
    std::cin>>LT_IFN_bound_prod_rate_;
    if (LT_IFN_bound_prod_rate_==-1)
        LT_IFN_bound_prod_rate_=200.0/1e5;

    ///TNF production rates
    ///APC
    std::cout<<"\n\n enter the following TNF production rates in APC cells\n";
    std::cout<<"of cells without receptor or free receptor default value is 0.5 pg per hour per 1e5 cells, enter -1 to keep this value \n";
    std::cin>>APC_TNF_free_prod_rate_;
    if (APC_TNF_free_prod_rate_==-1)
        APC_TNF_free_prod_rate_=5/1e5;
    std::cout<<"of cells with the internalized antigen  default value is 5 pg per hour per 1e5 cells, enter -1 to keep this value  \n";
    std::cin>>APC_TNF_Ag_prod_rate_;
    if (APC_TNF_Ag_prod_rate_==-1)
        APC_TNF_Ag_prod_rate_=570/1e5;
    std::cout<<"of cells with bound receptor  default value is 10 pg per hour per 1e5 cells, enter -1 to keep this value  \n";
    std::cin>>APC_TNF_bound_prod_rate_;
    if (APC_TNF_bound_prod_rate_==-1)
        APC_TNF_bound_prod_rate_=1110/1e5;

    ///NK
    std::cout<<"\n\n enter the following TNF production rates in APC cells\n";
    std::cout<<"of cells without receptor or free receptor default value is 0.5 pg per hour per 1e5 cells, enter -1 to keep this value \n";
    std::cin>>NK_TNF_free_prod_rate_;
    if (NK_TNF_free_prod_rate_==-1)
        NK_TNF_free_prod_rate_=5/1e5;
    std::cout<<"of cells with the internalized antigen  default value is 5 pg per hour per 1e5 cells, enter -1 to keep this value  \n";
    std::cin>>NK_TNF_Ag_prod_rate_;
    if (NK_TNF_Ag_prod_rate_==-1)
        NK_TNF_Ag_prod_rate_=570/1e5;
    std::cout<<"of cells with bound receptor  default value is 10 pg per hour per 1e5 cells, enter -1 to keep this value  \n";
    std::cin>>NK_TNF_bound_prod_rate_;
    if (NK_TNF_bound_prod_rate_==-1)
        NK_TNF_bound_prod_rate_=1110/1e5;
    std::cout<<"of cells with bound receptor  default value is 10 pg per hour per 1e5 cells, enter -1 to keep this value  \n";
    std::cin>>NK_TNF_blocked_prod_rate_;
    if (NK_TNF_blocked_prod_rate_==-1)
        NK_TNF_blocked_prod_rate_=1110/1e5;

    ///LT
    std::cout<<"\n\n enter the following TNF production rates in LT cells, enter -1 to keep this value\n";
    std::cout<<"of cells without receptor default value is 0.001 pg per hour per 1e5 cells, enter -1 to keep this value \n";
    std::cin>>LT_TNF_no_rec_prod_rate_;
    if (LT_TNF_no_rec_prod_rate_==-1)
        LT_TNF_no_rec_prod_rate_=0.001/1e5;
    std::cout<<"of cells with free receptor  default value is 500 pg per hour per 1e5 cells, enter -1 to keep this value\n";
    std::cin>>LT_TNF_free_prod_rate_;
    if (LT_TNF_free_prod_rate_==-1)
        LT_TNF_free_prod_rate_=10.0/1e5;
    std::cout<<"of cells with bound receptor,  default value is 1000 pg per hour per 1e5 cells, enter -1 to keep this value \n";
    std::cin>>LT_IFN_bound_prod_rate_;
    if (LT_TNF_bound_prod_rate_==-1)
        LT_TNF_bound_prod_rate_=20.0/1e5;

    IFN_deg=0;
    LT_TNF_bound_prod_rate_=0;
    APC_IFN_blocked_prod_rate_=0;
    APC_TNF_blocked_prod_rate_=0;
    NK_IFN_blocked_prod_rate_=0;
    LT_max_blocked_prol_rate_=0;
    LT_IFN_blocked_prod_rate_=0;
    LT_TNF_blocked_prod_rate_=0;
    LT_mAb_binding_rate_=0;
    LT_TNF_bound_prod_rate_=0;



    std::cout<<"Thanks";



    m=Media(max_num_cells_,
            init_num_APC_cells+init_num_LT_cells+init_num_NK_cells,
            0,
            0,
            0,
            0,
            TNF_deg,
            IFN_deg
            //  Ag_internalization_rate
            );

    APC=APC_cells(init_num_APC_cells,
                  APC_free_proliferation_rate_,
                  APC_Ag_proliferation_rate_,
                  APC_bound_proliferation_rate_,
                  APC_blocked_proliferation_rate_,
                  APC_exh_proliferation_rate_,
                  APC_free_apop_rate_,
                  APC_Ag_apop_rate_,
                  APC_bound_apop_rate_,
                  APC_blocked_apop_rate_,
                  APC_exh_apop_rate_,
                  APC_no_to_free_rate_per_Ag_ ,
                  APC_free_to_bound_rate_per_LT_,
                  APC_Ab_binding_rate_,
                  APC_exh_rate,
                  APC_IFN_free_prod_rate_,
                  APC_IFN_Ag_prod_rate_,
                  APC_IFN_bound_prod_rate_,
                  APC_IFN_blocked_prod_rate_,
                  APC_IFN_exh_prod_rate_,
                  percentage_APC_IFN_free_prod_rate_,
                  percentage_APC_IFN_Ag_prod_rate_,
                  percentage_APC_IFN_bound_prod_rate_,
                  percentage_APC_IFN_blocked_prod_rate_,
                  percentage_APC_IFN_exh_prod_rate_,

                  APC_TNF_free_prod_rate_,
                  APC_TNF_Ag_prod_rate_,
                  APC_TNF_bound_prod_rate_,
                  APC_TNF_blocked_prod_rate_,
                  APC_TNF_exh_prod_rate_,
                  percentage_APC_TNF_free_prod_rate_,
                  percentage_APC_TNF_Ag_prod_rate_,
                  percentage_APC_TNF_bound_prod_rate_,
                  percentage_APC_TNF_blocked_prod_rate_,
                  percentage_APC_TNF_exh_prod_rate_,
                  );

    NK=NK_cells(init_num_NK_cells,
                NK_max_proliferation_rate_,
                NK_no_to_free_rate_per_Ag_ ,
                NK_free_to_bound_rate_per_LT_,
                NK_Ab_binding_rate,
                NK_exh_rate,
                NK_IFN_free_prod_rate_,
                NK_IFN_Ag_prod_rate_,
                NK_IFN_bound_prod_rate_,
                NK_IFN_blocked_prod_rate_,
                NK_TNF_free_prod_rate_,
                NK_TNF_Ag_prod_rate_,
                NK_TNF_bound_prod_rate_,
                NK_TNF_blocked_prod_rate_
                );


    LT=LT_cells(init_num_LT_cells,
                LT_num_specific,
                LT_max_no_receptor_prol_rate_,
                LT_max_free_prol_rate_,
                LT_max_bound_prol_rate_,
                LT_max_blocked_prol_rate_,
                LT_IFN_no_rec_prod_rate_,
                LT_IFN_free_prod_rate_,
                LT_IFN_bound_prod_rate_,
                LT_IFN_blocked_prod_rate_,
                LT_TNF_no_rec_prod_rate_,
                LT_TNF_free_prod_rate_,
                LT_TNF_bound_prod_rate_,
                LT_TNF_blocked_prod_rate_,
                LT_no_to_free_rate_per_APC_,
                LT_free_to_bound_rate_per_APC_,
                LT_mAb_binding_rate_,
                LT_exh_rate_
                );

}*/
void Cell_simulator::run()
{
    filename="out.txt";
    std::ofstream f;
    f.open(filename.c_str(),std::ios_base::app);
    trun_d=0;
/*1*/   std::cout<<"round"<<"\t";
/*2*/    std::cout<<"IFN"<<"\t";
/*3*/    std::cout<<"TNF"<<"\t";
/*4*/    std::cout<<"APC"<<"\t";
/*5*/    std::cout<<"NK"<<"\t";
/*6*/    std::cout<<"LT()"<<"\n";
/*7*/    std::cout<<std::endl;



/*1*/    f<<"round"<<" , ";

/*2*/    f<<"IFNamma[]"<<" , ";
/*3*/    f<<"TNF[]"<<" , ";

/*4*/    f<<"Total APC"<<" , ";
/*5*/    f<<"Total NK"<<" , ";
/*6*/    f<<"Total LT"<<" , ";

/*7*/    f<<"APC0"<<" , ";
/*8*/    f<<"APCa "<<" , ";
/*9*/    f<<"APCbo"<<" , ";
/*10*/    f<<"APCbo_bl"<<" , ";
/*11*/    f<<"APC_bl"<<" , ";
/*12*/    f<<"APC exh"<<" , ";

/*13*/    f<<"%APC expresing receptor"<<" , ";
/*14*/    f<<"APC.IFNgamma_production_rate"<<" , ";
/*15*/    f<<"APC.percentage of IFN producing cells"<<" , ";
/*16*/    f<<"APC.TNF_production_rate"<<" , ";
/*17*/    f<<"APC.percentage of TNF producing cells";

/*18*/    f<<"NK0"<<" , ";
/*19*/    f<<"NKa"<<" , ";
/*20*/    f<<"NKbo"<<" , ";
/*21*/    f<<"NKbo_bl"<<" , ";
/*22*/    f<<"NKbl"<<" , ";
/*23*/    f<<"NK exh"<<" , ";

/*24*/    f<<"%NK expresing receptor"<<" , ";
/*25*/    f<<"NK.IFNgamma_production_rate"<<" , ";
/*26*/    f<<"NK.percentage of IFN producing cell"<<" , ";
/*27*/    f<<"NK.TNF_production_rate"<<" , ";
/*28*/    f<<"NK.percentage of TNF producing cell"<<" , ";

/*29*/    f<<"LT no Agsp"<<" , ";
/*30*/    f<<"LT0"<<" , ";
/*31*/    f<<"LTbo"<<" , ";
/*32*/    f<<"LTbl"<<" , ";
/*33*/    f<<"LTexh"<<" , ";

/*34*/    f<<"%LT expresing receptor"<<" , ";
/*35*/    f<<"LT.IFNgamma_production_rate"<<" , ";
/*36*/    f<<"LT.percentage of IFN producing cell"<<" , ";
/*37*/    f<<"LT.TNF_production_rate"<<" , ";
/*38*/    f<<"LT.percentage of TNF producing cell"<<" , ";
/*39*/    f<<"LT undergoing apoptosis"<<" , ";

/*40*/    f<<"Tymidine incorporated"<<" , ";
/*41*/    f<<"Ag"<<" , ";
/*42*/    f<<"Ab"<<" , ";
          f<<"\n";

    double eps=1e-7;
    while (trun_d<this->sim_duration_d)
    {

        if (trun_d+eps-floor(trun_d+eps)<time_step_d)
        {
            std::cout<<trun_d<<"\t";
            std::cout<<m.IFNgamma()<<"\t";
            std::cout<<m.TNF()<<"\t";
            std::cout<<APC.num_APC()<<"\t";
            std::cout<<NK.num_NK()<<"\t";
            std::cout<<LT.num_LT()<<"\n";


/*1*/            f<<trun_d<<" , ";
/*2*/            f<<m.IFNgamma()<<" , ";
/*3*/            f<<m.TNF()<<" , ";

/*4*/            f<<APC.num_APC()<<" , ";
/*5*/            f<<NK.num_NK()<<" , ";
/*6*/            f<<LT.num_LT()<<" , ";

/*7*/            f<<APC.APC0()<<" , ";
/*8*/            f<<APC.APCa()<<" , ";
/*9*/            f<<APC.APCbo()<<" , ";
/*10*/           f<<APC.APCbo_Ab()<<" , ";
/*11*/           f<<APC.APCbl()<<" , ";
/*12*/           f<<APC.APCexh()<<" , ";


/*13*/           f<<APC.percentage_cell_expressing_receptor()<<" , ";
/*14*/           f<<APC.APC_IFNgamma_production_rate()<<" , ";
/*15*/           f<<APC.percentage_APC_producing_IFN()<<" , ";
/*16*/           f<<APC.APC_TNF_production_rate()<<" , ";
/*17*/           f<<APC.percentage_APC_producing_TNF()<<" , ";


/*18*/           f<<NK.NK0()<<" , ";
/*19*/           f<<NK.NKa()<<" , ";
/*20*/           f<<NK.NKbo()<<" , ";
/*21*/           f<<NK.NKbo_Ab()<<" , ";
/*22*/           f<<NK.NKbl()<<" , ";
/*23*/           f<<NK.NKexh()<<" , ";

/*24*/           f<<NK.percentage_NK_expressing_receptor()<<" , ";
/*25*/           f<<NK.NK_IFNgamma_production_rate()<<" , ";
/*26*/           f<<NK.percentage_NK_producing_IFN()<<" , ";
/*27*/           f<<NK.NK_TNF_production_rate()<<" , ";
/*28*/           f<<NK.percentage_NK_producing_TNF()<<" , ";

/*29*/    f<<LT.LTns()<<" , ";
/*30*/    f<<LT.LT0()<<" , ";
/*31*/    f<<LT.LTbo()<<" , ";
/*32*/    f<<LT.LTbl()<<" , ";
/*33*/    f<<LT.LTexh()<<" , ";

/*34*/    f<<LT.LT_percentage_cell_expressing_receptor()<<" , ";
/*35*/    f<<LT.LT_IFNgamma_production_rate()<<" , ";
/*36*/    f<<LT.percentage_LT_IFN_production()<<" , ";
/*37*/    f<<LT.TNF_production_rate()<<" , ";
/*38*/    f<<LT.percentage_LT_TNF_production()<<" , ";

/*39*/    f<<LT.percentage_apoptotic_LT_cells()<<" , ";

/*40*/    f<<m.Tymidine_incorporated()<<" , ";
/*41*/    f<<m.Ag()<<" , ";
/*42*/    f<<m.Ab()<<"\n";

        };

        APC.update(time_step_d,m,NK,LT);
        NK.update(time_step_d,m,APC,LT);
        LT.update(time_step_d,trun_d,m,APC,NK);
        m.update(time_step_d,trun_d,APC,NK,LT);
        trun_d+=time_step_d;

    }
    f.close();
}


Cell_simulator::Cell_simulator(const SimParameters& sp,
                               const Experiment& E):

    m(),
    APC(),
    NK(),
    LT(),

    time_step_d(),
    sim_duration_d(),
    trun_d(),
    filename(),


    experiment_(E),
    fitExperiment_(),
    initialPar_(sp),
    fitPar_()
{}


Cell_simulator& Cell_simulator::applyParameters(const SimParameters& sp,
						const Treatment& tr)
{
    m=Media(0,
            0,0,sp.TNF_deg_,sp.IFN_deg_,tr.init_cells,tr.Ag,tr.Ab,sp.TymidineTriteate_,sp.Prol_TymTr_);


    APC=APC_cells(sp.init_ratio_APC_*tr.init_cells,
                  /*2*/ sp.IFN_APC0_prod_rate_,
                  /*3*/ sp.IFN_APCa_prod_rate_,
                  /*4*/ sp.IFN_APCbo_prod_rate_,


                  /// 3) TNF Poductions rates of each type of APC
                  /*5*/ sp.TNF_APC0_prod_rate_,
                  /*6*/ sp.TNF_APCa_prod_rate_,
                  /*7*/ sp.TNF_APCbo_prod_rate_,


                  /// 4) Percentages of IFN productions of each type of APC
                  /*8*/ sp.percentage_IFN_APC0_prod_rate_,
                  /*9*/ sp.percentage_IFN_APCa_prod_rate_,
                  /*10*/ sp.percentage_IFN_APCbo_prod_rate_,

                  /// 5)Percentages of TNF productions of each type of APC
                  /*11*/ sp.percentage_TNF_APC0_prod_rate_,
                  /*12*/ sp.percentage_TNF_APCa_prod_rate_,
                  /*13*/ sp.percentage_TNF_APCbo_prod_rate_,


                  /// 6) Proliferation rates
                  /*14*/ sp.APC_bound_proliferation_rate_,

                  /// 7) Apoptosis rates
                  /*15*/ sp.APC0_apop_rate_,
                  /*16*/ sp.APCa_apop_rate_,
                  /*17*/ sp.APCbo_apop_rate_,
                  /*18*/ sp.APCbl_apop_rate_,
                  /*19*/ sp.APCexh_apop_rate_,

                  /// 8) constant saturation of TNF for apoptosis
                  /*20*/ sp.Ks_APC_m_TNF_,

                  /// 9) conversion rates
                  /*21*/ sp.APC_Ag_,
                  /*22*/ sp.APC_APC_,
                  /*23*/ sp.APC_NK_,
                  /*24*/ sp.APC_LT_1_,
                  /*25*/ sp.APC_LT_2_,
                  /*26*/ sp.APC_Ab_,
                  /*27*/ sp.APC_exh_,

                  /// 10)Saturation constant of IFN and TNF for activation
                  /*28*/ sp.KsAPC_LT_,

                  /// 11)Saturation constant of APC_LT interaction
                  /*29*/ sp.APC_Ksi_,
                  /*30*/ sp.APC_Kst_,

                  /// 12) Percentages of cell expressing receptor
                  /*31*/ sp.APC0_expressing_receptor_,
                  /*32*/ sp.APCa_expressing_receptor_,
                  /// 13) Apoptosis rate for TNF
                  /*33*/ sp.u_APC_TNF_);

    NK=NK_cells (sp.init_ratio_NK_*tr.init_cells,
                 /// 2) IFN Poductions rates of each type of NK
                 /*2*/ sp.IFN_NK0_prod_rate_,
                 /*3*/ sp.IFN_NKa_prod_rate_,
                 /*4*/ sp.IFN_NKbo_prod_rate_,

                 /// 3) TNF Poductions rates of each type of NK
                 /*5*/ sp.TNF_NK0_prod_rate_,
                 /*6*/ sp.TNF_NKa_prod_rate_,
                 /*7*/ sp.TNF_NKbo_prod_rate_,

                 /// 4) Percentages of IFN productions of each type of NK
                 /*8*/ sp.percentage_IFN_NK0_prod_rate_,
                 /*9*/ sp.percentage_IFN_AgNKa_prod_rate_,
                 /*10*/ sp.percentage_IFN_NKbo_prod_rate_,

                 /// 5)Percentages of TNF productions of each type of NK
                 /*11*/ sp.percentage_TNF_NK0_prod_rate_,
                 /*12*/ sp.percentage_TNF_NKa_prod_rate_,
                 /*13*/ sp.percentage_TNF_NKbo_prod_rate_,

                 /// 6) Proliferation rates
                 /*13.5*/ sp.NK0_proliferation_rate_,
                 /*14*/ sp.NKa_proliferation_rate_,
                 /*15*/ sp.NKbo_proliferation_rate_,
                 /*16*/ sp.NKbl_proliferation_rate_,

                 /// 7) Apoptosis rates
                 /*17*/ sp.NK0_apop_rate_,
                 /*18*/ sp.NKa_apop_rate_,
                 /*19*/ sp.NKbo_apop_rate_,
                 /*20*/ sp.NKbl_apop_rate_,
                 /*21*/ sp.NKexh_apop_rate_,

                 /// 8) constant saturation of TNF for apoptosis
                 /*22*/ sp.Ks_NK_m_TNF_,

                 /// 9) conversion rates
                 /*23*/ sp.KaNK_,
                 /*24*/ sp.NK_NK_,
                 /*25*/ sp.NK_Ab_,
                 /*26*/ sp.NK_exh_,

                 /// 10)Saturation constant of NK interaction for activation
                 /*27*/ sp.KsAPC_NK_,

                 /// 11)Saturation constant of NK_LT interaction
                 /*28*/ sp.NK_Ksi_,
                 /*29*/ sp.NK_Kst_,

                 /// 12) Percentages of cell expressing receptor
                 /*30*/ sp.NK0_expressing_receptor_,
                 /*31*/ sp.NKa_expressing_receptor_,

                 /// 13) Apoptosis rate for TNF
                 /*32*/ sp.u_NK_TNF_);


    LT=LT_cells  (sp.ratio_init_LTns_*tr.init_cells,
                  sp.ratio_initLTspecific_*tr.init_cells,
                  /// 2) IFN Poductions rates of each type of LT
                     /*3*/ sp.IFN_LTns_prod_rate_,
                     /*4*/ sp.IFN_LTbo_prod_rate_,
                     /*5*/ sp.IFN_LTbl_prod_rate_,

                 /// 3) TNF Poductions rates of each type of LT
                     /*6*/ sp.TNF_LTns_prod_rate_,
                     /*7*/ sp.TNF_LTbo_prod_rate_,
                     /*8*/ sp.TNF_LTbl_prod_rate_,

                 /// 4) Percentages of IFN productions of each type of LT
                     /*9*/ sp.percentage_IFN_LTns_prod_rate_,
                     /*10*/ sp.percentage_IFN_LTbo_prod_rate_,
                     /*11*/ sp.percentage_IFN_LTbl_prod_rate_,


                 /// 5)Percentages of TNF productions of each type of LT
                     /*12*/ sp.percentage_TNF_LTns_prod_rate_,
                     /*13*/ sp.percentage_TNF_LTbo_prod_rate_,
                     /*14*/ sp.percentage_TNF_LTbl_prod_rate_,

                 /// 6) Proliferation rates
                     /*15*/ sp.LTns_proliferation_rate_,
                     /*16*/ sp.LTbo_proliferation_rate_,
                     /*17*/ sp.LTbl_proliferation_rate_,

                 /// 7) Apoptosis rates
                     /*18*/ sp.LTns_apop_rate_,
                     /*19*/ sp.LTbo_apop_rate_,
                     /*20*/ sp.LTbl_apop_rate_,
                     /*21*/ sp.LTexh_apop_rate_,

                 /// 8) constant saturation of TNF for apoptosis
                     /*22*/ sp.Ks_LT_m_TNF_,

                 /// 9) Percentages of cell expressing receptor
                     /*23*/ sp.LTns_expressing_receptor_,

                 /// 10) Apoptosis rate for TNF
                     /*24*/ sp.u_LT_TNF_,

                  /// 11) LT exh rate
                      /*25*/ sp.LT_exh_rate_,

                  /// 12) apoptosis related parameters
                      /*26*/ sp.t_apop_meas_,
                      /*27*/ sp.t_duration_apoptosis_
        );

    time_step_d=tr.time_step_d;
    sim_duration_d=tr.sim_duration_d;
    trun_d=0;

return *this;
}



Cell_simulator::Cell_simulator(){}




Cell_simulator::Cell_simulator(const Cell_simulator& other):
    m(other.m),
    APC(other.APC),
    NK(other.NK),
    LT(other.LT),

    time_step_d(other.time_step_d),
    sim_duration_d(other.sim_duration_d),
    trun_d(other.trun_d),
    filename (other.filename),


    experiment_(other.experiment_),
    fitExperiment_(other.fitExperiment_),
    initialPar_(other.initialPar_),
    fitPar_(fitPar_)
    {}


Cell_simulator&
Cell_simulator::operator=(const Cell_simulator& other)
{
    if (this!=&other)
    {
        Cell_simulator tmp(other);
        swap(*this,tmp);
    }
    return *this;
}

void swap(Cell_simulator& one, Cell_simulator& other)
{
    std::swap(one.m,other.m);
    std::swap(one.APC,other.APC);
    std::swap(one.NK,other.NK);
    std::swap(one.LT,other.LT);
    std::swap(one.time_step_d,other.time_step_d);
    std::swap(one.sim_duration_d,other.sim_duration_d);
    std::swap(one.trun_d,other.trun_d);
    std::swap(one.filename,other.filename);
    std::swap(one.experiment_,other.experiment_);
    std::swap(one.fitExperiment_,other.fitExperiment_);
    std::swap(one.initialPar_,other.initialPar_);
    std::swap(one.fitPar_,other.fitPar_);
    }

/**
  @param results the method Simulate only uses  the Time() of
  each measurement.
  */


Results Cell_simulator::Simulate(const SimParameters& simPar,
                                 const Treatment& tr,
                                 const Results& results)
{
    *this=applyParameters(simPar, tr);

    double Duratione=results.Duration();


    std::vector<Measurement> TNFs(results.TNF());
    std::size_t iTNFs=0;

    double tTNFs;
    if (!TNFs.empty())
        tTNFs=TNFs[iTNFs].Time();
    else
        tTNFs=Duratione+1;

    std::vector<Measurement> IFNs=results.IFN();
    std::size_t iIFNs=0;
    double tIFNs;
    if (!IFNs.empty())
        tIFNs=IFNs[iIFNs].Time();
    else
        tIFNs=Duratione+1;

    std::vector<Measurement> APC_exp=results.APC_expression();
    std::size_t iAPC_exp=0;
    double tAPC_exp;
    if (!APC_exp.empty())
        tAPC_exp=APC_exp[iAPC_exp].Time();
    else
        tAPC_exp=Duratione+1;

    std::vector<Measurement> NK_exp=results.NK_expression();
    std::size_t iNK_exp=0;
    double tNK_exp;
    if (!NK_exp.empty())
        tNK_exp=NK_exp[iNK_exp].Time();
    else
        tNK_exp=Duratione+1;

    std::vector<Measurement> LT_exp=results.LT_expression();
    std::size_t iLT_exp=0;
    double tLT_exp;
    if (!LT_exp.empty())
        tLT_exp=LT_exp[iLT_exp].Time();
    else tLT_exp=Duratione+1;

    std::vector<Measurement> APC_IFNs=results.APC_IFNg();
    std::size_t iAPC_IFN=0;
    double tAPC_IFN;
    if (!APC_IFNs.empty())
        tAPC_IFN=APC_IFNs[iAPC_IFN].Time();
    else tAPC_IFN=Duratione+1;

    std::vector<Measurement> APC_TNFs=results.APC_TNFa();
    std::size_t iAPC_TNF=0;
    double tAPC_TNF;
    if (!APC_TNFs.empty())
        tAPC_TNF=APC_TNFs[iAPC_TNF].Time();
    else tAPC_TNF=Duratione+1;

    std::vector<Measurement> NK_IFNs=results.NK_IFNg();
    std::size_t iNK_IFN=0;
    double tNK_IFN;
    if (!NK_IFNs.empty())
        tNK_IFN=NK_IFNs[iNK_IFN].Time();
    else tNK_IFN=Duratione+1;

    std::vector<Measurement> NK_TNFs=results.NK_TNFa();
    std::size_t iNK_TNF=0;
    double tNK_TNF;
    if (!NK_TNFs.empty())
        tNK_TNF=NK_TNFs[iNK_TNF].Time();
    else tNK_TNF=Duratione+1;

    std::vector<Measurement> LT_IFNs=results.LT_IFNg();
    std::size_t iLT_IFN=0;
    double tLT_IFN;
    if (!LT_IFNs.empty())
        tLT_IFN=LT_IFNs[iLT_IFN].Time();
    else tLT_IFN=Duratione+1;

    std::vector<Measurement> LT_TNFs=results.LT_TNFa();
    std::size_t iLT_TNF=0;
    double tLT_TNF;
    if (!LT_TNFs.empty())
        tLT_TNF=LT_TNFs[iLT_TNF].Time();
    else tLT_TNF=Duratione+1;

    std::vector<Measurement> LT_apops=results.LT_Apoptosis();
    std::size_t iLT_Apoptosis=0;
    double tLT_Apoptosis;
    if (!LT_apops.empty())
        tLT_Apoptosis=LT_apops[iLT_Apoptosis].Time();
    else tLT_Apoptosis=Duratione+1;

    std::vector<Measurement> Prols=results.Proliferation();
    std::size_t i_Proliferation=0;
    double t_Proliferation;
    if (!Prols.empty())
        t_Proliferation=Prols[i_Proliferation].Time();
    else t_Proliferation=Duratione+1;

    std::vector<Measurement> num_cellss=results.num_cells();
    std::size_t inum_cells=0;
    double t_num_cells;
    if (!num_cellss.empty())
        t_num_cells=num_cellss[inum_cells].Time();
    else t_num_cells=Duratione+1;

    double eps=1e-7;

    while (trun_d+eps<=results.Duration())
    {

        if(trun_d+eps>=tTNFs)
        {
            TNFs[iTNFs].setMeasurement(m.TNF());
            ++iTNFs;
            if (iTNFs<TNFs.size())
            {
                tTNFs=TNFs[iTNFs].Time();
            }
            else
            {
                tTNFs=results.Duration()+1;
            }
        };

        if(trun_d+eps>=tIFNs)
        {
            IFNs[iIFNs].setMeasurement(m.IFNgamma());
            ++iIFNs;
            if (iIFNs<IFNs.size())
            {
                tIFNs=IFNs[iIFNs].Time();
            }
            else
            {
                tIFNs=results.Duration()+1;
            }
        };


      if(trun_d+eps>=tAPC_exp)
        {
            APC_exp[iAPC_exp].setMeasurement(APC.percentage_cell_expressing_receptor());
            ++iAPC_exp;
            if (iAPC_exp<APC_exp.size())
            {
                tAPC_exp=APC_exp[iAPC_exp].Time();
            }
            else
            {
                tAPC_exp=results.Duration()+1;
            }
        };



        if(trun_d+eps>=tNK_exp)
        {
            NK_exp[iNK_exp].setMeasurement (
                        NK.percentage_NK_expressing_receptor());
            ++iNK_exp;
            if (iNK_exp<NK_exp.size())
            {
                tNK_exp=NK_exp[iNK_exp].Time();
            }
            else
            {
                tNK_exp=results.Duration()+1;
            }
        };


        if(trun_d+eps>=tLT_exp)
        {
            LT_exp[iLT_exp].setMeasurement (
                        LT.LT_percentage_cell_expressing_receptor());
            ++iLT_exp;
            if (iLT_exp<LT_exp.size())
            {
                tLT_exp=LT_exp[iLT_exp].Time();
            }
            else
            {
                tLT_exp=results.Duration()+1;
            }



        };

        if(trun_d+eps>=tAPC_IFN)
        {
            APC_IFNs[iAPC_IFN].setMeasurement (
                        APC.percentage_APC_producing_IFN());
            ++iAPC_IFN;
            if (iAPC_IFN<APC_IFNs.size())
            {
                tAPC_IFN=APC_IFNs[iAPC_IFN].Time();
            }
            else
            {
                tAPC_IFN=results.Duration()+1;
            }

          };

        if(trun_d+eps>=tAPC_TNF)
        {
            APC_TNFs[iAPC_TNF].setMeasurement (
                         APC.percentage_APC_producing_TNF());
            ++iAPC_TNF;
            if (iAPC_TNF<APC_TNFs.size())
            {
                tAPC_TNF=APC_TNFs[iAPC_TNF].Time();
            }
            else
            {
                tAPC_TNF=results.Duration()+1;
            }

             };

        if(trun_d+eps>=tNK_IFN)
        {
          NK_IFNs[iNK_IFN].setMeasurement (
          NK.percentage_NK_producing_IFN());
          ++iNK_IFN;
               if (iNK_IFN<NK_IFNs.size())
                  {
                    tNK_IFN=NK_IFNs[iNK_IFN].Time();
                  }
               else
                  {
                   tNK_IFN=results.Duration()+1;
                  }

                };

        if(trun_d+eps>=tNK_TNF)
        {
         NK_TNFs[iNK_TNF].setMeasurement (
         NK.percentage_NK_producing_TNF());
         ++iNK_TNF;
          if (iNK_TNF<NK_TNFs.size())
          {
           tNK_TNF=NK_TNFs[iNK_TNF].Time();
          }
          else
          {
          tNK_TNF=results.Duration()+1;
          }

         };

        if(trun_d+eps>=tLT_IFN)
        {
         LT_IFNs[iLT_IFN].setMeasurement (
         LT.percentage_LT_IFN_production());
         ++iLT_IFN;
         if (iLT_IFN<LT_IFNs.size())
            {
              tLT_IFN=LT_IFNs[iLT_IFN].Time();
            }
         else
            {
              tLT_IFN=results.Duration()+1;
            }

                    };

         if(trun_d+eps>=tLT_TNF)
         {
           LT_TNFs[iLT_TNF].setMeasurement (
           LT.percentage_LT_TNF_production());
           ++iLT_TNF;
           if (iLT_TNF<LT_TNFs.size())
              {
               tLT_TNF=LT_TNFs[iLT_TNF].Time();
              }
          else
              {
               tLT_TNF=results.Duration()+1;
              }

                    };

         if(trun_d+eps>=tLT_Apoptosis)
         {
          LT_apops[iLT_Apoptosis].setMeasurement (
          LT.percentage_apoptotic_LT_cells());
          ++iLT_Apoptosis;
          if (iLT_Apoptosis<LT_apops.size())
             {
               tLT_Apoptosis=LT_apops[iLT_Apoptosis].Time();
             }
          else
             {
               tLT_Apoptosis=results.Duration()+1;
             }

                     };

          if(trun_d+eps>=t_Proliferation)
          {
            Prols[i_Proliferation].setMeasurement (
            m.Tymidine_incorporated());
            ++i_Proliferation;
            if (i_Proliferation<Prols.size())
               {
                t_Proliferation=Prols[i_Proliferation].Time();
               }
           else
               {
                t_Proliferation=results.Duration()+1;
               }

                     };

          if(trun_d+eps>=t_num_cells)
          {
              num_cellss[inum_cells].setMeasurement(m.TNF());
              ++inum_cells;
              if (inum_cells<num_cellss.size())
              {
                  t_num_cells=num_cellss[inum_cells].Time();
              }
              else
              {
                  t_num_cells=results.Duration()+1;
              }
          };


      //  char ch;

      //  std::cout<<"\n media\n"<<m<<"\nLT\n"<<APC<<"\nNK\n"<<NK<<"\nLT\n"<<LT;
      //  std::cout<<"\ntime \t"<<trun_d;
      //  std::cin>>ch;
	APC.update(time_step_d,m,NK,LT);
        NK.update(time_step_d,m,APC,LT);
        LT.update(time_step_d,trun_d,m,APC,NK);
        m.update(time_step_d,trun_d, APC,NK,LT);




        trun_d+=time_step_d;

    }

    Results SimRes(TNFs,IFNs,APC_exp,NK_exp,LT_exp,APC_IFNs,APC_TNFs,NK_IFNs,NK_TNFs, LT_IFNs, LT_TNFs, LT_apops,Prols, num_cellss, Duratione);
    return SimRes;
}

Experiment Cell_simulator::Simulate(const SimParameters& simPar,
                                    const Experiment& exp)
{
    Experiment sim;
    for (std::size_t i=0; i<exp.size(); i++)
    {
	Results r=Simulate(simPar,exp.Treatment_i(i),exp.Result_i(i));
        sim.push_back(exp.Treatment_i(i),
		      r);
    }

    return sim;

}


OptimizationResults Cell_simulator::Optimize(const SimParameters& priorPar,
                                             const SimParameters& initPar,
                                             const Experiment& experiment,
					     double range,
					     std::size_t numEval)
{
    experiment_=experiment;
    SimParameters runPar(initPar);
    LevenbergMarquardt LM0(this,getData(priorPar.getParameters()),
                           initPar.getParameters());
    LM0.optimize();
   runPar.applyParameters(LM0.OptimParameters());

   std::cout<<Simulate(runPar,experiment);
   char ch;
   std::cin>>ch;
   std::size_t iEvals=LM0.numIter();
    while (iEvals<numEval)
    {
	LevenbergMarquardt LM(this,
                              getData(priorPar.getParameters()),
			      runPar.getRandomParameters(range));
	LM.optimize();
	if (LM.SS()<LM0.SS())
	{
	    LM0=LM;
	    runPar.applyParameters(LM.OptimParameters());
	}
	iEvals+=LM.numIter();


    }
    fitPar_=priorPar;
    fitPar_.applyParameters(LM0.OptimParameters());

    OptimizationResults O;
    O.InitialParameters_=initPar;
    O.OptimalParameters_=fitPar_;
    O.Data_=experiment;
    O.FittedData_=Simulate(fitPar_,experiment);
    O.SS_=LM0.SS();
    return O;


}



std::vector<double> Cell_simulator::yfit (const std::vector<double>& param0)
{
    SimParameters par(initialPar_);
    par.applyParameters(param0);
    fitExperiment_=Simulate(par,experiment_);

//param [i], modificar el peso de los parametros dividiendo
    std::vector<double> f=fitExperiment_.getData();
    std::vector<double> param=param0;

  f.insert(f.end(),param.begin(),param.end());

    return f;
 }


std::vector<double> Cell_simulator::yfit (const Parameters& param0)
{
    Parameters par(prior_);
    par.applyParameters(param0);
    fitExperiment_=Simulate(par,experiment_);

//param [i], modificar el peso de los parametros dividiendo
    std::vector<double> f=fitExperiment_.getData();
     return f;
 }




std::vector<double> Cell_simulator::difParam(const std::vector<double>& param)
{
    std::vector<double> result;
    std::vector<double> p0=initialPar_.getParameters();
    for (std::size_t i=0; i<param.size();i++)
    {
	result.push_back(param[i]-p0[i]);
    }
    return result;
}
//param [i], modificar el peso de los parametros dividiendo
std::vector<double> Cell_simulator::getData(const std::vector<double>& param0)
{
    std::vector<double> f=experiment_.getData();
    for (std::size_t i=0; i<f.size();i++)
	f[i]=f[i];
    std::vector<double> param(param0);

    for (std::size_t i=0; i<param.size();i++)
	param[i]=param[i];

   f.insert(f.end(),param.begin(),param.end());
    return f;

}


Cell_simulator& Cell_simulator::applyParameters(const Parameters& sp,
                                const Treatment& tr)
{

    m=Media(0,
            0,
            0,
            sp.mean("TNF_deg"),
            sp.mean("IFN_deg"),
            tr.init_cells,
            tr.Ag,
            tr.Ab,
            sp.mean("TymidineTriteate"),
            sp.mean("Prol_TymTr"));


    APC=APC_cells(sp.mean("init_ratio_APC") * tr.init_cells,
                  /*2*/ sp.mean("IFN_APC0_prod_rate"),
                  /*3*/ sp.mean("IFN_APCa_prod_rate"),
                  /*4*/ sp.mean("IFN_APCbo_prod_rate"),


                  /// 3) TNF Poductions rates of each type of APC
                  /*5*/ sp.mean("TNF_APC0_prod_rate"),
                  /*6*/ sp.mean("TNF_APCa_prod_rate"),
                  /*7*/ sp.mean("TNF_APCbo_prod_rate"),


                  /// 4) Percentages of IFN productions of each type of APC
                  /*8*/ sp.mean("percentage_IFN_APC0_prod_rate"),
                  /*9*/ sp.mean("percentage_IFN_APCa_prod_rate"),
                  /*10*/ sp.mean("percentage_IFN_APCbo_prod_rate"),

                  /// 5)Percentages of TNF productions of each type of APC
                  /*11*/ sp.mean("percentage_TNF_APC0_prod_rate"),
                  /*12*/ sp.mean("percentage_TNF_APCa_prod_rate"),
                  /*13*/ sp.mean("percentage_TNF_APCbo_prod_rate"),


                  /// 6) Proliferation rates
                  /*14*/ sp.mean("APC_bound_proliferation_rate"),

                  /// 7) Apoptosis rates
                  /*15*/ sp.mean("APC0_apop_rate"),
                  /*16*/ sp.mean("APCa_apop_rate"),
                  /*17*/ sp.mean("APCbo_apop_rate"),
                  /*18*/ sp.mean("APCbl_apop_rate"),
                  /*19*/ sp.mean("APCexh_apop_rate"),

                  /// 8) constant saturation of TNF for apoptosis
                  /*20*/ sp.mean("Ks_APC_m_TNF"),

                  /// 9) conversion rates
                  /*21*/ sp.mean("APC_Ag"),
                  /*22*/ sp.mean("APC_APC"),
                  /*23*/ sp.mean("APC_NK"),
                  /*24*/ sp.mean("APC_LT_1"),
                  /*25*/ sp.mean("APC_LT_2"),
                  /*26*/ sp.mean("APC_Ab"),
                  /*27*/ sp.mean("APC_exh"),

                  /// 10)Saturation constant of IFN and TNF for activation
                  /*28*/ sp.mean("KsAPC_LT"),

                  /// 11)Saturation constant of APC_LT interaction
                  /*29*/ sp.mean("APC_Ksi"),
                  /*30*/ sp.mean("APC_Kst"),

                  /// 12) Percentages of cell expressing receptor
                  /*31*/ sp.mean("APC0_expressing_receptor"),
                  /*32*/ sp.mean("APCa_expressing_receptor"),
                  /// 13) Apoptosis rate for TNF
                  /*33*/ sp.mean("u_APC_TNF"));

    NK=NK_cells  (sp.mean("init_ratio_APC")*tr.init_cells,
                 /// 2) IFN Poductions rates of each type of NK
                 /*2*/ sp.mean("IFN_NK0_prod_rate"),
                 /*3*/ sp.mean("IFN_NKa_prod_rate"),
                 /*4*/ sp.mean("IFN_NKbo_prod_rate"),

                 /// 3) TNF Poductions rates of each type of NK
                 /*5*/ sp.mean("TNF_NK0_prod_rate"),
                 /*6*/ sp.mean("TNF_NKa_prod_rate"),
                 /*7*/ sp.mean("TNF_NKbo_prod_rate"),

                 /// 4) Percentages of IFN productions of each type of NK
                 /*8*/ sp.mean("percentage_IFN_NK0_prod_rate"),
                 /*9*/ sp.mean("percentage_IFN_AgNKa_prod_rate"),
                 /*10*/ sp.mean("percentage_IFN_NKbo_prod_rate"),

                 /// 5)Percentages of TNF productions of each type of NK
                 /*11*/ sp.mean("percentage_TNF_NK0_prod_rate"),
                 /*12*/ sp.mean("percentage_TNF_NKa_prod_rate"),
                 /*13*/ sp.mean("percentage_TNF_NKbo_prod_rate"),

                 /// 6) Proliferation rates
                 /*13.5*/ sp.mean("NK0_proliferation_rate"),
                 /*14*/ sp.mean("NKa_proliferation_rate"),
                 /*15*/ sp.mean("NKbo_proliferation_rate"),
                 /*16*/ sp.mean("NKbl_proliferation_rate"),

                 /// 7) Apoptosis rates
                 /*17*/ sp.mean("NK0_apop_rate"),
                 /*18*/ sp.mean("NKa_apop_rate"),
                 /*19*/ sp.mean("NKbo_apop_rate"),
                 /*20*/ sp.mean("NKbl_apop_rate"),
                 /*21*/ sp.mean("NKexh_apop_rate"),

                 /// 8) constant saturation of TNF for apoptosis
                 /*22*/ sp.mean("Ks_NK_m_TNF"),

                 /// 9) conversion rates
                 /*23*/ sp.mean("KaNK"),
                 /*24*/ sp.mean("NK_NK"),
                 /*25*/ sp.mean("NK_Ab"),
                 /*26*/ sp.mean("NK_exh"),

                 /// 10)Saturation constant of NK interaction for activation
                 /*27*/ sp.mean("KsAPC_NK"),

                 /// 11)Saturation constant of NK_LT interaction
                 /*28*/ sp.mean("NK_Ksi"),
                 /*29*/ sp.mean("NK_Kst"),

                 /// 12) Percentages of cell expressing receptor
                 /*30*/ sp.mean("NK0_expressing_receptor"),
                 /*31*/ sp.mean("NKa_expressing_receptor"),

                 /// 13) Apoptosis rate for TNF
                 /*32*/ sp.mean("u_NK_TNF_"));


    LT=LT_cells  (sp.mean("ratio_init_LTns_")*tr.init_cells,
                  sp.mean("ratio_initLTspecific_")*tr.init_cells,
                  /// 2) IFN Poductions rates of each type of LT
                     /*3*/ sp.mean("IFN_LTns_prod_rate"),
                     /*4*/ sp.mean("IFN_LTbo_prod_rate"),
                     /*5*/ sp.mean("IFN_LTbl_prod_rate"),

                 /// 3) TNF Poductions rates of each type of LT
                     /*6*/ sp.mean("TNF_LTns_prod_rate"),
                     /*7*/ sp.mean("TNF_LTbo_prod_rate"),
                     /*8*/ sp.mean("TNF_LTbl_prod_rate"),

                 /// 4) Percentages of IFN productions of each type of LT
                     /*9*/ sp.mean("percentage_IFN_LTns_prod_rate"),
                     /*10*/ sp.mean("percentage_IFN_LTbo_prod_rate"),
                     /*11*/ sp.mean("percentage_IFN_LTbl_prod_rate"),


                 /// 5)Percentages of TNF productions of each type of LT
                     /*12*/ sp.mean("percentage_TNF_LTns_prod_rate"),
                     /*13*/ sp.mean("percentage_TNF_LTbo_prod_rate"),
                     /*14*/ sp.mean("percentage_TNF_LTbl_prod_rate"),

                 /// 6) Proliferation rates
                     /*15*/ sp.mean("LTns_proliferation_rate"),
                     /*16*/ sp.mean("LTbo_proliferation_rate"),
                     /*17*/ sp.mean("LTbl_proliferation_rate"),

                 /// 7) Apoptosis rates
                     /*18*/ sp.mean("LTns_apop_rate"),
                     /*19*/ sp.mean("LTbo_apop_rate"),
                     /*20*/ sp.mean("LTbl_apop_rate"),
                     /*21*/ sp.mean("LTexh_apop_rate"),

                 /// 8) constant saturation of TNF for apoptosis
                     /*22*/ sp.mean("Ks_LT_m_TNF"),

                 /// 9) Percentages of cell expressing receptor
                     /*23*/ sp.mean("LTns_expressing_receptor"),

                 /// 10) Apoptosis rate for TNF
                     /*24*/ sp.mean("u_LT_TNF"),

                  /// 11) LT exh rate
                      /*25*/ sp.mean("LT_exh_rate"),

                  /// 12) apoptosis related parameters
                      /*26*/ sp.mean("t_apop_meas"),
                      /*27*/ sp.mean("t_duration_apoptosis")
        );

    time_step_d=tr.time_step_d;
    sim_duration_d=tr.sim_duration_d;
    trun_d=0;

return *this;
}

Cell_simulator::Cell_simulator(const Parameters& sp,
                               const Experiment& E):

    m(),
    APC(),
    NK(),
    LT(),

    time_step_d(),
    sim_duration_d(),
    trun_d(),
    filename(),


    experiment_(E),
    fitExperiment_(),
    prior_(sp),
    fitPar_()
{}


Results Cell_simulator::Simulate(const Parameters& simPar,
                                 const Treatment& tr,
                                 const Results& results)
{
    *this=applyParameters(simPar, tr);

    double Duratione=results.Duration();


    std::vector<Measurement> TNFs(results.TNF());
    std::size_t iTNFs=0;

    double tTNFs;
    if (!TNFs.empty())
        tTNFs=TNFs[iTNFs].Time();
    else
        tTNFs=Duratione+1;

    std::vector<Measurement> IFNs=results.IFN();
    std::size_t iIFNs=0;
    double tIFNs;
    if (!IFNs.empty())
        tIFNs=IFNs[iIFNs].Time();
    else
        tIFNs=Duratione+1;

    std::vector<Measurement> APC_exp=results.APC_expression();
    std::size_t iAPC_exp=0;
    double tAPC_exp;
    if (!APC_exp.empty())
        tAPC_exp=APC_exp[iAPC_exp].Time();
    else
        tAPC_exp=Duratione+1;

    std::vector<Measurement> NK_exp=results.NK_expression();
    std::size_t iNK_exp=0;
    double tNK_exp;
    if (!NK_exp.empty())
        tNK_exp=NK_exp[iNK_exp].Time();
    else
        tNK_exp=Duratione+1;

    std::vector<Measurement> LT_exp=results.LT_expression();
    std::size_t iLT_exp=0;
    double tLT_exp;
    if (!LT_exp.empty())
        tLT_exp=LT_exp[iLT_exp].Time();
    else tLT_exp=Duratione+1;

    std::vector<Measurement> APC_IFNs=results.APC_IFNg();
    std::size_t iAPC_IFN=0;
    double tAPC_IFN;
    if (!APC_IFNs.empty())
        tAPC_IFN=APC_IFNs[iAPC_IFN].Time();
    else tAPC_IFN=Duratione+1;

    std::vector<Measurement> APC_TNFs=results.APC_TNFa();
    std::size_t iAPC_TNF=0;
    double tAPC_TNF;
    if (!APC_TNFs.empty())
        tAPC_TNF=APC_TNFs[iAPC_TNF].Time();
    else tAPC_TNF=Duratione+1;

    std::vector<Measurement> NK_IFNs=results.NK_IFNg();
    std::size_t iNK_IFN=0;
    double tNK_IFN;
    if (!NK_IFNs.empty())
        tNK_IFN=NK_IFNs[iNK_IFN].Time();
    else tNK_IFN=Duratione+1;

    std::vector<Measurement> NK_TNFs=results.NK_TNFa();
    std::size_t iNK_TNF=0;
    double tNK_TNF;
    if (!NK_TNFs.empty())
        tNK_TNF=NK_TNFs[iNK_TNF].Time();
    else tNK_TNF=Duratione+1;

    std::vector<Measurement> LT_IFNs=results.LT_IFNg();
    std::size_t iLT_IFN=0;
    double tLT_IFN;
    if (!LT_IFNs.empty())
        tLT_IFN=LT_IFNs[iLT_IFN].Time();
    else tLT_IFN=Duratione+1;

    std::vector<Measurement> LT_TNFs=results.LT_TNFa();
    std::size_t iLT_TNF=0;
    double tLT_TNF;
    if (!LT_TNFs.empty())
        tLT_TNF=LT_TNFs[iLT_TNF].Time();
    else tLT_TNF=Duratione+1;

    std::vector<Measurement> LT_apops=results.LT_Apoptosis();
    std::size_t iLT_Apoptosis=0;
    double tLT_Apoptosis;
    if (!LT_apops.empty())
        tLT_Apoptosis=LT_apops[iLT_Apoptosis].Time();
    else tLT_Apoptosis=Duratione+1;

    std::vector<Measurement> Prols=results.Proliferation();
    std::size_t i_Proliferation=0;
    double t_Proliferation;
    if (!Prols.empty())
        t_Proliferation=Prols[i_Proliferation].Time();
    else t_Proliferation=Duratione+1;

    std::vector<Measurement> num_cellss=results.num_cells();
    std::size_t inum_cells=0;
    double t_num_cells;
    if (!num_cellss.empty())
        t_num_cells=num_cellss[inum_cells].Time();
    else t_num_cells=Duratione+1;

    double eps=1e-7;

    while (trun_d+eps<=results.Duration())
    {

        if(trun_d+eps>=tTNFs)
        {
            TNFs[iTNFs].setMeasurement(m.TNF());
            ++iTNFs;
            if (iTNFs<TNFs.size())
            {
                tTNFs=TNFs[iTNFs].Time();
            }
            else
            {
                tTNFs=results.Duration()+1;
            }
        };

        if(trun_d+eps>=tIFNs)
        {
            IFNs[iIFNs].setMeasurement(m.IFNgamma());
            ++iIFNs;
            if (iIFNs<IFNs.size())
            {
                tIFNs=IFNs[iIFNs].Time();
            }
            else
            {
                tIFNs=results.Duration()+1;
            }
        };


      if(trun_d+eps>=tAPC_exp)
        {
            APC_exp[iAPC_exp].setMeasurement(APC.percentage_cell_expressing_receptor());
            ++iAPC_exp;
            if (iAPC_exp<APC_exp.size())
            {
                tAPC_exp=APC_exp[iAPC_exp].Time();
            }
            else
            {
                tAPC_exp=results.Duration()+1;
            }
        };



        if(trun_d+eps>=tNK_exp)
        {
            NK_exp[iNK_exp].setMeasurement (
                        NK.percentage_NK_expressing_receptor());
            ++iNK_exp;
            if (iNK_exp<NK_exp.size())
            {
                tNK_exp=NK_exp[iNK_exp].Time();
            }
            else
            {
                tNK_exp=results.Duration()+1;
            }
        };


        if(trun_d+eps>=tLT_exp)
        {
            LT_exp[iLT_exp].setMeasurement (
                        LT.LT_percentage_cell_expressing_receptor());
            ++iLT_exp;
            if (iLT_exp<LT_exp.size())
            {
                tLT_exp=LT_exp[iLT_exp].Time();
            }
            else
            {
                tLT_exp=results.Duration()+1;
            }



        };

        if(trun_d+eps>=tAPC_IFN)
        {
            APC_IFNs[iAPC_IFN].setMeasurement (
                        APC.percentage_APC_producing_IFN());
            ++iAPC_IFN;
            if (iAPC_IFN<APC_IFNs.size())
            {
                tAPC_IFN=APC_IFNs[iAPC_IFN].Time();
            }
            else
            {
                tAPC_IFN=results.Duration()+1;
            }

          };

        if(trun_d+eps>=tAPC_TNF)
        {
            APC_TNFs[iAPC_TNF].setMeasurement (
                         APC.percentage_APC_producing_TNF());
            ++iAPC_TNF;
            if (iAPC_TNF<APC_TNFs.size())
            {
                tAPC_TNF=APC_TNFs[iAPC_TNF].Time();
            }
            else
            {
                tAPC_TNF=results.Duration()+1;
            }

             };

        if(trun_d+eps>=tNK_IFN)
        {
          NK_IFNs[iNK_IFN].setMeasurement (
          NK.percentage_NK_producing_IFN());
          ++iNK_IFN;
               if (iNK_IFN<NK_IFNs.size())
                  {
                    tNK_IFN=NK_IFNs[iNK_IFN].Time();
                  }
               else
                  {
                   tNK_IFN=results.Duration()+1;
                  }

                };

        if(trun_d+eps>=tNK_TNF)
        {
         NK_TNFs[iNK_TNF].setMeasurement (
         NK.percentage_NK_producing_TNF());
         ++iNK_TNF;
          if (iNK_TNF<NK_TNFs.size())
          {
           tNK_TNF=NK_TNFs[iNK_TNF].Time();
          }
          else
          {
          tNK_TNF=results.Duration()+1;
          }

         };

        if(trun_d+eps>=tLT_IFN)
        {
         LT_IFNs[iLT_IFN].setMeasurement (
         LT.percentage_LT_IFN_production());
         ++iLT_IFN;
         if (iLT_IFN<LT_IFNs.size())
            {
              tLT_IFN=LT_IFNs[iLT_IFN].Time();
            }
         else
            {
              tLT_IFN=results.Duration()+1;
            }

                    };

         if(trun_d+eps>=tLT_TNF)
         {
           LT_TNFs[iLT_TNF].setMeasurement (
           LT.percentage_LT_TNF_production());
           ++iLT_TNF;
           if (iLT_TNF<LT_TNFs.size())
              {
               tLT_TNF=LT_TNFs[iLT_TNF].Time();
              }
          else
              {
               tLT_TNF=results.Duration()+1;
              }

                    };

         if(trun_d+eps>=tLT_Apoptosis)
         {
          LT_apops[iLT_Apoptosis].setMeasurement (
          LT.percentage_apoptotic_LT_cells());
          ++iLT_Apoptosis;
          if (iLT_Apoptosis<LT_apops.size())
             {
               tLT_Apoptosis=LT_apops[iLT_Apoptosis].Time();
             }
          else
             {
               tLT_Apoptosis=results.Duration()+1;
             }

                     };

          if(trun_d+eps>=t_Proliferation)
          {
            Prols[i_Proliferation].setMeasurement (m.Tymidine_incorporated());
            ++i_Proliferation;
            if (i_Proliferation<Prols.size())
               {
                t_Proliferation=Prols[i_Proliferation].Time();
               }
           else
               {
                t_Proliferation=results.Duration()+1;
               }

                     };

          if(trun_d+eps>=t_num_cells)
          {
              num_cellss[inum_cells].setMeasurement (m.TNF());
              ++inum_cells;
              if (inum_cells<num_cellss.size())
              {
                  t_num_cells=num_cellss[inum_cells].Time();
              }
              else
              {
                  t_num_cells=results.Duration()+1;
              }
          };


        //char ch;

//	std::cout<<"\n media\n"<<m<<"\nLT\n"<<APC<<"\nNK\n"<<NK<<"\nLT\n"<<LT;
//	std::cout<<"\ntime \t"<<trun_d;
//	std::cin>>ch;
        APC.update(time_step_d,m,NK,LT);
        NK.update(time_step_d,m,APC,LT);
        LT.update(time_step_d,trun_d,m,APC,NK);
        m.update(time_step_d,trun_d, APC,NK,LT);
        trun_d+=time_step_d;

    }

    Results SimRes(TNFs,IFNs,APC_exp,NK_exp,LT_exp,APC_IFNs,APC_TNFs,NK_IFNs,NK_TNFs, LT_IFNs, LT_TNFs, LT_apops,Prols, num_cellss, Duratione);
    return SimRes;
}


Experiment Cell_simulator::Simulate(const Parameters& simPar,
                                    const Experiment& exp)
{
    Experiment sim;
    for (std::size_t i=0; i<exp.size(); i++)
    {
        Results r=Simulate(simPar,exp.Treatment_i(i),exp.Result_i(i));
        sim.push_back(exp.Treatment_i(i),
                      r);
    }

    return sim;

}

OptimizationResults Cell_simulator::Optimize(const Parameters& priorPar,
                                             const Parameters& initPar,
                                             const Experiment& experiment,
                                             double range,
                                             std::size_t numEval){}
