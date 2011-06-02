#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include "Cell_simulator.h"
#include "SimParameters.h"

void Cell_simulator::ask_parameters()
{
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
    std::cout<<"default is 5 ug, enter -1 to keep this value\n";
    std::cin>>Ag;
    if (Ag==-1)
        Ag=5;

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
        std::cin>>APC_max_proliferation_rate_;
    if (APC_max_proliferation_rate_==-1)
        APC_max_proliferation_rate_=1.0/240;

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
    std::cout<<"of exhausted APC";
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




    std::cout<<"Thaks";



    m=Media(max_num_cells_,
            Ag,
            Ag_internalization_rate,
            init_num_APC_cells+init_num_LT_cells+init_num_NK_cells,///num_cells
            TNF_deg);///Esto?TNF e IFN??

    APC=APC_cells(init_num_APC_cells,
                  APC_max_proliferation_rate_,
                  APC_no_to_free_rate_per_Ag_ ,
                  APC_free_to_bound_rate_per_LT_,
                  APC_exh_rate,
                  APC_IFN_free_prod_rate_,
                  APC_IFN_Ag_prod_rate_,
                  APC_IFN_bound_prod_rate_,
                  APC_TNF_free_prod_rate_,
                  APC_TNF_Ag_prod_rate_,
                  APC_TNF_bound_prod_rate_);

    NK=NK_cells(init_num_NK_cells,
                NK_max_proliferation_rate_,
                NK_no_to_free_rate_per_Ag_ ,
                NK_free_to_bound_rate_per_LT_,
                NK_exh_rate,
                NK_IFN_free_prod_rate_,
                NK_IFN_Ag_prod_rate_,
                NK_IFN_bound_prod_rate_,
                NK_TNF_free_prod_rate_,
                NK_TNF_Ag_prod_rate_,
                NK_TNF_bound_prod_rate_);

    LT=LT_cells(init_num_LT_cells,
                LT_num_specific,
                LT_max_no_receptor_prol_rate_,
                LT_max_free_prol_rate_,
                LT_max_bound_prol_rate_,
                LT_IFN_no_rec_prod_rate_,
                LT_IFN_free_prod_rate_,
                LT_IFN_bound_prod_rate_,
                LT_TNF_no_rec_prod_rate_,
                LT_TNF_free_prod_rate_,
                LT_TNF_bound_prod_rate_,
                LT_no_to_free_rate_per_APC_,
                LT_free_to_bound_rate_per_APC_);

};
void Cell_simulator::run()
{
    filename="out.txt";
    std::ofstream f;
    f.open(filename.c_str());
    trun_d=0;
    std::cout<<"round"<<"\t";
    std::cout<<"IFN"<<"\t";
    std::cout<<"TNF"<<"\t";
    std::cout<<"APC"<<"\t";
    std::cout<<"NK"<<"\t";
    std::cout<<"LT()"<<"\n";
//    std::cout<<"LT no receptor"<<"\t";
//    std::cout<<"LT free "<<"\t";
//    std::cout<<"LT bound"<<"\t";
//    std::cout<<"APC with Ag"<<"\t";
//    std::cout<<"APC bound"<<"\t";
//    std::cout<<"APC exhausted"<<"t";
//    std::cout<<"APC.IFNgamma_production_rate"<<"\t";
//    std::cout<<"LT.IFNgamma_production_rate"<<"\t";
//    std::cout<<"APC.TNF_production_rate"<<"\t";
//    std::cout<<"LT.TNF_production_rate"<<"\t";
//    std::cout<<"Ag"<<"\n",


    f<<"round"<<" , ";
    f<<"IFNamma[]"<<" , ";
    f<<"TNF[]"<<" , ";
    f<<"Total APC"<<" , ";
    f<<"Total NK"<<" , ";
    f<<"Total LT"<<" , ";
    f<<"APC with Ag "<<" , ";
    f<<"APC bound"<<" , ";
    f<<"APC exhausted"<<" , ";
    f<<"APC.IFNgamma_production_rate"<<" , ";
    f<<"APC.TNF_production_rate"<<" , ";
    f<<"NK with Ag "<<" , ";
    f<<"NK bound"<<" , ";
    f<<"NK exhausted"<<" , ";
    f<<"NK.IFNgamma_production_rate"<<" , ";
    f<<"NK.TNF_production_rate"<<" , ";
    f<<"LT no receptor"<<",";
    f<<"LT free "<<" , ";
    f<<"LT bound"<<" , ";
    f<<"LT.IFNgamma_production_rate"<<" , ";
    f<<"LT.TNF_production_rate"<<" , ";
    f<<"Ag"<<" , ";
    f<<"\n";


    while (trun_d<this->sim_duration_d)
    {
        if (trun_d-floor(trun_d)<time_step_d)
        {
            std::cout<<trun_d<<"\t";
            std::cout<<m.IFNgamma()<<"\t";
            std::cout<<m.TNF()<<"\t";
            std::cout<<APC.num()<<"\t";
            std::cout<<NK.num()<<"\t";
            std::cout<<LT.num()<<"\n";
//            std::cout<<APC.num_Ag()<<"\t";
//            std::cout<<APC.num_bound()<<"\t";
//            std::cout<<APC.num_exhausted()<<"t";
//            std::cout<<APC.IFNgamma_production_rate()<<"\t";
//            std::cout<<LT.num_cells_not_expressing_receptor()<<"\t";
//            std::cout<<LT.num_cells_expressing_receptor_and_free()<<"\t";
//            std::cout<<LT.num_cells_expressing_receptor_and_bound()<<"\t";
//            std::cout<<LT.IFNgamma_production_rate()<<"\t";
//            std::cout<<APC.TNF_production_rate()<<"\t";
//            std::cout<<LT.TNF_production_rate()<<"\t";
//            std::cout<<m.Ag()<<"\n";

            f<<trun_d<<" , ";
            f<<m.IFNgamma()<<" , ";
            f<<m.TNF()<<" , ";
            f<<APC.num()<<" , ";
            f<<NK.num()<<" , ";
            f<<LT.num()<<" , ";
            f<<APC.num_Ag()<<" , ";
            f<<APC.num_bound()<<" , ";
            f<<APC.num_exhausted()<<" , ";
            f<<APC.IFNgamma_production_rate()<<" , ";
            f<<APC.TNF_production_rate()<<" , ";
            f<<NK.NK_num_Ag()<<" , ";
            f<<NK.NK_num_bound()<<" , ";
            f<<NK.NK_exhausted()<<" , ";
            f<<NK.IFNgamma_production_rate()<<" , ";
            f<<NK.TNF_production_rate()<<" , ";
            f<<LT.num_cells_not_expressing_receptor()<<" , ";
            f<<LT.num_cells_expressing_receptor_and_free()<<" , ";
            f<<LT.num_cells_expressing_receptor_and_bound()<<" , ";
            f<<LT.IFNgamma_production_rate()<<" , ";
            f<<LT.TNF_production_rate()<<" , ";
            f<<m.Ag()<<"\n";

        };

        APC.update(time_step_d,m,NK,LT);
        NK.update(time_step_d,m,APC,LT);
        LT.update(time_step_d,m,APC,NK);
        m.update(time_step_d,APC,NK,LT);
        trun_d+=time_step_d;

    }
    f.close();
};

/// Media Main step of the simulation
/// it updates the state of the media and the cells populations for a given time period
void Media::update(double time_step,const APC_cells& APC_, const NK_cells& NK_, const LT_cells& LT_)
{
    /// IFN is increased by the production rate of each population;
    IFNgamma_d+=APC_.IFNgamma_production_rate()*time_step+
                NK_.IFNgamma_production_rate()*time_step +
                LT_.IFNgamma_production_rate()*time_step ;
    /// TNF is increased by the production rate of each population;
    TNF_d+=APC_.TNF_production_rate()*time_step+
           NK_.TNF_production_rate()*time_step +
           LT_.TNF_production_rate()*time_step -
           TNF_d*time_step*TNF_degradation();
    /// The total number of cells is the adittion of APC + NK + LT
    num_cells_d=APC_.num()+NK_.num()+LT_.num();
//    /// The Ag is internalized
//    Ag()-=time_step*APC_.num()*APC_.no_to_free_rate_per_Ag()*Ag_internalization_rate;
};

void APC_cells::update(double time_step,const Media& m, const NK_cells& NK, const LT_cells& LT)
{
    /// the proliferation is one for cero cells,
    /// 0 when the number of cells equals the maximum
    /// and negative (apoptosis) when there are more cells than the maximum
    double proliferation_ratio=(m.Max_num_cells()-m.num_cells())/m.Max_num_cells();

    /// we update each subpopulation of cells independently and we take into account the transition from one state to the other

    /// the number of free cells (no Ag) they proliferate according to the cell concentration (factor proliferation ratio)
    /// and some of them are "lost" since they internalize the Ag
    num_free_d+=num_free_d*time_step*
                (proliferation_ratio*APC_max_proliferation_rate_d - APC_no_to_free_rate_per_Ag_d*m.Ag());

    /** the cells that have internalize the Ag proliferate in the same way than the free
    they grow also by the free cells that internalize the Ag
    they shrink by the cells that interact with the LT cells*/
    num_Ag_d+=num_free_d*time_step*APC_no_to_free_rate_per_Ag_d*m.Ag()+
              num_Ag_d*time_step*
              (proliferation_ratio*APC_max_proliferation_rate_d
               -APC_free_to_bound_rate_per_LT_d*LT.num_cells_expressing_receptor());

    /// the cells that have interacted with LT grow accordingly with the number of cells that have internalized the Ag and the
    /// number of APC cells expressing the receptor and bound with LT cells
    num_LT_bound_d+=num_Ag_d*time_step*APC_free_to_bound_rate_per_LT_d*LT.num_cells_expressing_receptor()+
                    num_LT_bound_d*time_step*(proliferation_ratio*APC_max_proliferation_rate_d -
                            APC_exh_rate_d);

    /// the cells that have interacted with LT get exhausted
    num_exhausted_d+=num_LT_bound_d*APC_exh_rate_d*time_step+num_exhausted_d*time_step*proliferation_ratio*APC_max_proliferation_rate_d

                     ;
};
/// main step for the NK cells
void NK_cells::update(double time_step,const Media& m,const APC_cells& APC, const LT_cells& LT)
{
    /// the proliferation is one for cero cells,
    /// 0 when the number of cells equals the maximum
    /// and negative (apoptosis) when there are more cells than the maximum
    double proliferation_ratio=(m.Max_num_cells()-m.num_cells())/m.Max_num_cells();

    /// we update each subpopulation of cells independently and we take into account the transition from one state to the other

    /// the number of free cells (no Ag) they proliferate according to the cell concentration (factor proliferation ratio)
    /// and some of them are "lost" since they internalize the Ag
    NK_num_free_d+=NK_num_free_d*time_step*
                (proliferation_ratio*NK_max_proliferation_rate_d-NK_no_to_free_rate_per_Ag_d*m.Ag());

    /** the cells that have internalize the Ag proliferate in the same way than the free
    they grow also by the free cells that internalize the Ag
    they shrink by the cells that interact with the LT cells*/
    NK_num_Ag_d+=NK_num_free_d*time_step*NK_no_to_free_rate_per_Ag_d*m.Ag()+
              NK_num_Ag_d*time_step*
              (proliferation_ratio*NK_max_proliferation_rate_d
               -NK_free_to_bound_rate_per_LT_d*LT.num_cells_expressing_receptor());

    /// the cells that have interacted with LT grow accordingly with the number of cells that have internalized the Ag and the
    /// number of LT cells expressing the receptor
    NK_num_LT_bound_d+=NK_num_Ag_d*time_step*NK_free_to_bound_rate_per_LT_d*LT.num_cells_expressing_receptor()+
                    NK_num_LT_bound_d*time_step*(proliferation_ratio*NK_max_proliferation_rate_d -
                            NK_exh_rate_d);

    /// the cells that have interacted with LT get exhausted
    NK_num_exhausted_d+=NK_num_LT_bound_d*NK_exh_rate_d*time_step+NK_num_exhausted_d*time_step*proliferation_ratio*NK_max_proliferation_rate_d

                     ;
};
/// Main step for LT
void LT_cells::update(double time_step,const Media& m,const APC_cells& a,const NK_cells& NK)

{
    double proliferation_ratio=(m.Max_num_cells()-m.num_cells())/m.Max_num_cells();


    /// cells not sensitive to the Ag proliferate passively
    num_non_Agsp_d+=time_step*num_non_Agsp_d*proliferation_ratio*LT_max_no_receptor_prol_rate_d;


    /// Ag specific cells proliferate and som of them start to express the receptor
    num_Agsp_no_receptor_d+=time_step*num_Agsp_no_receptor_d*
                            (proliferation_ratio*LT_max_no_receptor_prol_rate_d-LT_no_to_free_rate_per_APC_d*a.num_Ag());


    num_Agsp_free_receptor_d+=time_step*num_Agsp_no_receptor_d*LT_no_to_free_rate_per_APC_d*a.num_Ag()+
                              time_step*num_Agsp_free_receptor_d*
                              (proliferation_ratio*LT_max_free_prol_rate_d-LT_free_to_bound_rate_per_APC_d*a.num_Ag());


    num_Agsp_bound_receptor_d+=time_step*num_Agsp_free_receptor_d*LT_free_to_bound_rate_per_APC_d*a.num_Ag()+
                               time_step*num_Agsp_bound_receptor_d*proliferation_ratio*LT_max_bound_prol_rate_d;

}






SimParameters::SimParameters():
    Ag(10),
    max_num_cells(2e6),
    init_num_LT_cells (9e5),
    LT_num_specific (1000),
    init_num_APC_cells (1e5),
    init_num_NK_cells(1e5),
    sim_duration_d (120),
    time_step_d (0.01),
    APC_max_proliferation_rate_ (1.0/240),
    LT_max_no_receptor_prol_rate_(1.0/96),
    LT_max_free_prol_rate_(1.0/2),
    LT_max_bound_prol_rate_(1.0/1.2),
    APC_no_to_free_rate_per_Ag_(1.0/30),
    APC_free_to_bound_rate_per_LT_(1.0/1e5),
    LT_no_to_free_rate_per_APC_(1.0/6e5),
    LT_free_to_bound_rate_per_APC_(1.0/1e5),
    APC_exh_rate(1/1e5),
    LT_IFN_no_rec_prod_rate_(0.001/1e5),
    LT_IFN_free_prod_rate_(101.0/1e5),
    LT_IFN_bound_prod_rate_(200.0/1e5),
    LT_TNF_no_rec_prod_rate_(0.001/1e5),
    LT_TNF_free_prod_rate_(10.0/1e5),
    LT_TNF_bound_prod_rate_(20.0/1e5),
    APC_IFN_free_prod_rate_(0.5/1e5),
    APC_IFN_Ag_prod_rate_(5.0/1e5),
    APC_IFN_bound_prod_rate_(10.0/1e5),
    APC_TNF_free_prod_rate_(5/1e5),
    APC_TNF_Ag_prod_rate_(570/1e5),
    APC_TNF_bound_prod_rate_(1110/1e5),
    NK_IFN_free_prod_rate_(0.5/1e5),
    NK_IFN_Ag_prod_rate_(5.0/1e5),
    NK_IFN_bound_prod_rate_(10.0/1e5),
    NK_TNF_free_prod_rate_(5/1e5),
    NK_TNF_Ag_prod_rate_(570/1e5),
    NK_TNF_bound_prod_rate_(1110/1e5),
    Ag_internalization_rate (0.5),
    TNF_deg (0.5)
    {};



Cell_simulator::Cell_simulator(const SimParameters& sp):
    sim_duration_d(sp.sim_duration_d),
    time_step_d(sp.time_step_d),
    m(sp.max_num_cells_,sp.init_num_APC_cells+sp.init_num_LT_cells+sp.init_num_NK_cells,sp.Ag,sp.TNF_deg,Ag_internalization_rate,0.0),
    APC(APC_cells(sp.init_num_APC_cells,
                  sp.APC_max_proliferation_rate_,
                  sp.APC_no_to_free_rate_per_Ag_ ,
                  sp.APC_free_to_bound_rate_per_LT_,
                  sp.APC_exh_rate,
                  sp.APC_IFN_free_prod_rate_,
                  sp.APC_IFN_Ag_prod_rate_,
                  sp.APC_IFN_bound_prod_rate_,
                  sp.APC_TNF_free_prod_rate_,
                  sp.APC_TNF_Ag_prod_rate_,
                  sp.APC_TNF_bound_prod_rate_)),

     NK(NK_cells (sp.init_num_NK_cells,
                  sp.NK_max_proliferation_rate_,
                  sp.NK_no_to_free_rate_per_Ag_ ,
                  sp.NK_free_to_bound_rate_per_LT_,
                  sp.NK_exh_rate,
                  sp.NK_IFN_free_prod_rate_,
                  sp.NK_IFN_Ag_prod_rate_,
                  sp.NK_IFN_bound_prod_rate_,
                  sp.NK_TNF_free_prod_rate_,
                  sp.NK_TNF_Ag_prod_rate_,
                  sp.NK_TNF_bound_prod_rate_)),

    LT(LT_cells (sp.init_num_LT_cells,
                 sp.LT_num_specific,
                 sp.LT_max_no_receptor_prol_rate_,
                 sp.LT_max_free_prol_rate_,
                 sp.LT_max_bound_prol_rate_,
                 sp.LT_IFN_no_rec_prod_rate_,
                 sp.LT_IFN_free_prod_rate_,
                 sp.LT_IFN_bound_prod_rate_,
                 sp.LT_TNF_no_rec_prod_rate_,
                 sp.LT_TNF_free_prod_rate_,
                 sp.LT_TNF_bound_prod_rate_,
                 sp.LT_no_to_free_rate_per_APC_,
                 sp.LT_free_to_bound_rate_per_APC_)) {};

