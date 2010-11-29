#include "Cell_simulator.h"
#include <iostream>
#include <string>
#include <cmath>
#include <fstream>

void Cell_simulator::ask_parameters()
{
   // std::cout<<"enter the filename where you want to store the values";
   // std::cin>>this->filename;
   filename="out.txt";

    double max_num_cells_;
  // std::cout<<"enter the number of cells that the media supports\n";
  //  std::cout<<"defaut is 2e6, enter 0 to keep this value\n";
  //  std::cin>>max_num_cells_;
  //  if (max_num_cells_==0)
        max_num_cells_=2e6;


    double init_num_LT_cells;
    std::cout<<"enter the intialnumber of LT cells\n";
    std::cout<<"defaut is 9e5, enter -1 to keep this value\n";
    std::cin>>init_num_LT_cells;
    if (init_num_LT_cells==-1)
        init_num_LT_cells=9e5;

    double num_specific;
    std::cout<<"enter the intial number of LT cells specific for the antigen \n";
    std::cout<<"defaut is 1000, enter -1  to keep this value\n";
   std::cin>>num_specific;
   if (num_specific==-1)
        num_specific=1000;


    double init_num_APC_cells;
    std::cout<<"enter the intialnumber of APC cells\n";
    std::cout<<"defaut is 1e5, enter -1 to keep this value\n";
    std::cin>>init_num_APC_cells;
    if (init_num_APC_cells==-1)
        init_num_APC_cells=1e5;


    double AG;
    std::cout<<"enter the applied concentration of AG\n";
    std::cout<<"default is 5 ug, enter -1 to keep this value\n";
    std::cin>>AG;
    if (AG==-1)
        AG=5;



    //std::cout<<"enter the duration of the simulation in hours\n";
    //std::cout<<"default is 120, enter 0 to keep this value\n";
    //std::cin>>this->sim_duration_d;
    //if (sim_duration_d==0)
      sim_duration_d=120;


    //std::cout<<"enter the duration of the time step of the simulation in hours\n";
    //std::cout<<"defaut is 0.01, enter 0 to keep this value\n";
    //std::cin>>this->time_step_d;
    //if(time_step_d==0)
    time_step_d=0.01;

    double max_proliferation_rate_, max_no_receptor_prol_rate_;

    double max_free_prol_rate_;
    double max_bound_prol_rate_;
    //std::cout<<"enter the following proliferation rates\n";
   // std::cout<<"of APC cells default value =1/240 h\n";
   // std::cin>>max_proliferation_rate_;
   // if (max_proliferation_rate_==0)
   max_proliferation_rate_=1.0/240;
    //std::cout<<"of LT cells that do not express the receptor, default value  1/96 1/h \n";
    //std::cin>>max_no_receptor_prol_rate_;
    //if (max_no_receptor_prol_rate_==0)
    max_no_receptor_prol_rate_=1.0/96;

    //std::cout<<"of LT cells that do express the receptor and it it is free, default value 1/2 \n";
    //std::cin>>max_free_prol_rate_;
   // if (max_free_prol_rate_==0)
   max_free_prol_rate_=1.0/2;

    //std::cout<<"of LT cells where the receptor have been bound, default value 1/1.2 \n";
   // std::cin>>max_bound_prol_rate_;
   // if (max_bound_prol_rate_==0)
   max_bound_prol_rate_=1.0/1.2;

    double no_to_free_rate_per_AG_ ,free_to_bound_rate_per_LT_;
   // std::cout<<"\n\n enter the following conversion rates in APC cells\n";
   // std::cout<<"of Antigen internalization default it needs 30 h for 1 ug/ml of antigen \n";
//    std::cin>>no_to_free_rate_per_AG_;
   // if (no_to_free_rate_per_AG_==0)
   no_to_free_rate_per_AG_=1.0/30;

   // std::cout<<"of lingand receptor binding default value it needs 1 h for a population of 100e3 \n";
   // std::cin>>free_to_bound_rate_per_LT_;
   // if (free_to_bound_rate_per_LT_==0)
    free_to_bound_rate_per_LT_=1.0/1e5;

    double no_to_free_rate_per_APC_;
    double free_to_bound_rate_per_APC_;
    //std::cout<<"\n\n enter the following conversion rates in LT cells\n";
    //std::cout<<"of receptor expression per APC cell default value is 6 hours for 1e5 cells\n";
    //std::cin>>no_to_free_rate_per_APC_;
   // if (no_to_free_rate_per_APC_==0)
    no_to_free_rate_per_APC_=1.0/6e5;
    //std::cout<<"of lingand receptor binding per APC cell  default value is 1h for a population of 1e5 cells\n";
    //std::cin>>free_to_bound_rate_per_APC_;
   // if (free_to_bound_rate_per_APC_==0)
   free_to_bound_rate_per_APC_=1.0/1e5;

    double LT_IFN_no_rec_prod_rate_;
    double LT_IFN_free_prod_rate_;
    double LT_IFN_bound_prod_rate_;
  //  std::cout<<"\n\n enter the following IFN production rates in LT cells\n";
   // std::cout<<"of cells without receptor default value is 0.001 pg per hour per 1e5 cells \n";
   // std::cin>>LT_IFN_no_rec_prod_rate_;
   // if (LT_IFN_no_rec_prod_rate_==0)
   LT_IFN_no_rec_prod_rate_=0.001/1e5;

   // std::cout<<"of cells with free receptor  default value is 500 pg per hour per 1e5 cells\n";
   // std::cin>>LT_IFN_free_prod_rate_;
    //if (LT_IFN_free_prod_rate_==0)
    LT_IFN_free_prod_rate_=101.0/1e5;


    //std::cout<<"of cells with bound receptor,  default value is 1000 pg per hour per 1e5 cells \n";
    //std::cin>>LT_IFN_bound_prod_rate_;
    //if (LT_IFN_bound_prod_rate_==0)
    LT_IFN_bound_prod_rate_=200.0/1e5;

    double LT_TNF_no_rec_prod_rate_;
    double LT_TNF_free_prod_rate_;
    double LT_TNF_bound_prod_rate_;
  //  std::cout<<"\n\n enter the following IFN production rates in LT cells\n";
   // std::cout<<"of cells without receptor default value is 0.001 pg per hour per 1e5 cells \n";
   // std::cin>>LT_IFN_no_rec_prod_rate_;
   // if (LT_IFN_no_rec_prod_rate_==0)
   LT_TNF_no_rec_prod_rate_=0.001/1e5;

   // std::cout<<"of cells with free receptor  default value is 500 pg per hour per 1e5 cells\n";
   // std::cin>>LT_IFN_free_prod_rate_;
    //if (LT_IFN_free_prod_rate_==0)
    LT_TNF_free_prod_rate_=10.0/1e5;


    //std::cout<<"of cells with bound receptor,  default value is 1000 pg per hour per 1e5 cells \n";
    //std::cin>>LT_IFN_bound_prod_rate_;
    //if (LT_IFN_bound_prod_rate_==0)
    LT_TNF_bound_prod_rate_=20.0/1e5;

    double APC_IFN_free_prod_rate_,
    APC_IFN_AG_prod_rate_,
    APC_IFN_bound_prod_rate_;
    //std::cout<<"\n\n enter the following IFN production rates in LT cells\n";
    //std::cout<<"of cells without receptor or free receptor default value is 0.5 pg per hour per 1e5 cells \n";
    //std::cin>>APC_IFN_free_prod_rate_;
    //if (APC_IFN_free_prod_rate_==0)
    APC_IFN_free_prod_rate_=0.5/1e5;

    //std::cout<<"of cells with the internalized antigen  default value is 5 pg per hour per 1e5 cells  \n";
    //std::cin>>APC_IFN_AG_prod_rate_;
    //if (APC_IFN_AG_prod_rate_==0)
    APC_IFN_AG_prod_rate_=5.0/1e5;

    //std::cout<<"of cells with bound receptor  default value is 10 pg per hour per 1e5 cells  \n";
    //std::cin>>APC_IFN_bound_prod_rate_;
    //if (APC_IFN_bound_prod_rate_==0)
    APC_IFN_bound_prod_rate_=10.0/1e5;

double APC_TNF_free_prod_rate_,
    APC_TNF_AG_prod_rate_,
    APC_TNF_bound_prod_rate_;
    //std::cout<<"\n\n enter the following IFN production rates in LT cells\n";
    //std::cout<<"of cells without receptor or free receptor default value is 0.5 pg per hour per 1e5 cells \n";
    //std::cin>>APC_IFN_free_prod_rate_;
    //if (APC_IFN_free_prod_rate_==0)
    APC_TNF_free_prod_rate_=5/1e5;

    //std::cout<<"of cells with the internalized antigen  default value is 5 pg per hour per 1e5 cells  \n";
    //std::cin>>APC_IFN_AG_prod_rate_;
    //if (APC_IFN_AG_prod_rate_==0)
    APC_TNF_AG_prod_rate_=570/1e5;

    //std::cout<<"of cells with bound receptor  default value is 10 pg per hour per 1e5 cells  \n";
    //std::cin>>APC_IFN_bound_prod_rate_;
    //if (APC_IFN_bound_prod_rate_==0)
    APC_TNF_bound_prod_rate_=1110/1e5;




    m=Media(max_num_cells_,init_num_APC_cells+init_num_LT_cells,AG,0.0);

    APC=APC_cells(init_num_APC_cells,
                  max_proliferation_rate_,
                  no_to_free_rate_per_AG_ ,
                  free_to_bound_rate_per_LT_,
                  APC_IFN_free_prod_rate_,
                  APC_IFN_AG_prod_rate_,
                  APC_IFN_bound_prod_rate_,
                  APC_TNF_free_prod_rate_,
                  APC_TNF_AG_prod_rate_,
                  APC_TNF_bound_prod_rate_);

    LT=LT_cells( init_num_LT_cells,
                 num_specific,
                 max_no_receptor_prol_rate_,
                 max_free_prol_rate_,
                 max_bound_prol_rate_,
                 LT_IFN_no_rec_prod_rate_,
                 LT_IFN_free_prod_rate_,
                 LT_IFN_bound_prod_rate_,
                 LT_TNF_no_rec_prod_rate_,
                 LT_TNF_free_prod_rate_,
                 LT_TNF_bound_prod_rate_,
                 no_to_free_rate_per_APC_,
                 free_to_bound_rate_per_APC_);

};
void Cell_simulator::run()
{
    std::ofstream f;
    f.open(filename.c_str());
    trun_d=0;
    std::cout<<"trun_d"<<"\t"<<"m.IFNgamma()"<<"\t"<<"m.TNF()"<<"\t"<<"LT.num()"<<"\t";
    std::cout<<"no receptor"<<"\t";
    std::cout<<"free "<<"\t";
    std::cout<<"bound"<<"\t";
    std::cout<<"APC"<<"\n";
    std::cout<<"AG"<<"\n";
    std::cout<<"bound"<<"\n";
    std::cout<<"APC.IFNgamma_production_rate"<<"\n";
    std::cout<<"LT.IFNgamma_production_rate"<<"\n";
    std::cout<<"APC.TNF_production_rate"<<"\n";
    std::cout<<"LT.TNF_production_rate"<<"\n"

,
    f<<"trun_d"<<" , "<<"m.IFNamma()"<<" , "<<"m.TNF"<<" , "<<"LT.num()"<<" , ";
    f<<"no receptor"<<",";
    f<<"free "<<" , ";
    f<<"bound"<<" , ";
    f<<"APC"<<" , ";
    f<<"AG "<<" , ";
    f<<"bound"<<" , ";
    f<<"APC.IFNgamma_production_rate"<<" , ";
    f<<"LT.IFNgamma_production_rate"<<" , ";
    f<<"APC.TNF_production_rate"<<" , ";
    f<<"LT.TNF_production_rate"<<" , ";
    f<<"\n";

    while (trun_d<this->sim_duration_d)
    {
        if (trun_d-floor(trun_d)<time_step_d)
        {
            std::cout<<trun_d<<"\t"<<m.IFNgamma()<<m.TNF()<<"\t"<<LT.num()<<"\t";
            std::cout<<LT.num_cells_not_expressing_receptor()<<"\t";
            std::cout<<LT.num_cells_expressing_receptor_and_free()<<"\t";
            std::cout<<LT.num_cells_expressing_receptor_and_bound()<<"\t";
            std::cout<<APC.num()<<"\t";
            std::cout<<APC.num_AG()<<"\t";
            std::cout<<APC.num_bound()<<"\t";
            std::cout<<APC.IFNgamma_production_rate()<<"\t";
            std::cout<<LT.IFNgamma_production_rate()<<"\t";
            std::cout<<APC.TNF_production_rate()<<"\t";
            std::cout<<LT.TNF_production_rate()<<"\n";

       f<<trun_d<<" , "<<m.IFNgamma()<<" , "<<m.TNF()<<" , "<<LT.num()<<" , ";
            f<<LT.num_cells_not_expressing_receptor()<<" , ";
            f<<LT.num_cells_expressing_receptor_and_free()<<" , ";
            f<<LT.num_cells_expressing_receptor_and_bound()<<" , ";
            f<<APC.num()<<",";
            f<<APC.num_AG()<<",";
            f<<APC.num_bound()<<",";
            f<<APC.IFNgamma_production_rate()<<",";
            f<<LT.IFNgamma_production_rate()<<",";
            f<<APC.TNF_production_rate()<<",";
            f<<LT.TNF_production_rate()<<"\n";

        };
        APC.update(time_step_d,m,LT);
        LT.update(time_step_d,m,APC);
        m.update(time_step_d,APC,LT);
        trun_d+=time_step_d;

    }
    f.close();
};
//void Cell_simulator::show_results(){};





/// main step of the simulation
/// it updates the state of the media and the cells populations for a given time period
void Media::update(double time_step,const APC_cells& APC_,const LT_cells& LT_)
{
    /// IFN is increased by the production rate of each population;
    IFNgamma_d+=LT_.IFNgamma_production_rate()*time_step+APC_.IFNgamma_production_rate()*time_step;
    num_cells_d=APC_.num()+LT_.num();
    /// TNF is increased by the production rate of each population;
    TNF_d+=LT_.TNF_production_rate()*time_step+APC_.TNF_production_rate()*time_step;
    num_cells_d=APC_.num()+LT_.num();

};


/// main step for the APC cells
void APC_cells::update(double time_step,const Media& m, const LT_cells& LT)
{
    /// the proliferation is one for cero cells,
    /// 0 when the number of cells equals the maximum
    /// and negative (apoptosis) when there are more cells than the maximum
    double proliferation_ratio=(m.Max_num_cells()-m.num_cells())/m.Max_num_cells();

    /// we update each subpopulation of cells independently and we take into account the transition from one state to the other

    /// the number of free cells (no AG) they proliferate according to the cell concentration (factor proliferation ratio)
    /// and some of them are "lost" since they internalize the AG
    num_free_d+=num_free_d*time_step*
                (proliferation_ratio*max_proliferation_rate_d- no_to_free_rate_per_AG_d*m.AG());


    /** the cells that have internalize the AG proliferate in the same way than the free
    they grow also by the free cells that internalize the AG
    they shrink by the cells that interact with the LT cells

    */
    num_AG_d+=num_free_d*time_step*no_to_free_rate_per_AG_d*m.AG()+
              num_AG_d*time_step*
              (proliferation_ratio*max_proliferation_rate_d- free_to_bound_rate_per_LT_d*LT.num_cells_expressing_receptor());


    /// the cells that have interacted with LT grow accordingly with the number of cells that have internalized the AG and the
    /// number of LT cells expressing the receptor
    num_LT_bound_d+=num_AG_d*time_step*free_to_bound_rate_per_LT_d*LT.num_cells_expressing_receptor()+
                    num_LT_bound_d*time_step*proliferation_ratio*max_proliferation_rate_d;


};

void LT_cells::update(double time_step,const Media& m, const APC_cells& a)
{
    double proliferation_ratio=(m.Max_num_cells()-m.num_cells())/m.Max_num_cells();


    /// cells not sensitive to the AG proliferate passively
    num_non_AGsp_d+=time_step*num_non_AGsp_d*proliferation_ratio*max_no_receptor_prol_rate_d;


    /// AG specific cells proliferate and som of them start to express the receptor
    num_AGsp_no_receptor_d+=time_step*num_AGsp_no_receptor_d*
                            (proliferation_ratio*max_no_receptor_prol_rate_d-no_to_free_rate_per_APC_d*a.num_AG());


    num_AGsp_free_receptor_d+=time_step*num_AGsp_no_receptor_d*no_to_free_rate_per_APC_d*a.num_AG()+
                              time_step*num_AGsp_free_receptor_d*
                              (proliferation_ratio*max_free_prol_rate_d-free_to_bound_rate_per_APC_d*a.num_AG());


    num_AGsp_bound_receptor_d+=time_step*num_AGsp_free_receptor_d*free_to_bound_rate_per_APC_d*a.num_AG()+
                               time_step*num_AGsp_bound_receptor_d*proliferation_ratio*max_bound_prol_rate_d;

};
