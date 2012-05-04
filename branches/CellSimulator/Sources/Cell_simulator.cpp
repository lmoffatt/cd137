#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <cstddef>
#include <vector>
#include <sstream>
#include "Results.h"
#include "Includes/Cell_simulator.h"
#include "Includes/SimParameters.h"
#include "Includes/Media.h"
#include "Includes/APC.h"
#include "Includes/NK.h"
#include "Includes/LT.h"
#include "Includes/LevenbergMarquardt.h"
#include "Includes/OptimizationResults.h"
#include "Includes/BayesIteration.h"

std::ostream& Cell_simulator::run(std::ostream& f)
{
    f<<"Experiment run \n";

    trun_d=0;
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
    /*12*/    f<<"%APC expresing receptor"<<" , ";
    /*13*/    f<<"APC.IFNgamma_production_rate"<<" , ";
    /*14*/    f<<"APC.percentage of IFN producing cells"<<" , ";
    /*15*/    f<<"APC.TNF_production_rate"<<" , ";
    /*16*/    f<<"APC.percentage of TNF producing cells"<<" , ";

    /*17*/    f<<"NK0"<<" , ";
    /*18*/    f<<"NKa"<<" , ";
    /*19*/    f<<"NKbo"<<" , ";
    /*20*/    f<<"NKbo_bl"<<" , ";
    /*21*/    f<<"NKbl"<<" , ";
    /*22*/    f<<"%NK expresing receptor"<<" , ";
    /*23*/    f<<"NK.IFNgamma_production_rate"<<" , ";
    /*24*/    f<<"NK.percentage of IFN producing cell"<<" , ";
    /*25*/    f<<"NK.TNF_production_rate"<<" , ";
    /*26*/    f<<"NK.percentage of TNF producing cell"<<" , ";

    /*27*/    f<<"LT no Agsp"<<" , ";
    /*28*/    f<<"LT0"<<" , ";
    /*29*/    f<<"LTbo"<<" , ";
    /*30*/    f<<"LTbl"<<" , ";
    /*31*/    f<<"%LT expresing receptor"<<" , ";
    /*32*/    f<<"LT.IFNgamma_production_rate"<<" , ";
    /*33*/    f<<"LT.percentage of IFN producing cell"<<" , ";
    /*34*/    f<<"LT.TNF_production_rate"<<" , ";
    /*35*/    f<<"LT.percentage of TNF producing cell"<<" , ";
    /*36*/    f<<"LT undergoing apoptosis"<<" , ";
    /*37*/    f<<"Tymidine incorporated"<<" , ";

    /*38*/    f<<"Ag"<<" , ";
    /*39*/    f<<"Ab"<<" , ";
    f<<"\n";

    double eps=1e-7;
//    double tstart=0;
//    double tend=1;


    std::stringstream ss;
    while ((trun_d<this->sim_duration_d))
    {

        if ((trun_d+eps-floor(trun_d+eps)<time_step_d))/*||((trun_d>tstart)&&(trun_d<tend)))*/
        {

            if (m.TNF()==m.TNF())
                ss.str()="";
            else
                f<<ss.str();

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



            /*12*/           f<<APC.percentage_cell_expressing_receptor()<<" , ";
            /*13*/           f<<APC.APC_IFNgamma_production_rate()<<" , ";
            /*14*/           f<<APC.percentage_APC_producing_IFN()<<" , ";
            /*15*/           f<<APC.APC_TNF_production_rate()<<" , ";
            /*16*/           f<<APC.percentage_APC_producing_TNF()<<" , ";


            /*17*/           f<<NK.NK0()<<" , ";
            /*18*/           f<<NK.NKa()<<" , ";
            /*19*/           f<<NK.NKbo()<<" , ";
            /*20*/           f<<NK.NKbo_Ab()<<" , ";
            /*21*/           f<<NK.NKbl()<<" , ";
            /*22*/           f<<NK.percentage_NK_expressing_receptor()<<" , ";
            /*23*/           f<<NK.NK_IFNgamma_production_rate()<<" , ";
            /*24*/           f<<NK.percentage_NK_producing_IFN()<<" , ";
            /*25*/           f<<NK.NK_TNF_production_rate()<<" , ";
            /*26*/           f<<NK.percentage_NK_producing_TNF()<<" , ";

            /*27*/    f<<LT.LTns()<<" , ";
            /*28*/    f<<LT.LT0()<<" , ";
            /*29*/    f<<LT.LTbo()<<" , ";
            /*30*/    f<<LT.LTbl()<<" , ";

            /*31*/    f<<LT.LT_percentage_cell_expressing_receptor()<<" , ";
            /*32*/    f<<LT.LT_IFNgamma_production_rate()<<" , ";
            /*33*/    f<<LT.percentage_LT_IFN_production()<<" , ";
            /*34*/    f<<LT.TNF_production_rate()<<" , ";
            /*35*/    f<<LT.percentage_LT_TNF_production()<<" , ";
            /*36*/    f<<LT.percentage_apoptotic_LT_cells()<<" , ";
            /*37*/    f<<m.Tymidine_incorporated()<<" , ";

            /*38*/    f<<m.Ag()<<" , ";
            /*39*/    f<<m.Ab()<<"\n";


        }
        else
        {

            /*1*/            ss<<trun_d<<" , ";
            /*2*/            ss<<m.IFNgamma()<<" , ";
            /*3*/            ss<<m.TNF()<<" , ";

            /*4*/            ss<<APC.num_APC()<<" , ";
            /*5*/            ss<<NK.num_NK()<<" , ";
            /*6*/            ss<<LT.num_LT()<<" , ";

            /*7*/            ss<<APC.APC0()<<" , ";
            /*8*/            ss<<APC.APCa()<<" , ";
            /*9*/            ss<<APC.APCbo()<<" , ";
            /*10*/           ss<<APC.APCbo_Ab()<<" , ";
            /*11*/           ss<<APC.APCbl()<<" , ";



            /*12*/           ss<<APC.percentage_cell_expressing_receptor()<<" , ";
            /*13*/           ss<<APC.APC_IFNgamma_production_rate()<<" , ";
            /*14*/           ss<<APC.percentage_APC_producing_IFN()<<" , ";
            /*15*/           ss<<APC.APC_TNF_production_rate()<<" , ";
            /*16*/           ss<<APC.percentage_APC_producing_TNF()<<" , ";


            /*17*/           ss<<NK.NK0()<<" , ";
            /*18*/           ss<<NK.NKa()<<" , ";
            /*19*/           ss<<NK.NKbo()<<" , ";
            /*20*/           ss<<NK.NKbo_Ab()<<" , ";
            /*21*/           ss<<NK.NKbl()<<" , ";
            /*22*/           ss<<NK.percentage_NK_expressing_receptor()<<" , ";
            /*23*/           ss<<NK.NK_IFNgamma_production_rate()<<" , ";
            /*24*/           ss<<NK.percentage_NK_producing_IFN()<<" , ";
            /*25*/           ss<<NK.NK_TNF_production_rate()<<" , ";
            /*26*/           ss<<NK.percentage_NK_producing_TNF()<<" , ";

            /*27*/    ss<<LT.LTns()<<" , ";
            /*28*/    ss<<LT.LT0()<<" , ";
            /*29*/    ss<<LT.LTbo()<<" , ";
            /*30*/    ss<<LT.LTbl()<<" , ";

            /*31*/    ss<<LT.LT_percentage_cell_expressing_receptor()<<" , ";
            /*32*/    ss<<LT.LT_IFNgamma_production_rate()<<" , ";
            /*33*/    ss<<LT.percentage_LT_IFN_production()<<" , ";
            /*34*/    ss<<LT.TNF_production_rate()<<" , ";
            /*35*/    ss<<LT.percentage_LT_TNF_production()<<" , ";
            /*36*/    ss<<LT.percentage_apoptotic_LT_cells()<<" , ";
            /*37*/    ss<<m.Tymidine_incorporated()<<" , ";

            /*38*/    ss<<m.Ag()<<" , ";
            /*39*/    ss<<m.Ab()<<"\n";

        }
        if (m.TNF()!=m.TNF())
        {
            f<<ss.str()<<std::endl;
            break;
        }


        APC.update(time_step_d,m,NK,LT);
        NK.update(time_step_d,m,APC,LT);
        LT.update(time_step_d,trun_d,m,APC,NK);
        m.update(time_step_d,trun_d,APC,NK,LT);
        trun_d+=time_step_d;




    }
}


std::ostream& Cell_simulator::run(std::ostream& s, const Parameters& par) const
{
    Cell_simulator cell(*this);

    for (std::size_t i=0; i<this->experiment_.size(); i++)
    {
        cell.applyParameters(par,this->experiment_.Treatment_i(i));
        cell.run(s);
    }
    return s;
}

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
    /*12*/    f<<"%APC expresing receptor"<<" , ";
    /*13*/    f<<"APC.IFNgamma_production_rate"<<" , ";
    /*14*/    f<<"APC.percentage of IFN producing cells"<<" , ";
    /*15*/    f<<"APC.TNF_production_rate"<<" , ";
    /*16*/    f<<"APC.percentage of TNF producing cells"<<" , ";

    /*17*/    f<<"NK0"<<" , ";
    /*18*/    f<<"NKa"<<" , ";
    /*19*/    f<<"NKbo"<<" , ";
    /*20*/    f<<"NKbo_bl"<<" , ";
    /*21*/    f<<"NKbl"<<" , ";
    /*22*/    f<<"%NK expresing receptor"<<" , ";
    /*23*/    f<<"NK.IFNgamma_production_rate"<<" , ";
    /*24*/    f<<"NK.percentage of IFN producing cell"<<" , ";
    /*25*/    f<<"NK.TNF_production_rate"<<" , ";
    /*26*/    f<<"NK.percentage of TNF producing cell"<<" , ";

    /*27*/    f<<"LT no Agsp"<<" , ";
    /*28*/    f<<"LT0"<<" , ";
    /*29*/    f<<"LTbo"<<" , ";
    /*30*/    f<<"LTbl"<<" , ";
    /*31*/    f<<"%LT expresing receptor"<<" , ";
    /*32*/    f<<"LT.IFNgamma_production_rate"<<" , ";
    /*33*/    f<<"LT.percentage of IFN producing cell"<<" , ";
    /*34*/    f<<"LT.TNF_production_rate"<<" , ";
    /*35*/    f<<"LT.percentage of TNF producing cell"<<" , ";
    /*36*/    f<<"LT undergoing apoptosis"<<" , ";
    /*37*/    f<<"Tymidine incorporated"<<" , ";

    /*38*/    f<<"Ag"<<" , ";
    /*39*/    f<<"Ab"<<" , ";
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



            /*12*/           f<<APC.percentage_cell_expressing_receptor()<<" , ";
            /*13*/           f<<APC.APC_IFNgamma_production_rate()<<" , ";
            /*14*/           f<<APC.percentage_APC_producing_IFN()<<" , ";
            /*15*/           f<<APC.APC_TNF_production_rate()<<" , ";
            /*16*/           f<<APC.percentage_APC_producing_TNF()<<" , ";


            /*17*/           f<<NK.NK0()<<" , ";
            /*18*/           f<<NK.NKa()<<" , ";
            /*19*/           f<<NK.NKbo()<<" , ";
            /*20*/           f<<NK.NKbo_Ab()<<" , ";
            /*21*/           f<<NK.NKbl()<<" , ";
            /*22*/           f<<NK.percentage_NK_expressing_receptor()<<" , ";
            /*23*/           f<<NK.NK_IFNgamma_production_rate()<<" , ";
            /*24*/           f<<NK.percentage_NK_producing_IFN()<<" , ";
            /*25*/           f<<NK.NK_TNF_production_rate()<<" , ";
            /*26*/           f<<NK.percentage_NK_producing_TNF()<<" , ";

            /*27*/    f<<LT.LTns()<<" , ";
            /*28*/    f<<LT.LT0()<<" , ";
            /*29*/    f<<LT.LTbo()<<" , ";
            /*30*/    f<<LT.LTbl()<<" , ";

            /*31*/    f<<LT.LT_percentage_cell_expressing_receptor()<<" , ";
            /*32*/    f<<LT.LT_IFNgamma_production_rate()<<" , ";
            /*33*/    f<<LT.percentage_LT_IFN_production()<<" , ";
            /*34*/    f<<LT.TNF_production_rate()<<" , ";
            /*35*/    f<<LT.percentage_LT_TNF_production()<<" , ";
            /*36*/    f<<LT.percentage_apoptotic_LT_cells()<<" , ";
            /*37*/    f<<m.Tymidine_incorporated()<<" , ";

            /*38*/    f<<m.Ag()<<" , ";
            /*39*/    f<<m.Ab()<<"\n";

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
            0,
            0,
            sp.TNF_deg_,sp.IFN_deg_,sp.Ag_deg_,tr.init_cells,tr.Ag,tr.Ab,sp.Prol_TymTr_,sp.max_num_cells_);


    APC=APC_cells(sp.init_ratio_APC_*tr.init_cells,
                  /*2*/ sp.IFN_APC0_prod_rate_,
                  /*3*/ sp.IFN_APCa_prod_rate_,
                  /*4*/ sp.IFN_APCbo_prod_rate_,
                  sp.IFN_APC_generic_prod_rate_,
                  /// 3) TNF Poductions rates of each type of APC
                  /*5*/ sp.TNF_APC0_prod_rate_,
                  /*6*/ sp.TNF_APCa_prod_rate_,
                  /*7*/ sp.APC_TNF_Induction_CD137,
                  sp.TNF_APC_generic_prod_rate_,
                  /// 4) Percentages of IFN productions of each type of APC
                  /*8*/ sp.percentage_IFN_APC0_prod_rate_,
                  /*9*/ sp.percentage_IFN_APCa_prod_rate_,
                  //*10*/ sp.percentage_IFN_APCbo_prod_rate_,
                  /// 5)Percentages of TNF productions of each type of APC
                  /*11*/ sp.percentage_TNF_APC0_prod_rate_,
                  /// 6) Proliferation rates
                  /*12*/ sp.APC_bound_proliferation_rate_,
                  /// 7) Apoptosis rates
                  /*13*/ sp.APC0_apop_rate_,
                  /*14*/ sp.APCa_apop_rate_,
                  /*15*/ sp.APCbo_apop_rate_,
                  sp.APC_generic_apop_rate_,
                  /// 8) constant saturation of TNF for apoptosis
                  /*16*/ sp.Ks_APC_m_TNF_,
                  /// 9) conversion rates
                  /*17*/ sp.APC_Ag_,
                  /*18*/ sp.APC_APC_,
                  /*19*/ sp.APC_NK_,
                  /*20*/ sp.APC_LT_1_,
                  /*21*/ sp.APC_Ag_2_,
                  /*22*/ sp.APC_Ab_,
//                  /// 10)Saturation constant of IFN and TNF for activation
//                  /*23*/ sp.KsAPC_LT_,
                  /// 11)Saturation constant of APC_LT interaction
                  /*24*/ sp.APC_Ksi_,
                  /*25*/ sp.APC_Kst_,
                  /// 12) Percentages of cell expressing receptor
                  /*26*/ sp.APC0_expressing_receptor_,
                  /// 13) Apoptosis rate for TNF
                  /*27*/ sp.u_APC_TNF_);

    NK=NK_cells (sp.init_ratio_NK_*tr.init_cells,
                 /// 2) IFN Poductions rates of each type of NK
                 /*2*/ sp.IFN_NK0_prod_rate_,
                 /*3*/ sp.IFN_NKa_prod_rate_,
                 /*4*/ sp.IFN_NKbo_prod_rate_,
                       sp.IFN_NK_generic_prod_rate_,
                 /// 3) TNF Poductions rates of each type of NK
                 /*5*/ sp.TNF_NK0_prod_rate_,
                 /*6*/ sp.TNF_NKa_prod_rate_,
                 /*7*/ sp.TNF_NKbo_prod_rate_,
                       sp.TNF_NK_generic_prod_rate_,
                 /// 4) Percentages of IFN productions of each type of NK
                 /*8*/ sp.percentage_IFN_NK0_prod_rate_,
                 /// 5)Percentages of TNF productions of each type of NK
                 /*9*/ sp.percentage_TNF_NK0_prod_rate_,
                 /*10*/ sp.percentage_TNF_NKa_prod_rate_,
                 //*11*/ sp.percentage_TNF_NKbo_prod_rate_,
                 /// 6) Proliferation rates
                 /*12*/ sp.NK0_proliferation_rate_,
                 /*13*/ sp.NKa_proliferation_rate_,
                 sp.NK_generic_proliferation_rate_,
                 //*14*/ sp.NKbo_proliferation_rate_,
                 /// 7) Apoptosis rates
                 /*15*/ sp.NK0_apop_rate_,
                 /*16*/ sp.NKa_apop_rate_,
                 //*17*/ sp.NKbo_apop_rate_,
                 sp.NK_generic_apop_rate_,
                 /// 8) constant saturation of TNF for apoptosis
                 /*18*/ sp.Ks_NK_m_TNF_,
                 /// 9) conversion rates
                 /*19*/ sp.KaNK_,
                 /*20*/ sp.NK_NK_,
                 /*21*/ sp.NK_Ab_,
                 /// 10)Saturation constant of NK interaction for activation
                 /*22*/ sp.KsAPC_NK_,
                 /// 11)Saturation constant of NK_LT interaction
                 /*23*/ sp.NK_Ksi_,
                 /*24*/ sp.NK_Kst_,
                 /// 12) Percentages of cell expressing receptor
                 /*25*/ sp.NK0_expressing_receptor_,
                 /*26*/ sp.NKa_expressing_receptor_,
                 /// 13) Apoptosis rate for TNF
                 /*27*/ sp.u_NK_TNF_);


    LT=LT_cells  (sp.ratio_init_LTns_*tr.init_cells,
                  sp.ratio_initLTspecific_*tr.init_cells,
                  /// 2) IFN Poductions rates of each type of LT
                     /*3*/ sp.IFN_LTns_prod_rate_,
                     /*4*/ sp.IFN_LTbo_prod_rate_,
                     /*5*/ sp.IFN_LTbl_prod_rate_,
                  sp.IFN_LT_generic_prod_rate_,
                 /// 3) TNF Poductions rates of each type of LT
                     /*6*/ sp.TNF_LTns_prod_rate_,
                     /*7*/ sp.TNF_LTbo_prod_rate_,
                     /*8*/ sp.TNF_LTbl_prod_rate_,
                  sp.TNF_LT_generic_prod_rate_,
                 /// 4) Percentages of IFN productions of each type of LT
                     /*9*/ sp.percentage_IFN_LTns_prod_rate_,
                     /*10*/ sp.percentage_IFN_LTbo_prod_rate_,
                     //*11*/ sp.percentage_IFN_LTbl_prod_rate_,
                 /// 5)Percentages of TNF productions of each type of LT
                     /*12*/ sp.percentage_TNF_LTns_prod_rate_,
                     /*13*/ sp.percentage_TNF_LTbo_prod_rate_,
                     //*14*/ sp.percentage_TNF_LTbl_prod_rate_,
                 /// 6) Proliferation rates
                     /*15*/ sp.LTns_proliferation_rate_,
                     /*16*/ sp.LTbo_proliferation_rate_,
                     /*17*/ sp.LTbl_proliferation_rate_,
                            sp.LT_generic_proliferation_rate_,
                 /// 7) Apoptosis rates
                     /*18*/ sp.LTns_apop_rate_,
                     /*19*/ sp.LTbo_apop_rate_,
                     /*20*/ sp.LTbl_apop_rate_,
                            sp.LT_generic_apop_rate_,
                 /// 8) constant saturation of TNF for apoptosis
                     /*21*/ sp.Ks_LT_m_TNF_,
                 /// 9) Percentages of cell expressing receptor
                     /*22*/ sp.LTns_expressing_receptor_,
                 /// 10) Apoptosis rate for TNF
                     /*23*/ sp.u_LT_TNF_,
                  /// 12) apoptosis related parameters
                      /*24*/ sp.t_apop_meas_,
                      /*25*/ sp.t_duration_apoptosis_,
                  /// Lt Ab binding rate
                      /*26*/ sp.LT_Ab_
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
    fitPar_(other.fitPar_),
    prior_(other.prior_),
    current_(other.current_)



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

    std::swap(one.prior_,other.prior_);
    std::swap(one.current_,other.current_);

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

    while (trun_d+eps<=results.Duration()&&(!( m.TNF()!=m.TNF())))
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
              num_cellss[inum_cells].setMeasurement(m.num_cells());
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


std::vector<double> Cell_simulator::yfit (const Parameters& param0)const
{
    Parameters par(prior_);
    par.applyParameters(param0);
    Cell_simulator tmp(*this);


    std::vector<double> f=tmp.Simulate(par,experiment_).getData();
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


Cell_simulator& Cell_simulator::applyParameters(const Parameters& p,
                                const Treatment& tr)
{

    m=Media(p,tr);


    APC=APC_cells(p,tr);

    NK=NK_cells  (p,tr);


    LT=LT_cells  (p,tr);

    time_step_d=tr.time_step_d;
    sim_duration_d=tr.sim_duration_d;
    trun_d=0;

return *this;
}

Cell_simulator::Cell_simulator(const Parameters& prior,
                               const Parameters& current,
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
    prior_(prior),
    current_(current),
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

    RungeKutta4  RK(this,getState());

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
              num_cellss[inum_cells].setMeasurement (m.num_cells());
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


//	std::cout<<"\ntime \t"<<trun_d;
//	std::cout<<"\n media\n"<<m<<"\nLT\n"<<APC<<"\nNK\n"<<NK<<"\nLT\n"<<LT;
//	std::cin>>ch;
          if (0)
          {
          std::vector<double> y=RK.next(time_step_d);
        setState(y);
}
          else
{
        APC.update(time_step_d,m,NK,LT);
        NK.update(time_step_d,m,APC,LT);
        LT.update(time_step_d,trun_d,m,APC,NK);
        m.update(time_step_d,trun_d, APC,NK,LT);
        trun_d+=time_step_d;

 }

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

void Cell_simulator::Optimize(const Parameters& priorPar,
                                             const Experiment& experiment,
                              const std::string &filename)
{


    BayesIteration b(this,priorPar,&experiment,filename);
    b.getPosterior();



}


 Parameters Cell_simulator::getStandardParameters(){
     Parameters sp;
     sp.setMode("standard");
     /// 1) Init ratio of LT, NK, APC cells
         /*1*/ sp.push_back_1S("init_K_ratio_LT",0.89/(1-0.89),0.899/(1-0.899));// Fórmula leucocitaria
         /*2*/ sp.push_back_1S("init_K_ratio_APC_NK",0.8/(1-0.8),0.9/(1-0.9)); // Fórmula leucocitaria
         /*3*/ sp.push_back_1S("init_K_ratio_NK_APC",0.8/(1-0.8),0.9/(1-0.9));
     /// APC
         /// 2) IFN Poductions rates of each type of APC
         /*2*/ sp.push_back_1S("IFN_APC0_prod_rate",1.0e-8,1.0e-5);
         /*3*/ sp.push_back_1S("IFN_APCa_prod_rate",2.0e-6,2e-3);
         /*4*/ sp.push_back_1S("IFN_APCbo_prod_rate",1e-6,1e-3);
     sp.push_back_1S("IFN_APC_generic_prod_rate",0.2,10.0);
         /// 3) TNF Poductions rates of each type of APC
         /*5*/ sp.push_back_1S("TNF_APC0_prod_rate",1.0e-8,1.0e-5);
         /*6*/ sp.push_back_1S("TNF_APCa_prod_rate",2.0e-4,2.0);//k
         /*7*/ sp.push_back_1S("TNF_APCbo_prod_rate",1.0e-4,1.0);
     sp.push_back_1S("TNF_APC_generic_prod_rate",0.2,10.0);
         /// 4) Percentages of IFN productions of each type of APC
         /*8*/ sp.push_back_1S("Kpercentage_IFN_APC0_prod_rate",0.01,0.06);//oj
         /*9*/ sp.push_back_1S("Kpercentage_IFN_APCa_prod_rate",0.08,0.5);
         //*10*/ sp.push_back_1S("Kpercentage_IFN_APCbo_prod_rate",0.14,0.5);
         /// 5)Percentages of TNF productions of each type of APC
         /*11*/ sp.push_back_1S("Kpercentage_TNF_APC0_prod_rate",0.02,0.08);
         /// 6) Proliferation rates
         /*12*/ sp.push_back_1S("APC_bound_proliferation_rate",0.00001,0.001);
         /// 7) Apoptosis rates
         /*13*/ sp.push_back_1S("APC0_apop_rate",3.0e-6,3.0e-4);//k
         /*14*/ sp.push_back_1S("APCa_apop_rate",0.00001,0.001);//k
         /*15*/ sp.push_back_1S("APCbo_apop_rate",0.00001,0.001);//k
     sp.push_back_1S("APC_generic_apop_rate",0.2,10.0);
         /// 8) constant saturation of TNF for apoptosis
         /*16*/ sp.push_back_1S("Ks_APC_m_TNF",0.002,20);//K
         /// 9) conversion rates
         /*17*/ sp.push_back_1S("APC_Ag",1.0e-6,1e-1);//ojímetro (promedio de k)
         /*18*/ sp.push_back_1S("APC_APC",1.0e-10,1.0e-1);//ver
         /*19*/ sp.push_back_1S("APC_NK",1.0e-10,1.0e-1);//ver
         /*20*/ sp.push_back_1S("APC_LT_1",4.0e-7,4.0e-1);//K, multiplicar por la posibilidad de encuentro? No diferencio x que no tengo el dato, dejo que el programa modifique
         /*21*/ sp.push_back_1S("APC_Ag_2",4.0e-7,4.0e-1);//K, multiplicar por la posibilidad de encuentro?
         /*22*/ sp.push_back_1S("APC_Ab",1.0e-6,1e-1);
//         /// 10)Saturation constant of IFN and TNF for activation
//         /*23*/ sp.push_back_1S("KsAPC_LT",10.0e1,10.0e5);//k
         /// 11)Saturation constant of APC_LT interaction
         /*24*/ sp.push_back_1S("APC_Ksi",1.0e-2,1e2);//Kirschner
         /*25*/ sp.push_back_1S("APC_Kst",1.0e2,1.0e6);//Kirschner (promedio)
         /// 12) Percentages of cell expressing receptor
         /*26*/ sp.push_back_1S("APC0_Kratio_expressing_receptor",001,0.007);
         /// 13) Apoptosis rate for TNF
         /*27*/ sp.push_back_1S("u_APC_TNF",0.0000417,0.417);//K

         /// NK
         /// 2) IFN Poductions rates of each type of NK
         /*2*/  sp.push_back_1S("IFN_NK0_prod_rate",1.0e-8,1.0e-5);
         /*3*/  sp.push_back_1S("IFN_NKa_prod_rate",2.0e-4,2.0);
         /*4*/  sp.push_back_1S("IFN_NKbo_prod_rate",1.0e-4,1.0);
      sp.push_back_1S("IFN_NK_generic_prod_rate",0.2,10.0);
         /// 3) TNF Poductions rates of each type of NK
         /*5*/  sp.push_back_1S("TNF_NK0_prod_rate",1.0e-8,1.0e-5);
         /*6*/  sp.push_back_1S("TNF_NKa_prod_rate",2.0e-6,2e-3);
         /*7*/  sp.push_back_1S("TNF_NKbo_prod_rate",1e-6,1e-3);
       sp.push_back_1S("IFN_NK_generic_prod_rate",0.2,10.0);
         /// 4) Percentages of IFN productions of each type of NK
         /*8*/  sp.push_back_1S("Kpercentage_IFN_NK0_prod_rate",0.1,0.06);
         /// 5)Percentages of TNF productions of each type of NK
         /*9*/  sp.push_back_1S("Kpercentage_TNF_NK0_prod_rate",0.1,0.06);
         /*10*/  sp.push_back_1S("Kpercentage_TNF_NKa_prod_rate",0.05,0.5);
         //*11*/  sp.push_back_1S("Kpercentage_TNF_NKbo_prod_rate",0.04,0.5);
         /// 6) Proliferation rates
         /*12*/  sp.push_back_1S("NK0_proliferation_rate",3.0e-6,3.0e-4);
         /*13*/  sp.push_back_1S("NKa_proliferation_rate",0.00001,0.001);
        sp.push_back_1S("NK_generic_proliferation_rate",0.2,10.0);
         //*14*/  sp.push_back_1S("NKbo_proliferation_rate",0.00001,0.001);
         /// 7) Apoptosis rates
         /*15*/  sp.push_back_1S("NK0_apop_rate",3.0e-6,3.0e-4);
         /*16*/  sp.push_back_1S("NKa_apop_rate",0.00001,0.001);
         //*17*/  sp.push_back_1S("NKbo_apop_rate",0.00001,0.001);
         sp.push_back_1S("NK_generic_apop_rate",0.2,10.0);
         /// 8) constant saturation of TNF for apoptosis
         /*18*/  sp.push_back_1S("Ks_NK_m_TNF",0.002,20);
         /// 9) conversion rates
         /*19*/  sp.push_back_1S("KaNK",1e-6,1e-2);
         /*20*/  sp.push_back_1S("NK_NK",1e-8,1e-1);
         /*21*/  sp.push_back_1S("NK_Ab",1e-8,1e-1);
         /// 10)Saturation constant of APC NK interaction for activation
         /*22*/  sp.push_back_1S("KsAPC_NK",0.005,50);
         /// 11)Saturation constant of NK_LT interaction
         /*23*/  sp.push_back_1S("NK_Ksi",1.0e-2,1e6);
         /*24*/  sp.push_back_1S("NK_Kst",1.0e-2,1e6);
         /// 12) Percentages of cell expressing receptor
         /*25*/  sp.push_back_1S("NK0_Kratio_expressing_receptor",0.0,0.03);
         /*26*/  sp.push_back_1S("NKa_Kratio_expressing_receptor",0.1,0.5);
         /// 13) Apoptosis rate for TNF
         /*27*/  sp.push_back_1S("u_NK_TNF",0.0000417,0.417);
         /// LT
         /// 1) Init number of LT
            /*2*/  sp.push_back_1S("Kratio_initLTspecific",0.0,0.05);//K
         /// 2) IFN Poductions rates of each type of LT
            /*3*/  sp.push_back_1S("IFN_LTns_prod_rate",0.0000002,0.002);
            /*4*/  sp.push_back_1S("IFN_LTbo_prod_rate",0.000002,0.02);//k
            /*5*/  sp.push_back_1S("IFN_LTbl_prod_rate",0.000001,0.01);
          sp.push_back_1S("IFN_LT_generic_prod_rate",0.2,10.0);
        /// 3) TNF Poductions rates of each type of LT
            /*6*/  sp.push_back_1S("TNF_LTns_prod_rate",0.00000002,0.0002);
            /*7*/  sp.push_back_1S("TNF_LTbo_prod_rate",0.0001,0.01);
            /*8*/  sp.push_back_1S("TNF_LTbl_prod_rate",0.00005,0.005);
           sp.push_back_1S("IFN_LT_generic_prod_rate",0.2,10.0);
        /// 4) Percentages of IFN productions of each type of LT
            /*9*/  sp.push_back_1S("Kpercentage_IFN_LTns_prod_rate",0.01,0.06);
            /*10*/  sp.push_back_1S("Kpercentage_IFN_LTbo_prod_rate",0.05,0.5);
            //*11*/  sp.push_back_1S("Kpercentage_IFN_LTbl_prod_rate",0.01,0.25);
        /// 5)Percentages of TNF productions of each type of LT
            /*12*/  sp.push_back_1S("Kpercentage_TNF_LTns_prod_rate",0.0,0.05);
            /*13*/  sp.push_back_1S("Kpercentage_TNF_LTbo_prod_rate",0.2,0.25);
            //*14*/  sp.push_back_1S("Kpercentage_TNF_LTbl_prod_rate",0.01,0.125);
        /// 6) Proliferation rates
            /*15*/  sp.push_back_1S("LTns_proliferation_rate",1.0/6000.0,1.0/60);//oj
            /*16*/  sp.push_back_1S("LTbo_proliferation_rate",0.083,0.83);//K
            /*17*/  sp.push_back_1S("LTbl_proliferation_rate",0.041,0.41);//e
            sp.push_back_1S("LT_generic_proliferation_rate",0.2,10.0);
        /// 7) Apoptosis rates
            /*18*/  sp.push_back_1S("LTns_apop_rate",0.0001,1);
            /*19*/  sp.push_back_1S("LTbo_apop_rate",0.055,0.55);
            /*20*/  sp.push_back_1S("LTbl_apop_rate",0.11,1.1);
             sp.push_back_1S("LT_generic_apop_rate",0.2,10.0);
       /// 8) constant saturation of TNF for apoptosis
            /*21*/  sp.push_back_1S("Ks_LT_m_TNF",0.0004,4.0);//k, promedio de LN y lung

        /// 9) Percentages of cell expressing receptor
            /*22*/  sp.push_back_1S("LTns_Kratio_expressing_receptor",0.01,0.1);

        /// 10) Apoptosis rate for TNF
            /*23*/  sp.push_back_1S("u_LT_TNF",1.0/240.0,10);//k


        /// 12) apoptosis related parameters
            /*24*/ sp.push_back_1S("t_duration_apoptosis",0.1,20);
            /*25*/ sp.push_back_1S("LT_Ab",0.01,10);

         /// Media
         /*1*/ sp.push_back_1S("TNF_deg",1.0/18,1.0/6);//k
         /*2*/ sp.push_back_1S("IFN_deg",1.0/18,1.0/6);//k
         /*3*/ sp.push_back_1S("Ag_deg",1.0/18,1/6);//oj
         /*4*/ sp.push_back_1S("Prol_TymTr",0.001,10);
     /*5*/ sp.push_back_1S("max_num_cells",1.0e6,2.0e6);
 }



 Parameters Cell_simulator::getMinimalParameters(){
     Parameters sp;
     sp.setMode("minimal");
     /// 1) Init ratio of LT, NK, APC cells
         /*1*/ sp.push_back_1S("init_K_ratio_LT",0.89/(1-0.89),0.899/(1-0.899));// Fórmula leucocitaria
         /*2*/ sp.push_back_1S("init_K_ratio_APC_NK",0.8/(1-0.8),0.9/(1-0.9)); // Fórmula leucocitaria
              /*3*/ sp.push_back_1S("init_K_ratio_NK_APC",0.30,0.6);
     /// APC
         /// 2) IFN Poductions rates of each type of APC
         /*2*/ sp.push_back_1S("IFN_APC0_prod_rate",1.0e-8,1.0e-5);
         /*3*/ sp.push_back_dB("IFN_APCa_prod_rate",2.0e-6,2e-3);
         /*4*/ sp.push_back_dB("IFN_APCbo_prod_rate",1e-6,1e-3);
         /// 3) TNF Poductions rates of each type of APC
         /*5*/ sp.push_back_dB("TNF_APC0_prod_rate",1.0e-8,1.0e-5);
         /*6*/ sp.push_back_dB("TNF_APCa_prod_rate",2.0e-4,2.0);//k
         /*7*/ sp.push_back_dB("TNF_APCbo_prod_rate",1.0e-4,1.0);
         /// 4) Percentages of IFN productions of each type of APC
         /*8*/ sp.push_back_dB("Kpercentage_IFN_APC0_prod_rate",0.01,0.06);//oj
         /*9*/ sp.push_back_dB("Kpercentage_IFN_APCa_prod_rate",0.08,0.5);
         //*10*/ sp.push_back_dB("Kpercentage_IFN_APCbo_prod_rate",0.14,0.5);
         /// 5)Percentages of TNF productions of each type of APC
         /*11*/ sp.push_back_dB("Kpercentage_TNF_APC0_prod_rate",0.02,0.08);
         /// 6) Proliferation rates
         /*12*/ sp.push_back_dB("APC_bound_proliferation_rate",0.00001,0.001);
         /// 7) Apoptosis rates
         /*13*/ sp.push_back_dB("APC0_apop_rate",3.0e-6,3.0e-4);//k
         /*14*/ sp.push_back_dB("APCa_apop_rate",0.00001,0.001);//k
         /*15*/ sp.push_back_dB("APCbo_apop_rate",0.00001,0.001);//k
     sp.push_back_1S("APC_generic_apop:rate", 0.1,2.0);
         /// 8) constant saturation of TNF for apoptosis
         /*16*/ sp.push_back_dB("Ks_APC_m_TNF",0.002,20);//K
         /// 9) conversion rates
         /*17*/ sp.push_back_dB("APC_Ag",1.0e-6,1e-1);//ojímetro (promedio de k)
         /*18*/ sp.push_back_dB("APC_APC",1.0e-10,1.0e-1);//ver
         /*19*/ sp.push_back_dB("APC_NK",1.0e-10,1.0e-1);//ver
         /*20*/ sp.push_back_dB("APC_LT_1",4.0e-7,4.0e-1);//K, multiplicar por la posibilidad de encuentro? No diferencio x que no tengo el dato, dejo que el programa modifique
         /*21*/ sp.push_back_dB("APC_Ag_2",4.0e-7,4.0e-1);//K, multiplicar por la posibilidad de encuentro?
         /*22*/ sp.push_back_dB("APC_Ab",1.0e-6,1e-1);
//         /// 10)Saturation constant of IFN and TNF for activation
//         /*23*/ sp.push_back_dB("KsAPC_LT",10.0e1,10.0e5);//k
         /// 11)Saturation constant of APC_LT interaction
         /*24*/ sp.push_back_dB("APC_Ksi",1.0e-2,1e2);//Kirschner
         /*25*/ sp.push_back_dB("APC_Kst",1.0e2,1.0e6);//Kirschner (promedio)
         /// 12) Percentages of cell expressing receptor
         /*26*/ sp.push_back_dB("APC0_Kratio_expressing_receptor",001,0.007);
         /// 13) Apoptosis rate for TNF
         /*27*/ sp.push_back_dB("u_APC_TNF",0.0000417,0.417);//K

         /// NK
         /// 2) IFN Poductions rates of each type of NK
         /*2*/  sp.push_back_dB("IFN_NK0_prod_rate",1.0e-8,1.0e-5);
         /*3*/  sp.push_back_dB("IFN_NKa_prod_rate",2.0e-4,2.0);
         /*4*/  sp.push_back_dB("IFN_NKbo_prod_rate",1.0e-4,1.0);
         /// 3) TNF Poductions rates of each type of NK
         /*5*/  sp.push_back_dB("TNF_NK0_prod_rate",1.0e-8,1.0e-5);
         /*6*/  sp.push_back_dB("TNF_NKa_prod_rate",2.0e-6,2e-3);
         /*7*/  sp.push_back_dB("TNF_NKbo_prod_rate",1e-6,1e-3);
         /// 4) Percentages of IFN productions of each type of NK
         /*8*/  sp.push_back_dB("Kpercentage_IFN_NK0_prod_rate",0.1,0.06);
         /// 5)Percentages of TNF productions of each type of NK
         /*9*/  sp.push_back_dB("Kpercentage_TNF_NK0_prod_rate",0.1,0.06);
         /*10*/  sp.push_back_dB("Kpercentage_TNF_NKa_prod_rate",0.05,0.5);
         //*11*/  sp.push_back_dB("Kpercentage_TNF_NKbo_prod_rate",0.04,0.5);
         /// 6) Proliferation rates
         /*12*/  sp.push_back_dB("NK0_proliferation_rate",3.0e-6,3.0e-4);
         /*13*/  sp.push_back_dB("NKa_proliferation_rate",0.00001,0.001);
         //*14*/  sp.push_back_dB("NKbo_proliferation_rate",0.00001,0.001);
         /// 7) Apoptosis rates
         /*15*/  sp.push_back_dB("NK0_apop_rate",3.0e-6,3.0e-4);
         /*16*/  sp.push_back_dB("NKa_apop_rate",0.00001,0.001);
         //*17*/  sp.push_back_dB("NKbo_apop_rate",0.00001,0.001);
         /// 8) constant saturation of TNF for apoptosis
         /*18*/  sp.push_back_dB("Ks_NK_m_TNF",0.002,20);
         /// 9) conversion rates
         /*19*/  sp.push_back_dB("KaNK",1e-6,1e-2);
         /*20*/  sp.push_back_dB("NK_NK",1e-8,1e-1);
         /*21*/  sp.push_back_dB("NK_Ab",1e-8,1e-1);
         /// 10)Saturation constant of APC NK interaction for activation
         /*22*/  sp.push_back_dB("KsAPC_NK",0.005,50);
         /// 11)Saturation constant of NK_LT interaction
         /*23*/  sp.push_back_dB("NK_Ksi",1.0e-2,1e6);
         /*24*/  sp.push_back_dB("NK_Kst",1.0e-2,1e6);
         /// 12) Percentages of cell expressing receptor
         /*25*/  sp.push_back_1S("NK0_Kratio_expressing_receptor",0.0,0.03);
         /*26*/  sp.push_back_1S("NKa_Kratio_expressing_receptor",0.1,0.5);
         /// 13) Apoptosis rate for TNF
         /*27*/  sp.push_back_dB("u_NK_TNF",0.0000417,0.417);
         /// LT
         /// 1) Init number of LT
            /*2*/  sp.push_back_1S("Kratio_initLTspecific",0.0,0.05);//K
         /// 2) IFN Poductions rates of each type of LT
            /*3*/  sp.push_back_dB("IFN_LTns_prod_rate",0.0000002,0.002);
            /*4*/  sp.push_back_dB("IFN_LTbo_prod_rate",0.000002,0.02);//k
            /*5*/  sp.push_back_dB("IFN_LTbl_prod_rate",0.000001,0.01);
        /// 3) TNF Poductions rates of each type of LT
            /*6*/  sp.push_back_dB("TNF_LTns_prod_rate",0.00000002,0.0002);
            /*7*/  sp.push_back_dB("TNF_LTbo_prod_rate",0.0001,0.01);
            /*8*/  sp.push_back_dB("TNF_LTbl_prod_rate",0.00005,0.005);
        /// 4) Percentages of IFN productions of each type of LT
            /*9*/  sp.push_back_dB("Kpercentage_IFN_LTns_prod_rate",0.01,0.06);
            /*10*/  sp.push_back_dB("Kpercentage_IFN_LTbo_prod_rate",0.05,0.5);
            //*11*/  sp.push_back_dB("Kpercentage_IFN_LTbl_prod_rate",0.01,0.25);
        /// 5)Percentages of TNF productions of each type of LT
            /*12*/  sp.push_back_dB("Kpercentage_TNF_LTns_prod_rate",0.0,0.05);
            /*13*/  sp.push_back_dB("Kpercentage_TNF_LTbo_prod_rate",0.2,0.25);
            //*14*/  sp.push_back_dB("Kpercentage_TNF_LTbl_prod_rate",0.01,0.125);
        /// 6) Proliferation rates
            /*15*/  sp.push_back_dB("LTns_proliferation_rate",1.0/6000.0,1.0/60);//oj
            /*16*/  sp.push_back_dB("LTbo_proliferation_rate",0.083,0.83);//K
            /*17*/  sp.push_back_dB("LTbl_proliferation_rate",0.041,0.41);//e
        /// 7) Apoptosis rates
            /*18*/  sp.push_back_dB("LTns_apop_rate",0.0001,1);
            /*19*/  sp.push_back_dB("LTbo_apop_rate",0.055,0.55);
            /*20*/  sp.push_back_dB("LTbl_apop_rate",0.11,1.1);
       /// 8) constant saturation of TNF for apoptosis
            /*21*/  sp.push_back_dB("Ks_LT_m_TNF",0.0004,4.0);//k, promedio de LN y lung

        /// 9) Percentages of cell expressing receptor
            /*22*/  sp.push_back_dB("LTns_Kratio_expressing_receptor",0.01,0.1);

        /// 10) Apoptosis rate for TNF
            /*23*/  sp.push_back_dB("u_LT_TNF",1.0/240.0,10);//k


        /// 12) apoptosis related parameters
            /*24*/ sp.push_back_dB("t_duration_apoptosis",0.1,20);
            /*25*/ sp.push_back_dB("LT_Ab",0.01,10);

         /// Media
         /*1*/ sp.push_back_dB("TNF_deg",1.0/18,1.0/6);//k
         /*2*/ sp.push_back_dB("IFN_deg",1.0/18,1.0/6);//k
         /*3*/ sp.push_back_dB("Ag_deg",1.0/18,1/6);//oj
         /*4*/ sp.push_back_dB("Prol_TymTr",0.001,10);
     /*5*/ sp.push_back_dB("max_num_cells",1.0e6,2.0e6);


   return sp;
 }



 std::ostream& Cell_simulator::put(std::ostream &s,const Parameters& param0) const
 {
     Parameters par(prior_);
     par.applyParameters(param0);
     Cell_simulator tmp(*this);


     Experiment f=tmp.Simulate(par,experiment_);

     s<<param0;

     s<<f;

     run(s,param0);


      return s;
 }


  std::vector<double> Cell_simulator::Derivative(double t, std::vector<double> y)
 {
    trun_d=t;
    setState(y);
    std::vector<double> DMedia=m.Derivative(APC,NK,LT);

    std::vector<double> DAPC=APC.Derivative(m,NK,LT);
    std::vector<double> DNK=NK.Derivative(m,APC);
    std::vector<double> DLT=LT.Derivative(trun_d,m,APC);

    DMedia.insert(DMedia.end(),DAPC.begin(),DAPC.end());
    DMedia.insert(DMedia.end(),DNK.begin(),DNK.end());
    DMedia.insert(DMedia.end(),DLT.begin(),DLT.end());

    return DMedia;



 }

 std::vector<double> Cell_simulator::getState()const
 {
     std::vector<double> stateMedia=m.getState();

     std::vector<double> stateAPC=APC.getState();
     std::vector<double> stateNK=NK.getState();
     std::vector<double> stateLT=LT.getState();

     stateMedia.insert(stateMedia.end(),stateAPC.begin(),stateAPC.end());
     stateMedia.insert(stateMedia.end(),stateNK.begin(),stateNK.end());
     stateMedia.insert(stateMedia.end(),stateLT.begin(),stateLT.end());

     return stateMedia;

 }

 void Cell_simulator::setState(const std::vector<double>& y)
 {
     std::size_t sizeMedia=2;
     std::size_t sizeAPC=7;
     std::size_t sizeNK=7;
     std::size_t sizeLT=7;
     std::size_t k=0;

     std::vector<double> yMedia(sizeMedia);
     for (std::size_t i=0; i<sizeMedia;i++)
     {
         yMedia[i]=y[k];
         k++;
     }

     std::vector<double> yAPC(sizeAPC);
     for (std::size_t i=0; i<sizeAPC;i++)
     {
         yAPC[i]=y[k];
         k++;
     }

     std::vector<double> yNK(sizeNK);
     for (std::size_t i=0; i<sizeNK;i++)
     {
         yNK[i]=y[k];
         k++;
     }

     std::vector<double> yLT(sizeLT);
     for (std::size_t i=0; i<sizeLT;i++)
     {
         yLT[i]=y[k];
         k++;
     }

     m.setState(yMedia,trun_d,APC,NK,LT);
     APC.setState(yAPC);
     NK.setState(yNK);
     LT.setState(yLT);
 }

