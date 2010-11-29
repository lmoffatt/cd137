#ifndef CELL_SIMULATOR_H_INCLUDED
#define CELL_SIMULATOR_H_INCLUDED
#include <iostream>
#include <string>

class APC_cells;
class LT_cells;


/// Composition of the media that sorrounds the cells.
/// This version does not have volume, all variables expressed as concentrations
class Media
{
public:
    /// Interpheron gamma concentration in the media
    double& IFNgamma()
    {
        return IFNgamma_d;
    };
    const double& IFNgamma()const
    {
        return IFNgamma_d;
    };


    /// maximum number of cells tolerated by the media
    /// proliferation occurs when the number of cells is lower than this number
    double& Max_num_cells()
    {
        return max_num_cells_d;
    };
    const double& Max_num_cells() const
    {
        return max_num_cells_d;
    };

    /// number of cells contained by the media
    double& num_cells()
    {
        return num_cells_d;
    };
    const double& num_cells() const
    {
        return num_cells_d;
    };

    /// Antigen concentration
    double& AG()
    {
        return AG_d;
    };
    const double& AG()const
    {
        return AG_d;
    };


    void update(double time_step,const APC_cells& APC_,const LT_cells& LT_);

    Media(double& max_num_cells_,double init_num_cells_ ,double AG_=0,double IFNgamma_init=0):
        IFNgamma_d(IFNgamma_init),
        max_num_cells_d(max_num_cells_),
        num_cells_d(init_num_cells_),
        AG_d(AG_) {}

    Media(){}



private:
    double IFNgamma_d;
    double max_num_cells_d;
    double num_cells_d;
    double AG_d;
    double Ab_d;

};

class LT_cells;


/// Antigen presentation cells
class APC_cells
{
public:
    /// number of cells
    double num() const
    {
        return num_AG_d+num_free_d+num_LT_bound_d;
    };

    /// total production of interpheron gamma
    double IFNgamma_production_rate() const{

        double sum=IFN_free_prod_rate_d*num_free_d+
        IFN_AG_prod_rate_d*num_AG_d+
        IFN_bound_prod_rate_d*num_LT_bound_d;
        return sum;
    };

    /// number of cells that have bound the antigen
    double& num_AG()
    {
        return num_AG_d;
    };
    const double& num_AG()const
    {
        return num_AG_d;
    };

    double num_bound() const{
    return num_LT_bound_d;
    };




    void update(double time_step,const Media& m, const LT_cells& LT);


    APC_cells(double APC_init,
              double max_proliferation_rate_,
              double no_to_free_rate_per_AG_ ,
              double free_to_bound_rate_per_LT,
              double IFN_free_prod_rate_,
              double IFN_AG_prod_rate_,
              double IFN_bound_prod_rate_):

        num_free_d(APC_init),
        num_AG_d(0),
        num_LT_bound_d(0),
        IFN_free_prod_rate_d(IFN_free_prod_rate_),
        IFN_AG_prod_rate_d(IFN_AG_prod_rate_),
        IFN_bound_prod_rate_d(IFN_bound_prod_rate_),
        max_proliferation_rate_d(max_proliferation_rate_),
        no_to_free_rate_per_AG_d(no_to_free_rate_per_AG_),
        free_to_bound_rate_per_LT_d (free_to_bound_rate_per_LT) {};


APC_cells(){};

private:
/// this are variables that define the state of these cells

    /// number of free cells
    double num_free_d;

    /// number of cells that have internalized  the antigen (and therefore express the ligand)
    double num_AG_d;

    /// number of cells that have the receptor bound to its ligand
    double num_LT_bound_d;

    /// number of cells that have the receptor bound to the blocking antibody
    double num_Ab_d;

    double IFN_free_prod_rate_d;
    double IFN_AG_prod_rate_d;
    double IFN_bound_prod_rate_d;




/// those are parameters that do not vary

    double max_proliferation_rate_d;
    double no_to_free_rate_per_AG_d;
    double free_to_bound_rate_per_LT_d;

};




/// Lymphocytes T cells
class LT_cells
{
public:
/// number of cells
    double num() const
    {
        return num_AGsp_bound_receptor_d+num_AGsp_free_receptor_d+num_AGsp_no_receptor_d+num_non_AGsp_d;
    };

    /// total production of interpheron gamma
    double IFNgamma_production_rate() const
    {
        return (num_non_AGsp_d+num_AGsp_no_receptor_d)*IFN_no_rec_prod_rate_d+
               num_AGsp_free_receptor_d*IFN_free_prod_rate_d+
               num_AGsp_bound_receptor_d*IFN_bound_prod_rate_d;
    };

///number of non specific cells
   double num_cells_not_AG_specific()const
    {
        return num_non_AGsp_d;
    };



   ///number of cells without receptor
   double num_cells_not_expressing_receptor()const
    {
        return num_AGsp_no_receptor_d;
    };


/// number of cells that express the receptor and it is free
    double num_cells_expressing_receptor_and_free()const
    {
        return num_AGsp_free_receptor_d;
    };


    double num_cells_expressing_receptor()const
    {
        return num_AGsp_free_receptor_d+num_AGsp_bound_receptor_d;
    };

/// number of cells that express the receptor and it is bound
    double num_cells_expressing_receptor_and_bound()const
    {
        return num_AGsp_bound_receptor_d;
    };



    void update(double time_step,const Media& m, const APC_cells& a);

    LT_cells(double num_LT_init_,
             double num_specific_,
             double max_no_receptor_prol_rate_,
             double max_free_prol_rate_,
             double max_bound_prol_rate_,
             double IFN_no_rec_prod_rate_,
             double IFN_free_prod_rate_,
             double IFN_bound_prod_rate_,
             double no_to_free_rate_per_APC_,
             double free_to_bound_rate_per_APC_):
        num_non_AGsp_d(num_LT_init_),
        num_AGsp_bound_receptor_d(0),
        num_AGsp_free_receptor_d(0),
        num_AGsp_no_receptor_d(num_specific_),
        IFN_free_prod_rate_d(IFN_free_prod_rate_),
        IFN_no_rec_prod_rate_d(IFN_no_rec_prod_rate_),
        IFN_bound_prod_rate_d(IFN_bound_prod_rate_),
        max_no_receptor_prol_rate_d(max_no_receptor_prol_rate_),
        max_free_prol_rate_d(max_free_prol_rate_),
        max_bound_prol_rate_d(max_bound_prol_rate_),
        no_to_free_rate_per_APC_d(no_to_free_rate_per_APC_),
        free_to_bound_rate_per_APC_d (free_to_bound_rate_per_APC_) {};


        LT_cells(){};

private:
/// this are variables that define the state of these cells

    /// number of non Ag specific cells
    double num_non_AGsp_d;

    /// number of Ag specific cells that have no receptor
    double num_AGsp_no_receptor_d;

    /// number of Ag specific cells that have the receptor and it is free
    double num_AGsp_free_receptor_d;

    /// number of Ag specific cells that have the receptor but bound to its ligand
    double num_AGsp_bound_receptor_d;





/// those are parameters that do not vary

///IFN production rates
    double IFN_no_rec_prod_rate_d;
    double IFN_free_prod_rate_d;
    double IFN_bound_prod_rate_d;

///proliferation rates
double max_no_receptor_prol_rate_d;
  double max_free_prol_rate_d;
  double max_bound_prol_rate_d;


///conversion rates
  double no_to_free_rate_per_APC_d;
  double free_to_bound_rate_per_APC_d;

};


class Cell_simulator
{
public:

    void ask_parameters();
    void run();
    void show_results()const;
    void update(double time_step);

    Cell_simulator(){};

private:
    LT_cells LT;
    APC_cells APC;
    Media   m;

    /// time in hours
    double time_step_d;
    double trun_d;

    double sim_duration_d;

    std::string filename;

};




#endif // CELL_SIMULATOR_H_INCLUDE
