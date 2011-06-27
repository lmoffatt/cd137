#ifndef APC_H_INCLUDED
#define APC_H_INCLUDED

class Media;
class LT_cells;
class NK_cells;


/// Antigen presentation cells (Monocytes)
class APC_cells
{
    public:
        /// Total APC cells
        double num() const;

        /// Total production of interpheron gamma
        double IFNgamma_production_rate() const;

        /// Total production of Tumor Necrosis Factor
        double TNF_production_rate() const;

        /// Number of cells that have bound the antigen and express the receptor
        double& num_Ag();
        const double& num_Ag()const;

        /// Number of cells that have express the ligand and bound to LT
        double num_bound() const;

        /// Number of cells that are exhausted
        double& num_exhausted ();
        const double& num_exhausted()const;

        /// Ag internalization rate

        const double& no_to_free_rate_per_Ag()const;
        double& no_to_free_rate_per_Ag();


        void update(double time_step,const Media& m, const NK_cells& NK,const LT_cells& LT);

        APC_cells(double APC_init,
                  double APC_max_proliferation_rate_,
                  double APC_no_to_free_rate_per_Ag_ ,
                  double free_to_bound_rate_per_LT,
                  double APC_exh_rate,
                  double IFN_free_prod_rate_,
                  double IFN_Ag_prod_rate_,
                  double IFN_bound_prod_rate_,
                  double TNF_free_prod_rate_,
                  double TNF_Ag_prod_rate_,
                  double TNF_bound_prod_rate_
                  ):

                  num_free_d(APC_init),
                  num_Ag_d(0),
                  num_LT_bound_d(0),
                  num_exhausted_d(0),
                  IFN_free_prod_rate_d(IFN_free_prod_rate_),
                  IFN_Ag_prod_rate_d(IFN_Ag_prod_rate_),
                  IFN_bound_prod_rate_d(IFN_bound_prod_rate_),
                  TNF_free_prod_rate_d(TNF_free_prod_rate_),
                  TNF_Ag_prod_rate_d(TNF_Ag_prod_rate_),
                  TNF_bound_prod_rate_d(TNF_bound_prod_rate_),
                  APC_max_proliferation_rate_d(APC_max_proliferation_rate_),
                  APC_no_to_free_rate_per_Ag_d(APC_no_to_free_rate_per_Ag_),
                  APC_free_to_bound_rate_per_LT_d (free_to_bound_rate_per_LT),
                  APC_exh_rate_d (APC_exh_rate)
                  {};

        APC_cells(){};

        /// main step for the APC cells


    private:


        /// number of free cells
        double num_free_d;
        /// number of cells that have internalized  the antigen (and therefore express the ligand and receptor)
        double num_Ag_d;
        /// number of cells that have the receptor bound to its ligand
        double num_LT_bound_d;
        /// number of cells that are exhausted
        double num_exhausted_d;

        ///TNF and INF Poductions rates of each type of APC
        double IFN_free_prod_rate_d;
        double IFN_Ag_prod_rate_d;
        double IFN_bound_prod_rate_d;
        double TNF_free_prod_rate_d;
        double TNF_Ag_prod_rate_d;
        double TNF_bound_prod_rate_d;

        /// those are parameters that do not vary
        double APC_max_proliferation_rate_d;
        double APC_no_to_free_rate_per_Ag_d;
        double APC_free_to_bound_rate_per_LT_d;
        double APC_exh_rate_d;

};



#endif // APC_H_INCLUDED
