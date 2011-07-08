#ifndef NK_H_INCLUDED
#define NK_H_INCLUDED

class Media;
class APC_cells;
class LT_cells;

/// NK cells
class NK_cells
{
    public:
        /// Total NK cells
        double num() const;

        /// Total production of interpheron gamma
        double IFNgamma_production_rate() const;

        /// Total production of Tumor Necrosis Factor
        double TNF_production_rate() const;

        /// Number of cells that have bound the antigen and express the receptor
        double& NK_num_Ag();
        const double& NK_num_Ag()const;

        /// Number of cells that bound to the blocking mAb;
        double& NK_blocked();
        const double& NK_blocked()const;

        /// Number of cells that have express the ligand/receptor and bound
        double& NK_num_bound();
        const double& NK_num_bound() const;

        /// Percentage of cells exprssing receptor
        double percentage_cell_expressing_receptor () const;

        /// Number of cells that are exhausted
        double& NK_exhausted ();
        const double& NK_exhausted()const;

        /// Ag internalization rate
        const double& NK_no_to_free_rate_per_Ag()const;
        double& NK_no_to_free_rate_per_Ag();

        /// Receptor binding rate
        const double& NK_free_to_bound_rate_per_LT()const;
        double& NK_free_to_bound_rate_per_LT();

        /// Ab binding rate
        const double& NK_Ab_binding_rate () const;
        double& NK_Ab_binding_rate();

        /// Exhaustation rate
        const double NK_exh_rate () const;
        double& NK_exh_rate();


        void update(double time_step,const Media& m, const APC_cells& APC,const LT_cells& LT);

        NK_cells (double NK_init,
                  double NK_max_proliferation_rate_,
                  double NK_no_to_free_rate_per_Ag_ ,
                  double NK_free_to_bound_rate_per_LT,
                  double NK_Ab_binding_rate,
                  double NK_exh_rate,
                  double NK_IFN_free_prod_rate_,
                  double NK_IFN_Ag_prod_rate_,
                  double NK_IFN_bound_prod_rate_,
                  double NK_IFN_blocked_prod_rate_,
                  double NK_TNF_free_prod_rate_,
                  double NK_TNF_Ag_prod_rate_,
                  double NK_TNF_bound_prod_rate_,
                  double NK_TNF_blocked_prod_rate_
                  ):

                  NK_num_free_d(NK_init),
                  NK_num_Ag_d(0),
                  NK_num_LT_bound_d(0),
                  NK_blocked_d(0),
                  NK_num_exhausted_d(0),
                  NK_exh_rate_d (NK_exh_rate),
                  NK_IFN_free_prod_rate_d(NK_IFN_free_prod_rate_),
                  NK_IFN_Ag_prod_rate_d(NK_IFN_Ag_prod_rate_),
                  NK_IFN_bound_prod_rate_d(NK_IFN_bound_prod_rate_),
                  NK_IFN_blocked_prod_rate_d(NK_IFN_blocked_prod_rate_),
                  NK_TNF_free_prod_rate_d(NK_TNF_free_prod_rate_),
                  NK_TNF_Ag_prod_rate_d(NK_TNF_Ag_prod_rate_),
                  NK_TNF_bound_prod_rate_d(NK_TNF_bound_prod_rate_),
                  NK_TNF_blocked_prod_rate_d(NK_TNF_blocked_prod_rate_),
                  NK_max_proliferation_rate_d(NK_max_proliferation_rate_),
                  NK_no_to_free_rate_per_Ag_d(NK_no_to_free_rate_per_Ag_),
                  NK_free_to_bound_rate_per_LT_d (NK_free_to_bound_rate_per_LT),
                  NK_Ab_binding_rate_d (NK_Ab_binding_rate)

                  {};

        NK_cells(){};

        /// main step for the NK cells


    private:


        /// number of free cells
        double NK_num_free_d;
        /// number of cells that have internalized  the antigen (and therefore express the ligand)
        double NK_num_Ag_d;
        /// number of cells that are blocked
        double NK_blocked_d;
        /// number of cells that have the receptor bound to its ligand
        double NK_num_LT_bound_d;
        /// number of cells that are exhausted
        double NK_num_exhausted_d;

        ///TNF and INF Poductions rates of each type of APC
        double NK_IFN_free_prod_rate_d;
        double NK_IFN_Ag_prod_rate_d;
        double NK_IFN_bound_prod_rate_d;
        double NK_IFN_blocked_prod_rate_d;
        double NK_TNF_free_prod_rate_d;
        double NK_TNF_Ag_prod_rate_d;
        double NK_TNF_bound_prod_rate_d;
        double NK_TNF_blocked_prod_rate_d;

        /// those are parameters that do not vary
        double NK_max_proliferation_rate_d;
        double NK_no_to_free_rate_per_Ag_d;
        double NK_free_to_bound_rate_per_LT_d;
        double NK_Ab_binding_rate_d;
        double NK_exh_rate_d;

};
#endif // NK_H_INCLUDED
