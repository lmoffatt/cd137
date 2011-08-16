#ifndef NK_H_INCLUDED
#define NK_H_INCLUDED
#include "Includes/SimParameters.h"
#include "Includes/Treatment.h"
#include <iostream>

class Media;
class APC_cells;
class LT_cells;

/// NK cells
class NK_cells
{
    public:
    ~NK_cells(){}
        /// Total NK cells
        double num() const;

        /// Total production of interpheron gamma
        double IFNgamma_production_rate() const;

        /// Total production of Tumor Necrosis Factor
        double TNF_production_rate() const;

        /// Number of cells free
        double& NK_num_free();
        const double& NK_num_free()const;

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
                  double NK_IFN_free_prod_rate_,
                  double NK_IFN_Ag_prod_rate_,
                  double NK_IFN_bound_prod_rate_,
                  double NK_IFN_blocked_prod_rate_,
                  double NK_TNF_free_prod_rate_,
                  double NK_TNF_Ag_prod_rate_,
                  double NK_TNF_bound_prod_rate_,
                  double NK_TNF_blocked_prod_rate_,
                  double NK_exh_rate
                  );

        NK_cells (const SimParameters& sp,
                  const Treatment& tr);

        NK_cells& operator=(const NK_cells& other);

        friend void swap(NK_cells& one, NK_cells& other);

        NK_cells(const NK_cells& other);

        NK_cells();

        void reset(const SimParameters& sp,const Treatment& tr);

        /// main step for the NK cells


	friend std::ostream& operator<<(std::ostream& s, const NK_cells& c);


    private:


        /// number of free cells
        double NK_num_free_d;
        /// number of cells that have internalized  the antigen (and therefore express the ligand)
        double NK_num_Ag_d;
        /// number of cells that have the receptor bound to its ligand
        double NK_num_LT_bound_d;
        /// number of cells that are blocked
        double NK_blocked_d;
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
