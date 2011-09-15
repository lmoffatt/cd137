#ifndef LT_H_INCLUDED
#define LT_H_INCLUDED
#include "Includes/SimParameters.h"
#include "Includes/Treatment.h"

class Media;
class APC_cells;
class NK_cells;

/// Lymphocytes T cells
class LT_cells
{
    public:

    ~LT_cells(){}
    /// Number of cells
        double num() const;

    /// Total production of interpheron gamma
        double IFNgamma_production_rate() const;

    /// Total production of Tumor Necrosis alpha
        double TNF_production_rate() const;

    /// Number of non specific cells
        double num_cells_not_Ag_specific()const;

    /// Number of cells without receptor
        double num_cells_not_expressing_receptor()const;

    /// Number of cells that express the receptor and it is free
        double num_cells_expressing_receptor_and_free()const;

    /// Number of cells that express the receptor and it is blocked
        double num_blocked () const;

    /// Number of LT exhausted cells
        double num_exhausted () const;

    /// Number of cell that express the receptor
        double num_cells_expressing_receptor()const;

    /// Percentage of cell expressing the receptor
        double LT_percentage_cell_expressing_receptor () const;

    /// Number of cells that express the receptor and it is bound
        double num_cells_expressing_receptor_and_bound ()const;


        void update(double time_step,const Media& m, const APC_cells& a, const NK_cells& NK);

        LT_cells(double num_LT_init_,
                 double LT_num_specific_,
                 double LT_max_no_receptor_prol_rate_,
                 double LT_max_free_prol_rate_,
                 double LT_max_bound_prol_rate_,
                 double LT_max_blocked_prol_rate_,
                 double IFN_no_rec_prod_rate_,
                 double IFN_free_prod_rate_,
                 double IFN_bound_prod_rate_,
                 double IFN_blocked_prod_rate_,
                 double TNF_no_rec_prod_rate_,
                 double TNF_free_prod_rate_,
                 double TNF_bound_prod_rate_,
                 double TNF_blocked_prod_rate_,
                 double LT_no_to_free_rate_per_APC_,
                 double LT_free_to_bound_rate_per_APC_,
                 double LT_mAb_binding_rate_,
                 double LT_exh_rate_);


        LT_cells(const SimParameters& sp,
                 const Treatment& tr):
            num_non_Agsp_d(sp.init_ratio_LT_cells_*tr.init_cells),
            num_Agsp_no_receptor_d(sp.LT_ratio_specific_*tr.init_cells),
            num_Agsp_free_receptor_d(0),
            num_Agsp_bound_receptor_d(0),
            num_blocked_d(0),
            IFN_no_rec_prod_rate_d(sp.LT_IFN_no_rec_prod_rate_),
            IFN_free_prod_rate_d(sp.LT_IFN_free_prod_rate_),
            IFN_bound_prod_rate_d(sp.LT_IFN_bound_prod_rate_),
            IFN_blocked_prod_rate_d(sp.LT_IFN_blocked_prod_rate_),
            TNF_no_rec_prod_rate_d(sp.LT_TNF_no_rec_prod_rate_),
            TNF_free_prod_rate_d(sp.LT_TNF_free_prod_rate_),
            TNF_bound_prod_rate_d(sp.LT_TNF_bound_prod_rate_),
            TNF_blocked_prod_rate_d (sp.LT_TNF_blocked_prod_rate_),
            LT_max_no_receptor_prol_rate_d(sp.LT_max_no_receptor_prol_rate_),
            LT_max_free_prol_rate_d(sp.LT_max_free_prol_rate_),
            LT_max_bound_prol_rate_d(sp.LT_max_bound_prol_rate_),
            LT_max_blocked_prol_rate_d(sp.LT_max_blocked_prol_rate_),
            LT_no_to_free_rate_per_APC_d(sp.LT_no_to_free_rate_per_APC_),
            LT_free_to_bound_rate_per_APC_d (sp.LT_free_to_bound_rate_per_APC_),
            LT_mAb_binding_rate_d (sp.LT_mAb_binding_rate_),
            exh_rate_d (sp.LT_exh_rate_){}

            LT_cells(){}

            LT_cells(const LT_cells& other);

            LT_cells& operator=(const LT_cells& other);

            friend void swap(LT_cells& one, LT_cells& other);
	    friend std::ostream& operator<<(std::ostream& s, const LT_cells& c);


            void reset(const SimParameters& sp,const Treatment& tr);

    private:

        /// number of non Ag specific cells
        double num_non_Agsp_d;

        /// number of Ag specific cells that have no receptor
        double num_Agsp_no_receptor_d;

        /// number of Ag specific cells that have the receptor and it is free
        double num_Agsp_free_receptor_d;

        /// number of Ag specific cells that have the receptor but bound to its ligand
        double num_Agsp_bound_receptor_d;

        /// number of cells that have bound to the mAb
        double num_blocked_d;

        /// number of LT exhausted
        double num_exhausted_d;

        /// those are parameters that do not vary

        ///IFN production rates
        double IFN_no_rec_prod_rate_d;
        double IFN_free_prod_rate_d;
        double IFN_bound_prod_rate_d;
        double IFN_blocked_prod_rate_d;

        ///TNF production rates
        double TNF_no_rec_prod_rate_d;
        double TNF_free_prod_rate_d;
        double TNF_bound_prod_rate_d;
        double TNF_blocked_prod_rate_d;

        ///proliferation rates
        double LT_max_no_receptor_prol_rate_d;
        double LT_max_free_prol_rate_d;
        double LT_max_bound_prol_rate_d;
        double LT_max_blocked_prol_rate_d;


        ///conversion rates
        double LT_no_to_free_rate_per_APC_d;
        double LT_free_to_bound_rate_per_APC_d;
        double LT_mAb_binding_rate_d;
        double exh_rate_d;

};

#endif // LT_H_INCLUDED
