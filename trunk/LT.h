#ifndef LT_H_INCLUDED
#define LT_H_INCLUDED
class Media;
class APC_cells;
class NK_cells;

/// Lymphocytes T cells
class LT_cells
{
    public:
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

    /// Number of cell that express the receptor
        double num_cells_expressing_receptor()const;

    /// Number of cells that express the receptor and it is bound
        double num_cells_expressing_receptor_and_bound ()const;

        void update(double time_step,const Media& m, const APC_cells& a, const NK_cells& NK);

        LT_cells(double num_LT_init_,
                 double LT_num_specific_,
                 double LT_max_no_receptor_prol_rate_,
                 double LT_max_free_prol_rate_,
                 double LT_max_bound_prol_rate_,
                 double IFN_no_rec_prod_rate_,
                 double IFN_free_prod_rate_,
                 double IFN_bound_prod_rate_,
                 double TNF_no_rec_prod_rate_,
                 double TNF_free_prod_rate_,
                 double TNF_bound_prod_rate_,
                 double LT_no_to_free_rate_per_APC_,
                 double LT_free_to_bound_rate_per_APC_):
            num_non_Agsp_d(num_LT_init_),
            num_Agsp_bound_receptor_d(0),
            num_Agsp_free_receptor_d(0),
            num_Agsp_no_receptor_d(LT_num_specific_),
            IFN_free_prod_rate_d(IFN_free_prod_rate_),
            IFN_no_rec_prod_rate_d(IFN_no_rec_prod_rate_),
            IFN_bound_prod_rate_d(IFN_bound_prod_rate_),
            TNF_free_prod_rate_d(TNF_free_prod_rate_),
            TNF_no_rec_prod_rate_d(TNF_no_rec_prod_rate_),
            TNF_bound_prod_rate_d(TNF_bound_prod_rate_),
            LT_max_no_receptor_prol_rate_d(LT_max_no_receptor_prol_rate_),
            LT_max_free_prol_rate_d(LT_max_free_prol_rate_),
            LT_max_bound_prol_rate_d(LT_max_bound_prol_rate_),
            LT_no_to_free_rate_per_APC_d(LT_no_to_free_rate_per_APC_),
            LT_free_to_bound_rate_per_APC_d (LT_free_to_bound_rate_per_APC_) {};


            LT_cells(){};

    private:

        /// number of non Ag specific cells
        double num_non_Agsp_d;

        /// number of Ag specific cells that have no receptor
        double num_Agsp_no_receptor_d;

        /// number of Ag specific cells that have the receptor and it is free
        double num_Agsp_free_receptor_d;

        /// number of Ag specific cells that have the receptor but bound to its ligand
        double num_Agsp_bound_receptor_d;

        /// those are parameters that do not vary

        ///IFN production rates
        double IFN_no_rec_prod_rate_d;
        double IFN_free_prod_rate_d;
        double IFN_bound_prod_rate_d;

        ///TNF production rates
        double TNF_no_rec_prod_rate_d;
        double TNF_free_prod_rate_d;
        double TNF_bound_prod_rate_d;

        ///proliferation rates
        double LT_max_no_receptor_prol_rate_d;
        double LT_max_free_prol_rate_d;
        double LT_max_bound_prol_rate_d;


        ///conversion rates
        double LT_no_to_free_rate_per_APC_d;
        double LT_free_to_bound_rate_per_APC_d;

};

#endif // LT_H_INCLUDED
