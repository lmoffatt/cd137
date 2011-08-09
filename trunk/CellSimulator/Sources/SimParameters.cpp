#include "Includes/SimParameters.h"


SimParameters::SimParameters():


    max_num_cells(2e6),
    init_ratio_APC_cells (1e5),
    init_ratio_NK_cells(1e5),
    init_ratio_LT_cells (9e5),
    LT_ratio_specific (1000),
    APC_max_proliferation_rate_ (1.0/240),
    NK_max_proliferation_rate_(1.0/240),
    LT_max_no_receptor_prol_rate_(1.0/240),
    LT_max_free_prol_rate_(1.0/2),
    LT_max_bound_prol_rate_(1.0/1.2),
    LT_max_blocked_prol_rate_ (1.0/2),
    APC_no_to_free_rate_per_Ag_(1.0/30),
    APC_free_to_bound_rate_per_LT_(1.0/1e5),
    APC_Ab_binding_rate_(1.0/30),
    APC_exh_rate(1.0/1e5),
    NK_no_to_free_rate_per_Ag_(1.0/30),
    NK_free_to_bound_rate_per_LT_(1.0/1e5),
    NK_Ab_binding_rate(9/1e5),
    NK_exh_rate(1.0/1e5),
    LT_no_to_free_rate_per_APC_(1.0/6e5),
    LT_free_to_bound_rate_per_APC_(1.0/1e5),
    LT_mAb_binding_rate_(1.0/1e5),
    APC_IFN_free_prod_rate_(0.5/1e5),
    APC_IFN_Ag_prod_rate_(5.0/1e5),
    APC_IFN_bound_prod_rate_(10.0/1e5),
    APC_IFN_blocked_prod_rate_(5/1e5),
    NK_IFN_free_prod_rate_(0.5/1e5),
    NK_IFN_Ag_prod_rate_(5.0/1e5),
    NK_IFN_bound_prod_rate_(10.0/1e5),
    NK_IFN_blocked_prod_rate_ (0.5/1e5),
    LT_IFN_no_rec_prod_rate_(0.001/1e5),
    LT_IFN_free_prod_rate_(101.0/1e5),
    LT_IFN_bound_prod_rate_(200.0/1e5),
    LT_IFN_blocked_prod_rate_(101.0/1e5),
    APC_TNF_free_prod_rate_(5/1e5),
    APC_TNF_Ag_prod_rate_(570/1e5),
    APC_TNF_bound_prod_rate_(1110/1e5),
    APC_TNF_blocked_prod_rate_(570/1e5),
    NK_TNF_free_prod_rate_(5/1e5),
    NK_TNF_Ag_prod_rate_(570/1e5),
    NK_TNF_bound_prod_rate_(1110/1e5),
    NK_TNF_blocked_prod_rate_ (5/1e5),
    LT_TNF_no_rec_prod_rate_(0.001/1e5),
    LT_TNF_free_prod_rate_(10.0/1e5),
    LT_TNF_bound_prod_rate_(20.0/1e5),
    LT_TNF_blocked_prod_rate_(10.0/1e5),




    // Ag_internalization_rate (0.5),
    TNF_deg (0.5),
    IFN_deg (0.5)
    {};
