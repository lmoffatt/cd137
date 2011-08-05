
#ifndef SIMPARAMETERS_H_INCLUDED
#define SIMPARAMETERS_H_INCLUDED

struct SimParameters{
   double max_num_cells_;
   double max_num_cells;
   double init_num_APC_cells;
   double init_num_NK_cells;
   double init_num_LT_cells;
   double LT_num_specific;

   double sim_duration_d;
   double time_step_d;

   double Ag;
   double Ab;
  // double Ag_internalization_rate;

   double APC_max_proliferation_rate_;
   double NK_max_proliferation_rate_;
   double LT_max_no_receptor_prol_rate_;
   double LT_max_free_prol_rate_;
   double LT_max_bound_prol_rate_;
   double LT_max_blocked_prol_rate_;

   double APC_no_to_free_rate_per_Ag_;
   double APC_free_to_bound_rate_per_LT_;
   double APC_Ab_binding_rate_;
   double APC_exh_rate;
   double NK_no_to_free_rate_per_Ag_;
   double NK_free_to_bound_rate_per_LT_;
   double NK_Ab_binding_rate;
   double NK_exh_rate;
   double LT_no_to_free_rate_per_APC_;
   double LT_free_to_bound_rate_per_APC_;
   double LT_mAb_binding_rate_;

   double APC_IFN_free_prod_rate_;
   double APC_IFN_Ag_prod_rate_;
   double APC_IFN_bound_prod_rate_;
   double APC_IFN_blocked_prod_rate_;
   double NK_IFN_free_prod_rate_;
   double NK_IFN_Ag_prod_rate_;
   double NK_IFN_bound_prod_rate_;
   double NK_IFN_blocked_prod_rate_;
   double LT_IFN_no_rec_prod_rate_;
   double LT_IFN_free_prod_rate_;
   double LT_IFN_bound_prod_rate_;
   double LT_IFN_blocked_prod_rate_;

   double APC_TNF_free_prod_rate_;
   double APC_TNF_Ag_prod_rate_;
   double APC_TNF_bound_prod_rate_;
   double APC_TNF_blocked_prod_rate_;
   double NK_TNF_free_prod_rate_;
   double NK_TNF_Ag_prod_rate_;
   double NK_TNF_bound_prod_rate_;
   double NK_TNF_blocked_prod_rate_;
   double LT_TNF_no_rec_prod_rate_;
   double LT_TNF_free_prod_rate_;
   double LT_TNF_bound_prod_rate_;
   double LT_TNF_blocked_prod_rate_;



   double TNF_deg;
   double IFN_deg;



SimParameters();
};

#endif // SIMPARAMETERS_H_INCLUDED
