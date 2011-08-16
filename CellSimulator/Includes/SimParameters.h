
#ifndef SIMPARAMETERS_H_INCLUDED
#define SIMPARAMETERS_H_INCLUDED

#include <vector>
#include<string>
#include <iostream>

class SimParameters
{
public:
   std::vector<double> getParameters()const;

   SimParameters& applyParameters(const std::vector<double>& param);

   std::string mode_;

   double max_num_cells_;
   double init_ratio_APC_cells_;
   double init_ratio_NK_cells_;
   double init_ratio_LT_cells_;
   double LT_ratio_specific_;


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

   SimParameters(const SimParameters& other);

   friend void swap(SimParameters& one, SimParameters& other);

   friend std::ostream& operator<<(std::ostream& s,SimParameters p);

   SimParameters& operator=(const SimParameters& other);

   SimParameters();

   void reset(const SimParameters& sp);
   ~SimParameters(){}
};



#endif // SIMPARAMETERS_H_INCLUDED
