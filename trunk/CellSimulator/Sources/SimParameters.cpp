#include "Includes/SimParameters.h"


SimParameters::SimParameters():

    mode_("FULL"),
    max_num_cells_(2e6),
    init_ratio_APC_cells_ (1e5),
    init_ratio_NK_cells_(1e5),
    init_ratio_LT_cells_ (9e5),
    LT_ratio_specific_ (1000),
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
    {}

std::vector<double> SimParameters::getParameters()const
{
    std::vector<double> par;

    if (mode_=="FULL")
    {
	par.push_back(max_num_cells_);
	par.push_back(init_ratio_APC_cells_);
	par.push_back(init_ratio_NK_cells_);
	par.push_back(init_ratio_LT_cells_);
	par.push_back(LT_ratio_specific_);


	// par.push_back(Ag_internalization_rate);

	par.push_back(APC_max_proliferation_rate_);
	par.push_back(NK_max_proliferation_rate_);
	par.push_back(LT_max_no_receptor_prol_rate_);
	par.push_back(LT_max_free_prol_rate_);
	par.push_back(LT_max_bound_prol_rate_);
	par.push_back(LT_max_blocked_prol_rate_);

	par.push_back(APC_no_to_free_rate_per_Ag_);
	par.push_back(APC_free_to_bound_rate_per_LT_);
	par.push_back(APC_Ab_binding_rate_);
	par.push_back(APC_exh_rate);
	par.push_back(NK_no_to_free_rate_per_Ag_);
	par.push_back(NK_free_to_bound_rate_per_LT_);
	par.push_back(NK_Ab_binding_rate);
	par.push_back(NK_exh_rate);
	par.push_back(LT_no_to_free_rate_per_APC_);
	par.push_back(LT_free_to_bound_rate_per_APC_);
	par.push_back(LT_mAb_binding_rate_);

	par.push_back(APC_IFN_free_prod_rate_);
	par.push_back(APC_IFN_Ag_prod_rate_);
	par.push_back(APC_IFN_bound_prod_rate_);
	par.push_back(APC_IFN_blocked_prod_rate_);
	par.push_back(NK_IFN_free_prod_rate_);
	par.push_back(NK_IFN_Ag_prod_rate_);
	par.push_back(NK_IFN_bound_prod_rate_);
	par.push_back(NK_IFN_blocked_prod_rate_);
	par.push_back(LT_IFN_no_rec_prod_rate_);
	par.push_back(LT_IFN_free_prod_rate_);
	par.push_back(LT_IFN_bound_prod_rate_);
	par.push_back(LT_IFN_blocked_prod_rate_);

	par.push_back(APC_TNF_free_prod_rate_);
	par.push_back(APC_TNF_Ag_prod_rate_);
	par.push_back(APC_TNF_bound_prod_rate_);
	par.push_back(APC_TNF_blocked_prod_rate_);
	par.push_back(NK_TNF_free_prod_rate_);
	par.push_back(NK_TNF_Ag_prod_rate_);
	par.push_back(NK_TNF_bound_prod_rate_);
	par.push_back(NK_TNF_blocked_prod_rate_);
	par.push_back(LT_TNF_no_rec_prod_rate_);
	par.push_back(LT_TNF_free_prod_rate_);
	par.push_back(LT_TNF_bound_prod_rate_);
	par.push_back(LT_TNF_blocked_prod_rate_);



	par.push_back(TNF_deg);
	par.push_back(IFN_deg);
    }
    else if(mode_=="PARTIAL")
    {

    }

    return par;


}

SimParameters& SimParameters::applyParameters(const std::vector<double>& param)
{

    max_num_cells_=param[0];
    init_ratio_APC_cells_=param[1];
    init_ratio_NK_cells_=param[2];
    init_ratio_LT_cells_=param[3];
    LT_ratio_specific_=param[4];


   // Ag_internalization_rate=param[0];

    APC_max_proliferation_rate_=param[5];
    NK_max_proliferation_rate_=param[6];
    LT_max_no_receptor_prol_rate_=param[7];
    LT_max_free_prol_rate_=param[8];
    LT_max_bound_prol_rate_=param[9];
    LT_max_blocked_prol_rate_=param[10];

    APC_no_to_free_rate_per_Ag_=param[11];
    APC_free_to_bound_rate_per_LT_=param[12];
    APC_Ab_binding_rate_=param[13];
    APC_exh_rate=param[14];
    NK_no_to_free_rate_per_Ag_=param[15];
    NK_free_to_bound_rate_per_LT_=param[16];
    NK_Ab_binding_rate=param[17];
    NK_exh_rate=param[18];
    LT_no_to_free_rate_per_APC_=param[19];
    LT_free_to_bound_rate_per_APC_=param[20];
    LT_mAb_binding_rate_=param[21];

    APC_IFN_free_prod_rate_=param[22];
    APC_IFN_Ag_prod_rate_=param[23];
    APC_IFN_bound_prod_rate_=param[24];
    APC_IFN_blocked_prod_rate_=param[25];
    NK_IFN_free_prod_rate_=param[26];
    NK_IFN_Ag_prod_rate_=param[27];
    NK_IFN_bound_prod_rate_=param[28];
    NK_IFN_blocked_prod_rate_=param[29];
    LT_IFN_no_rec_prod_rate_=param[30];
    LT_IFN_free_prod_rate_=param[31];
    LT_IFN_bound_prod_rate_=param[32];
    LT_IFN_blocked_prod_rate_=param[33];

    APC_TNF_free_prod_rate_=param[34];
    APC_TNF_Ag_prod_rate_=param[35];
    APC_TNF_bound_prod_rate_=param[36];
    APC_TNF_blocked_prod_rate_=param[37];
    NK_TNF_free_prod_rate_=param[38];
    NK_TNF_Ag_prod_rate_=param[39];
    NK_TNF_bound_prod_rate_=param[40];
    NK_TNF_blocked_prod_rate_=param[41];
    LT_TNF_no_rec_prod_rate_=param[42];
    LT_TNF_free_prod_rate_=param[43];
    LT_TNF_bound_prod_rate_=param[44];
    LT_TNF_blocked_prod_rate_=param[45];



    TNF_deg=param[46];
    IFN_deg=param[47];

return *this;
}
