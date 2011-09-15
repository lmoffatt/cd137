#include "Includes/SimParameters.h"
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cstring>

SimParameters::SimParameters():
    mode_("FULL"),
    max_num_cells_(2e6),
    init_ratio_APC_cells_ (1e5/1e6),
    init_ratio_NK_cells_(1e5/1e6),
    init_ratio_LT_cells_ (7.99e5/1e6),
    LT_ratio_specific_ (1000/1e6),

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
    LT_exh_rate_(1.0/1e4),

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
        par.push_back(10*log(max_num_cells_)); // multiplicar si estoy seguro del valor
        par.push_back(10*log(init_ratio_APC_cells_));
        par.push_back(10*log(init_ratio_NK_cells_));
        par.push_back(10*log(init_ratio_LT_cells_));
        par.push_back(log(LT_ratio_specific_));


        // par.push_back(log(Ag_internalization_rate));

        par.push_back(log(APC_max_proliferation_rate_));
        par.push_back(log(NK_max_proliferation_rate_));
        par.push_back(log(LT_max_no_receptor_prol_rate_));
        par.push_back(log(LT_max_free_prol_rate_));
        par.push_back(log(LT_max_bound_prol_rate_));
        par.push_back(log(LT_max_blocked_prol_rate_));

        par.push_back(log(APC_no_to_free_rate_per_Ag_));
        par.push_back(log(APC_free_to_bound_rate_per_LT_));
        par.push_back(log(APC_Ab_binding_rate_));
        par.push_back(log(APC_exh_rate));
        par.push_back(log(NK_no_to_free_rate_per_Ag_));
        par.push_back(log(NK_free_to_bound_rate_per_LT_));
        par.push_back(log(NK_Ab_binding_rate));
        par.push_back(log(NK_exh_rate));
        par.push_back(log(LT_no_to_free_rate_per_APC_));
        par.push_back(log(LT_free_to_bound_rate_per_APC_));
        par.push_back(log(LT_mAb_binding_rate_));
        par.push_back(log(LT_exh_rate_));

        par.push_back(log(APC_IFN_free_prod_rate_));
        par.push_back(log(APC_IFN_Ag_prod_rate_));
        par.push_back(log(APC_IFN_bound_prod_rate_));
        par.push_back(log(APC_IFN_blocked_prod_rate_));
        par.push_back(log(NK_IFN_free_prod_rate_));
        par.push_back(log(NK_IFN_Ag_prod_rate_));
        par.push_back(log(NK_IFN_bound_prod_rate_));
        par.push_back(log(NK_IFN_blocked_prod_rate_));
        par.push_back(log(LT_IFN_no_rec_prod_rate_));
        par.push_back(log(LT_IFN_free_prod_rate_));
        par.push_back(log(LT_IFN_bound_prod_rate_));
        par.push_back(log(LT_IFN_blocked_prod_rate_));

        par.push_back(log(APC_TNF_free_prod_rate_));
        par.push_back(log(APC_TNF_Ag_prod_rate_));
        par.push_back(log(APC_TNF_bound_prod_rate_));
        par.push_back(log(APC_TNF_blocked_prod_rate_));
        par.push_back(log(NK_TNF_free_prod_rate_));
        par.push_back(log(NK_TNF_Ag_prod_rate_));
        par.push_back(log(NK_TNF_bound_prod_rate_));
        par.push_back(log(NK_TNF_blocked_prod_rate_));
        par.push_back(log(LT_TNF_no_rec_prod_rate_));
        par.push_back(log(LT_TNF_free_prod_rate_));
        par.push_back(log(LT_TNF_bound_prod_rate_));
        par.push_back(log(LT_TNF_blocked_prod_rate_));
        par.push_back(log(TNF_deg));
        par.push_back(log(IFN_deg));
    }

    else if (mode_=="PARTIAL")
    {
        par.push_back(10*log(max_num_cells_)); // multiplicar si estoy seguro del valor
        par.push_back(10*log(init_ratio_APC_cells_));
        par.push_back(10*log(init_ratio_NK_cells_));
        par.push_back(10*log(init_ratio_LT_cells_));
        par.push_back(log(LT_ratio_specific_));


        // par.push_back(log(Ag_internalization_rate));

        par.push_back(log(APC_max_proliferation_rate_));
        par.push_back(log(NK_max_proliferation_rate_));
        par.push_back(log(LT_max_no_receptor_prol_rate_));
        par.push_back(log(LT_max_free_prol_rate_));
        par.push_back(log(LT_max_bound_prol_rate_));
        par.push_back(log(LT_max_blocked_prol_rate_));

        par.push_back(log(APC_no_to_free_rate_per_Ag_));
        par.push_back(log(APC_free_to_bound_rate_per_LT_));
        par.push_back(log(APC_Ab_binding_rate_));
        par.push_back(log(APC_exh_rate));
        par.push_back(log(NK_no_to_free_rate_per_Ag_));
        par.push_back(log(NK_free_to_bound_rate_per_LT_));
        par.push_back(log(NK_Ab_binding_rate));
        par.push_back(log(NK_exh_rate));
        par.push_back(log(LT_no_to_free_rate_per_APC_));
        par.push_back(log(LT_free_to_bound_rate_per_APC_));
        par.push_back(log(LT_mAb_binding_rate_));
        par.push_back(log(LT_exh_rate_));

        par.push_back(log(APC_IFN_free_prod_rate_));
        par.push_back(log(APC_IFN_Ag_prod_rate_));
        par.push_back(log(APC_IFN_bound_prod_rate_));
        //par.push_back(log(APC_IFN_blocked_prod_rate_));
        par.push_back(log(NK_IFN_free_prod_rate_));
        par.push_back(log(NK_IFN_Ag_prod_rate_));
        par.push_back(log(NK_IFN_bound_prod_rate_));
        //par.push_back(log(NK_IFN_blocked_prod_rate_));
        par.push_back(log(LT_IFN_no_rec_prod_rate_));
        par.push_back(log(LT_IFN_free_prod_rate_));
        par.push_back(log(LT_IFN_bound_prod_rate_));
        par.push_back(log(LT_IFN_blocked_prod_rate_));

        par.push_back(log(APC_TNF_free_prod_rate_));
        par.push_back(log(APC_TNF_Ag_prod_rate_));
        par.push_back(log(APC_TNF_bound_prod_rate_));
        //par.push_back(log(APC_TNF_blocked_prod_rate_));
        par.push_back(log(NK_TNF_free_prod_rate_));
        par.push_back(log(NK_TNF_Ag_prod_rate_));
        par.push_back(log(NK_TNF_bound_prod_rate_));
        //par.push_back(log(NK_TNF_blocked_prod_rate_));
        par.push_back(log(LT_TNF_no_rec_prod_rate_));
        par.push_back(log(LT_TNF_free_prod_rate_));
        par.push_back(log(LT_TNF_bound_prod_rate_));
        par.push_back(log(LT_TNF_blocked_prod_rate_));



        par.push_back(log(TNF_deg));
        par.push_back(log(IFN_deg));
    }

    return par;


}

SimParameters& SimParameters::applyParameters(const std::vector<double>& param)
{
    if (mode_=="FULL")
    {
        std::size_t i=0;
        max_num_cells_=exp(param[i++]/10);//Divdir si estoy seguro de su valor
        init_ratio_APC_cells_=exp(param[i++]/10);
        init_ratio_NK_cells_=exp(param[i++]/10);
        init_ratio_LT_cells_=exp(param[i++]/10);
        LT_ratio_specific_=exp(param[i++]);


        // Ag_internalization_rate=exp(param[0]);

        APC_max_proliferation_rate_=exp(param[i++]);
        NK_max_proliferation_rate_=exp(param[i++]);
        LT_max_no_receptor_prol_rate_=exp(param[i++]);
        LT_max_free_prol_rate_=exp(param[i++]);
        LT_max_bound_prol_rate_=exp(param[i++]);
        LT_max_blocked_prol_rate_=exp(param[i++]);


        APC_no_to_free_rate_per_Ag_=exp(param[i++]);
        APC_free_to_bound_rate_per_LT_=exp(param[i++]);
        APC_Ab_binding_rate_=exp(param[i++]);
        APC_exh_rate=exp(param[i++]);
        NK_no_to_free_rate_per_Ag_=exp(param[i++]);
        NK_free_to_bound_rate_per_LT_=exp(param[i++]);
        NK_Ab_binding_rate=exp(param[i++]);
        NK_exh_rate=exp(param[i++]);
        LT_no_to_free_rate_per_APC_=exp(param[i++]);
        LT_free_to_bound_rate_per_APC_=exp(param[i++]);
        LT_mAb_binding_rate_=exp(param[i++]);
        LT_exh_rate_=exp (param[i++]);

        APC_IFN_free_prod_rate_=exp(param[i++]);
        APC_IFN_Ag_prod_rate_=exp(param[i++]);
        APC_IFN_bound_prod_rate_=exp(param[i++]);
        APC_IFN_blocked_prod_rate_=exp(param[i++]);
        NK_IFN_free_prod_rate_=exp(param[i++]);
        NK_IFN_Ag_prod_rate_=exp(param[i++]);
        NK_IFN_bound_prod_rate_=exp(param[i++]);
        NK_IFN_blocked_prod_rate_=exp(param[i++]);
        LT_IFN_no_rec_prod_rate_=exp(param[i++]);
        LT_IFN_free_prod_rate_=exp(param[i++]);
        LT_IFN_bound_prod_rate_=exp(param[i++]);
        LT_IFN_blocked_prod_rate_=exp(param[i++]);

        APC_TNF_free_prod_rate_=exp(param[i++]);
        APC_TNF_Ag_prod_rate_=exp(param[i++]);
        APC_TNF_bound_prod_rate_=exp(param[i++]);
        APC_TNF_blocked_prod_rate_=exp(param[i++]);
        NK_TNF_free_prod_rate_=exp(param[i++]);
        NK_TNF_Ag_prod_rate_=exp(param[i++]);
        NK_TNF_bound_prod_rate_=exp(param[i++]);
        NK_TNF_blocked_prod_rate_=exp(param[i++]);
        LT_TNF_no_rec_prod_rate_=exp(param[i++]);
        LT_TNF_free_prod_rate_=exp(param[i++]);
        LT_TNF_bound_prod_rate_=exp(param[i++]);
        LT_TNF_blocked_prod_rate_=exp(param[i++]);
        TNF_deg=exp(param[i++]);
        IFN_deg=exp(param[i++]);
    }
    else if(mode_=="PARTIAL")
    {
        std::size_t i=0;
        max_num_cells_=exp(param[i++]/10);
        init_ratio_APC_cells_=exp(param[i++]/10);
        init_ratio_NK_cells_=exp(param[i++]/10);
        init_ratio_LT_cells_=exp(param[i++]/10);
        LT_ratio_specific_=exp(param[i++]);


        // Ag_internalization_rate=exp(param[0]);

        APC_max_proliferation_rate_=exp(param[i++]);
        NK_max_proliferation_rate_=exp(param[i++]);
        LT_max_no_receptor_prol_rate_=exp(param[i++]);
        LT_max_free_prol_rate_=exp(param[i++]);
        LT_max_bound_prol_rate_=exp(param[i++]);
        LT_max_blocked_prol_rate_=exp(param[i++]);


        APC_no_to_free_rate_per_Ag_=exp(param[i++]);
        APC_free_to_bound_rate_per_LT_=exp(param[i++]);
        APC_Ab_binding_rate_=exp(param[i++]);
        APC_exh_rate=exp(param[i++]);
        NK_no_to_free_rate_per_Ag_=exp(param[i++]);
        NK_free_to_bound_rate_per_LT_=exp(param[i++]);
        NK_Ab_binding_rate=exp(param[i++]);
        NK_exh_rate=exp(param[i++]);
        LT_no_to_free_rate_per_APC_=exp(param[i++]);
        LT_free_to_bound_rate_per_APC_=exp(param[i++]);
        LT_mAb_binding_rate_=exp(param[i++]);
        LT_exh_rate_=exp (param[i++]);

        APC_IFN_free_prod_rate_=exp(param[i++]);
        APC_IFN_Ag_prod_rate_=exp(param[i++]);
        APC_IFN_bound_prod_rate_=exp(param[i++]);
        // APC_IFN_blocked_prod_rate_=exp(param[i++]);
        APC_IFN_blocked_prod_rate_=APC_IFN_Ag_prod_rate_;
        NK_IFN_free_prod_rate_=exp(param[i++]);
        NK_IFN_Ag_prod_rate_=exp(param[i++]);
        NK_IFN_bound_prod_rate_=exp(param[i++]);
        //NK_IFN_blocked_prod_rate_=exp(param[i++]);
        NK_IFN_blocked_prod_rate_=NK_IFN_free_prod_rate_;
        LT_IFN_no_rec_prod_rate_=exp(param[i++]);
        LT_IFN_free_prod_rate_=exp(param[i++]);
        LT_IFN_bound_prod_rate_=exp(param[i++]);
        LT_IFN_blocked_prod_rate_=exp(param[i++]);

        APC_TNF_free_prod_rate_=exp(param[i++]);
        APC_TNF_Ag_prod_rate_=exp(param[i++]);
        APC_TNF_bound_prod_rate_=exp(param[i++]);
        //APC_TNF_blocked_prod_rate_=exp(param[i++]);
        APC_TNF_blocked_prod_rate_=APC_TNF_Ag_prod_rate_;
        NK_TNF_free_prod_rate_=exp(param[i++]);
        NK_TNF_Ag_prod_rate_=exp(param[i++]);
        NK_TNF_bound_prod_rate_=exp(param[i++]);
        //NK_TNF_blocked_prod_rate_=exp(param[i++]);
        NK_TNF_blocked_prod_rate_=NK_TNF_Ag_prod_rate_;
        LT_TNF_no_rec_prod_rate_=exp(param[i++]);
        LT_TNF_free_prod_rate_=exp(param[i++]);
        LT_TNF_bound_prod_rate_=exp(param[i++]);
        LT_TNF_no_rec_prod_rate_=exp(param[i++]);
        LT_TNF_free_prod_rate_=exp(param[i++]);
        LT_TNF_blocked_prod_rate_=exp(param[i++]);



        TNF_deg=exp(param[i++]);
        IFN_deg=exp(param[i++]);


    }
    return *this;
}





SimParameters::SimParameters(const SimParameters& other):
    mode_(other.mode_),
    max_num_cells_(other.max_num_cells_),
    init_ratio_APC_cells_(other.init_ratio_APC_cells_),
    init_ratio_NK_cells_ (other.init_ratio_NK_cells_),
    init_ratio_LT_cells_(other.init_ratio_LT_cells_),
    LT_ratio_specific_ (other.LT_ratio_specific_),


    //  (Ag_internalization_rate),

    APC_max_proliferation_rate_  (other.APC_max_proliferation_rate_),
    NK_max_proliferation_rate_ (other.NK_max_proliferation_rate_),
    LT_max_no_receptor_prol_rate_ (other.LT_max_no_receptor_prol_rate_),
    LT_max_free_prol_rate_ (other.LT_max_free_prol_rate_),
    LT_max_bound_prol_rate_ (other.LT_max_bound_prol_rate_),
    LT_max_blocked_prol_rate_ (other.LT_max_blocked_prol_rate_),

    APC_no_to_free_rate_per_Ag_ (other.APC_no_to_free_rate_per_Ag_),
    APC_free_to_bound_rate_per_LT_ (other.APC_free_to_bound_rate_per_LT_),
    APC_Ab_binding_rate_ (other.APC_Ab_binding_rate_),
    APC_exh_rate (other.APC_exh_rate),
    NK_no_to_free_rate_per_Ag_ (other.NK_no_to_free_rate_per_Ag_),
    NK_free_to_bound_rate_per_LT_ (other.NK_free_to_bound_rate_per_LT_),
    NK_Ab_binding_rate (other.NK_Ab_binding_rate),
    NK_exh_rate  (other.NK_exh_rate),
    LT_no_to_free_rate_per_APC_ (other.LT_no_to_free_rate_per_APC_),
    LT_free_to_bound_rate_per_APC_ (other.LT_free_to_bound_rate_per_APC_),
    LT_mAb_binding_rate_ (other.LT_mAb_binding_rate_),
    LT_exh_rate_ (other.LT_exh_rate_),

    APC_IFN_free_prod_rate_  (other.APC_IFN_free_prod_rate_),
    APC_IFN_Ag_prod_rate_ (other.APC_IFN_Ag_prod_rate_),
    APC_IFN_bound_prod_rate_  (other.APC_IFN_bound_prod_rate_),
    APC_IFN_blocked_prod_rate_ (other.APC_IFN_blocked_prod_rate_),
    NK_IFN_free_prod_rate_ (other.NK_IFN_free_prod_rate_),
    NK_IFN_Ag_prod_rate_(other.NK_IFN_Ag_prod_rate_),
    NK_IFN_bound_prod_rate_(other.NK_IFN_bound_prod_rate_),
    NK_IFN_blocked_prod_rate_(other.NK_IFN_blocked_prod_rate_),
    LT_IFN_no_rec_prod_rate_ (other.LT_IFN_no_rec_prod_rate_),
    LT_IFN_free_prod_rate_(other.LT_IFN_free_prod_rate_),
    LT_IFN_bound_prod_rate_(other.LT_IFN_bound_prod_rate_),
    LT_IFN_blocked_prod_rate_(other.LT_IFN_blocked_prod_rate_),

    APC_TNF_free_prod_rate_(other.APC_TNF_free_prod_rate_),
    APC_TNF_Ag_prod_rate_(other.APC_TNF_Ag_prod_rate_),
    APC_TNF_bound_prod_rate_(other.APC_TNF_bound_prod_rate_),
    APC_TNF_blocked_prod_rate_(other.APC_TNF_blocked_prod_rate_),
    NK_TNF_free_prod_rate_(other.NK_TNF_free_prod_rate_),
    NK_TNF_Ag_prod_rate_(other.NK_TNF_Ag_prod_rate_),
    NK_TNF_bound_prod_rate_(other.NK_TNF_bound_prod_rate_),
    NK_TNF_blocked_prod_rate_(other.NK_TNF_blocked_prod_rate_),
    LT_TNF_no_rec_prod_rate_(other.LT_TNF_no_rec_prod_rate_),
    LT_TNF_free_prod_rate_(other.LT_TNF_free_prod_rate_),
    LT_TNF_bound_prod_rate_(other.LT_TNF_bound_prod_rate_),
    LT_TNF_blocked_prod_rate_(other.LT_TNF_blocked_prod_rate_),



    TNF_deg (other.TNF_deg),
    IFN_deg (other.IFN_deg)
{}

SimParameters&
SimParameters::operator=(const SimParameters& other)
{
    if (this!=&other)
    {
        SimParameters tmp(other);
        swap(*this,tmp);
    }
    return *this;
}

void swap(SimParameters& one, SimParameters& other)
{
    std::swap(one.mode_,other.mode_);
    std::swap(one.max_num_cells_,other.max_num_cells_);
    std::swap(one.init_ratio_APC_cells_,other.init_ratio_APC_cells_);
    std::swap(one.init_ratio_NK_cells_ ,other.init_ratio_NK_cells_);
    std::swap(one.init_ratio_LT_cells_,other.init_ratio_LT_cells_);
    std::swap(one.LT_ratio_specific_ ,other.LT_ratio_specific_);


    //  (Ag_internalization_rate),

    std::swap(one.APC_max_proliferation_rate_  ,other.APC_max_proliferation_rate_);
    std::swap(one.NK_max_proliferation_rate_ ,other.NK_max_proliferation_rate_);
    std::swap(one.LT_max_no_receptor_prol_rate_ ,other.LT_max_no_receptor_prol_rate_);
    std::swap(one.LT_max_free_prol_rate_ ,other.LT_max_free_prol_rate_);
    std::swap(one.LT_max_bound_prol_rate_ ,other.LT_max_bound_prol_rate_);
    std::swap(one.LT_max_blocked_prol_rate_ ,other.LT_max_blocked_prol_rate_);

    std::swap(one.APC_no_to_free_rate_per_Ag_ ,other.APC_no_to_free_rate_per_Ag_);
    std::swap(one.APC_free_to_bound_rate_per_LT_ ,other.APC_free_to_bound_rate_per_LT_);
    std::swap(one.APC_Ab_binding_rate_ ,other.APC_Ab_binding_rate_);
    std::swap(one.APC_exh_rate ,other.APC_exh_rate);
    std::swap(one.NK_no_to_free_rate_per_Ag_ ,other.NK_no_to_free_rate_per_Ag_);
    std::swap(one.NK_free_to_bound_rate_per_LT_ ,other.NK_free_to_bound_rate_per_LT_);
    std::swap(one.NK_Ab_binding_rate ,other.NK_Ab_binding_rate);
    std::swap(one.NK_exh_rate  ,other.NK_exh_rate);
    std::swap(one.LT_no_to_free_rate_per_APC_ ,other.LT_no_to_free_rate_per_APC_);
    std::swap(one.LT_free_to_bound_rate_per_APC_ ,other.LT_free_to_bound_rate_per_APC_);
    std::swap(one.LT_mAb_binding_rate_ ,other.LT_mAb_binding_rate_);
    std::swap(one.LT_exh_rate_,other.LT_exh_rate_);

    std::swap(one.APC_IFN_free_prod_rate_  ,other.APC_IFN_free_prod_rate_);
    std::swap(one.APC_IFN_Ag_prod_rate_ ,other.APC_IFN_Ag_prod_rate_);
    std::swap(one.APC_IFN_bound_prod_rate_ ,other.APC_IFN_bound_prod_rate_);
    std::swap(one.APC_IFN_blocked_prod_rate_ ,other.APC_IFN_blocked_prod_rate_);
    std::swap(one.NK_IFN_free_prod_rate_ ,other.NK_IFN_free_prod_rate_);
    std::swap(one.NK_IFN_Ag_prod_rate_,other.NK_IFN_Ag_prod_rate_);
    std::swap(one.NK_IFN_bound_prod_rate_,other.NK_IFN_bound_prod_rate_);
    std::swap(one.NK_IFN_blocked_prod_rate_,other.NK_IFN_blocked_prod_rate_);
    std::swap(one.LT_IFN_no_rec_prod_rate_ ,other.LT_IFN_no_rec_prod_rate_);
    std::swap(one.LT_IFN_free_prod_rate_,other.LT_IFN_free_prod_rate_);
    std::swap(one.LT_IFN_bound_prod_rate_,other.LT_IFN_bound_prod_rate_);
    std::swap(one.LT_IFN_blocked_prod_rate_,other.LT_IFN_blocked_prod_rate_);

    std::swap(one.APC_TNF_free_prod_rate_,other.APC_TNF_free_prod_rate_);
    std::swap(one.APC_TNF_Ag_prod_rate_,other.APC_TNF_Ag_prod_rate_);
    std::swap(one.APC_TNF_bound_prod_rate_,other.APC_TNF_bound_prod_rate_);
    std::swap(one.APC_TNF_blocked_prod_rate_,other.APC_TNF_blocked_prod_rate_);
    std::swap(one.NK_TNF_free_prod_rate_,other.NK_TNF_free_prod_rate_);
    std::swap(one.NK_TNF_Ag_prod_rate_,other.NK_TNF_Ag_prod_rate_);
    std::swap(one.NK_TNF_bound_prod_rate_,other.NK_TNF_bound_prod_rate_);
    std::swap(one.NK_TNF_blocked_prod_rate_,other.NK_TNF_blocked_prod_rate_);
    std::swap(one.LT_TNF_no_rec_prod_rate_,other.LT_TNF_no_rec_prod_rate_);
    std::swap(one.LT_TNF_free_prod_rate_,other.LT_TNF_free_prod_rate_);
    std::swap(one.LT_TNF_bound_prod_rate_,other.LT_TNF_bound_prod_rate_);
    std::swap(one.LT_TNF_blocked_prod_rate_,other.LT_TNF_blocked_prod_rate_);



    std::swap(one.TNF_deg,other.TNF_deg);
    std::swap(one.IFN_deg,other.IFN_deg);
}

std::ostream& operator<<(std::ostream& s,SimParameters p)
{
    s<<"\n max_num_cells_ \t"<<p.max_num_cells_;
    s<<"\n init_ratio_APC_cells_ \t"<<p.init_ratio_APC_cells_;
    s<<"\n init_ratio_NK_cells_ \t"<<p.init_ratio_NK_cells_;
    s<<"\n init_ratio_LT_cells_ \t"<<p.init_ratio_LT_cells_;
    s<<"\n LT_ratio_specific_ \t"<<p.LT_ratio_specific_;
    s<<"\n APC_max_proliferation_rate_ \t"<<p.APC_max_proliferation_rate_;
    s<<"\n NK_max_proliferation_rate_ \t"<<p.NK_max_proliferation_rate_;
    s<<"\n LT_max_no_receptor_prol_rate_ \t"<<p.LT_max_no_receptor_prol_rate_;
    s<<"\n LT_max_free_prol_rate_ \t"<<p.LT_max_free_prol_rate_;
    s<<"\n LT_max_bound_prol_rate_ \t"<<p.LT_max_bound_prol_rate_;

    s<<"\n LT_max_blocked_prol_rate_ \t"<<p.LT_max_blocked_prol_rate_;
    s<<"\n APC_no_to_free_rate_per_Ag_ \t"<<p.APC_no_to_free_rate_per_Ag_;
    s<<"\n APC_free_to_bound_rate_per_LT_ \t"<<p.APC_free_to_bound_rate_per_LT_;
    s<<"\n APC_Ab_binding_rate_ \t"<<p.APC_Ab_binding_rate_;
    s<<"\n APC_exh_rate \t"<<p.APC_exh_rate;
    s<<"\n NK_no_to_free_rate_per_Ag_ \t"<<p.NK_no_to_free_rate_per_Ag_;
    s<<"\n NK_free_to_bound_rate_per_LT_ \t"<<p.NK_free_to_bound_rate_per_LT_;
    s<<"\n NK_Ab_binding_rate \t"<<p.NK_Ab_binding_rate;
    s<<"\n NK_exh_rate \t"<<p.NK_exh_rate;
    s<<"\n LT_no_to_free_rate_per_APC_ \t"<<p.LT_no_to_free_rate_per_APC_;
    s<<"\n LT_exh_rate_\t"<<p.LT_exh_rate_;
    s<<"\n LT_free_to_bound_rate_per_APC_ \t"<<p.LT_free_to_bound_rate_per_APC_;
    s<<"\n LT_mAb_binding_rate_ \t"<<p.LT_mAb_binding_rate_;
    s<<"\n APC_IFN_free_prod_rate_ \t"<<p.APC_IFN_free_prod_rate_;
    s<<"\n APC_IFN_Ag_prod_rate_ \t"<<p.APC_IFN_Ag_prod_rate_;
    s<<"\n APC_IFN_bound_prod_rate_ \t"<<p.APC_IFN_bound_prod_rate_;
    s<<"\n APC_IFN_blocked_prod_rate_ \t"<<p.APC_IFN_blocked_prod_rate_;
    s<<"\n NK_IFN_free_prod_rate_ \t"<<p.NK_IFN_free_prod_rate_;
    s<<"\n NK_IFN_Ag_prod_rate_ \t"<<p.NK_IFN_Ag_prod_rate_;
    s<<"\n NK_IFN_bound_prod_rate_ \t"<<p.NK_IFN_bound_prod_rate_;
    s<<"\n NK_IFN_blocked_prod_rate_ \t"<<p.NK_IFN_blocked_prod_rate_;
    s<<"\n LT_IFN_no_rec_prod_rate_ \t"<<p.LT_IFN_no_rec_prod_rate_;
    s<<"\n LT_IFN_free_prod_rate_ \t"<<p.LT_IFN_free_prod_rate_;
    s<<"\n LT_IFN_bound_prod_rate_ \t"<<p.LT_IFN_bound_prod_rate_;

    s<<"\n LT_IFN_blocked_prod_rate_ \t"<<p.LT_IFN_blocked_prod_rate_;
    s<<"\n APC_TNF_free_prod_rate_ \t"<<p.APC_TNF_free_prod_rate_;
    s<<"\n APC_TNF_Ag_prod_rate_ \t"<<p.APC_TNF_Ag_prod_rate_;
    s<<"\n APC_TNF_bound_prod_rate_ \t"<<p.APC_TNF_bound_prod_rate_;
    s<<"\n APC_TNF_blocked_prod_rate_ \t"<<p.APC_TNF_blocked_prod_rate_;
    s<<"\n NK_TNF_free_prod_rate_ \t"<<p.NK_TNF_free_prod_rate_;
    s<<"\n NK_TNF_Ag_prod_rate_ \t"<<p.NK_TNF_Ag_prod_rate_;
    s<<"\n NK_TNF_bound_prod_rate_ \t"<<p.NK_TNF_bound_prod_rate_;
    s<<"\n NK_TNF_blocked_prod_rate_ \t"<<p.NK_TNF_blocked_prod_rate_;

    s<<"\n LT_TNF_no_rec_prod_rate_ \t"<<p.LT_TNF_no_rec_prod_rate_;
    s<<"\n LT_TNF_free_prod_rate_ \t"<<p.LT_TNF_free_prod_rate_;
    s<<"\n LT_TNF_bound_prod_rate_ \t"<<p.LT_TNF_bound_prod_rate_;
    s<<"\n LT_TNF_blocked_prod_rate_ \t"<<p.LT_TNF_blocked_prod_rate_;
    s<<"\n TNF_deg \t"<<p.TNF_deg;
    s<<"\n IFN_deg \t"<<p.IFN_deg;

    return s;

}
std::vector<double> SimParameters::getRandomParameters(double range)const
{
    std::vector<double> result=getParameters();
    srand ( time(NULL) );

    for (std::size_t i=0; i<result.size();++i)
    {
        result[i]+=2*range*((1.0*rand())/RAND_MAX-0.5);
    }
    return result;
}
