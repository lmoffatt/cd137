#include <cmath>
#include<limits>
#include "Includes/LevenbergMarquardt.h"
#include "Includes/MatrixInverse.h"
#include "LevenbergMarquardtParameters.h"




LevenbergMarquardtParameters::LevenbergMarquardtParameters(
    ABC_function* f,
    const std::vector<double>& data,
    const Parameters& initialParam,
    const std::vector<double>& weigth):
    f_(f),
    data_(data),
    w_(weigth),
    initialParam_(initialParam),
    nPar_(initialParam_.size()),
    nData_(data.size()),
    dx_(1e-9),
    maxIter_(30),
    maxFeval_(10000),
    minParamChange_(1e-12),
    minSSChange_(1e-12),
    minGradient_(0.0000001),
    maxLanda_(1e11),
    landa_(1000),
    v_(10),
    nIter_(0),
    nFeval_(0),
    currSS_(std::numeric_limits<double>::infinity()),
    newSS0_(std::numeric_limits<double>::infinity()),
    currParam_(initialParam_),
    newParam_(currParam_),
    newParam0_(currParam_),
    currYfit_(std::vector<double>(nData_)),
    newYfit_(currYfit_),
    newYfit0_(currYfit_),
    J_(std::vector< std::vector< double> >(nData_,std::vector<double>(nPar_))),
    G_(currParam_.size()),
    JTJ_(std::vector< std::vector<double> > (nPar_,std::vector<double>(nPar_))),
    JTJinv_(JTJ_),
    d_(currParam_.size()),
    optimParam_(currParam_),
    surpassIter_(false),
    surpassFeval_(false),
    ParamChange_(std::numeric_limits<double>::infinity()),
    SSChange_(ParamChange_),
    NormGrad_(ParamChange_),
    smallParamChange_(false),
    smallSSChange_(false),
    smallGradient_(false){}



LevenbergMarquardtParameters::LevenbergMarquardtParameters(){}

LevenbergMarquardtParameters::LevenbergMarquardtParameters (const LevenbergMarquardtParameters& other):
    data_(other.data_),
    w_(other.w_),
    initialParam_(other.initialParam_),
    nPar_(other.nPar_),
    nData_(other.nData_),
    dx_(other.dx_),
    maxIter_(other.maxIter_),
    maxFeval_(other.maxFeval_),
    minParamChange_(other.minParamChange_),
    minSSChange_(other.minSSChange_),
    minGradient_(other.minGradient_),
    maxLanda_(other.maxLanda_),
    landa_(other.landa_),
    landa0_(other.landa0_),
    v_(other.v_),
    nIter_(other.nIter_),
    nFeval_(other.nFeval_),



    currSS_(other.currSS_),
    newSS0_ (other.newSS0_),
    currParam_ (other.currParam_),
    newParam_ (other.newParam_),
    newParam0_ (other.newParam0_),
            currYfit_(other.currYfit_),
            newYfit_(other.newYfit_),
    newYfit0_(other.newYfit0_),

            J_(other.J_),
            G_(other.G_),
            JTJ_(other.JTJ_),
            JTJinv_(other.JTJinv_),

            d_(other.d_),

            optimParam_(other.optimParam_),

            surpassIter_(other.surpassIter_),
            surpassFeval_(other.surpassFeval_),

            ParamChange_(other.ParamChange_),
            SSChange_(other.SSChange_),
            NormGrad_(other.NormGrad_),

            smallParamChange_(other.smallParamChange_),

            smallSSChange_(other.smallSSChange_),

            smallGradient_(other.smallGradient_)

    {}



LevenbergMarquardtParameters&
LevenbergMarquardtParameters::operator=(const LevenbergMarquardtParameters& other)
{
    if (this!=&other)
    {
        LevenbergMarquardtParameters tmp(other);
        swap(*this,tmp);
    }
    return *this;
}

void swap(LevenbergMarquardtParameters& one, LevenbergMarquardtParameters& other)
{
    std::swap(one.data_,other.data_);
    std::swap(one.w_,other.w_);

    std::swap(one.initialParam_,other.initialParam_);
    std::swap(one.nPar_,other.nPar_);
    std::swap(one.nData_,other.nData_);
    std::swap(one.dx_,other.dx_);
    std::swap(one.maxIter_,other.maxIter_);
    std::swap(one.maxFeval_,other.maxFeval_);
    std::swap(one.minParamChange_,other.minParamChange_);
    std::swap(one.minSSChange_,other.minSSChange_);
    std::swap(one.minGradient_,other.minGradient_);
    std::swap(one.maxLanda_,other.maxLanda_);
    std::swap(one.landa_,other.landa_);
    std::swap(one.landa0_,other.landa0_);
    std::swap(one.v_,other.v_);
    std::swap(one.nIter_,other.nIter_);
    std::swap(one.nFeval_,other.nFeval_);



    std::swap(one.currSS_,other.currSS_);
    std::swap(one.newSS0_,other.newSS0_);
    std::swap(one.currParam_,other.currParam_);
    std::swap(one.newParam_,other.newParam_);
    std::swap(one.newParam0_,other.newParam0_);
    std::swap(one.currYfit_,other.currYfit_);
    std::swap(one.newYfit_,other.newYfit_);
    std::swap(one.newYfit_,other.newYfit_);
std::swap(one.J_,other.J_);
    std::swap(one.G_,other.G_);
    std::swap(one.JTJ_,other.JTJ_);
    std::swap(one.JTJinv_,other.JTJinv_);

    std::swap(one.d_,other.d_);

    std::swap(one.optimParam_,other.optimParam_);

    std::swap(one.surpassIter_,other.surpassIter_);
    std::swap(one.surpassFeval_,other.surpassFeval_);

    std::swap(one.ParamChange_,other.ParamChange_);
    std::swap(one.SSChange_,other.SSChange_);
    std::swap(one.NormGrad_,other.NormGrad_);

    std::swap(one.smallParamChange_,other.smallParamChange_);

    std::swap(one.smallSSChange_,other.smallSSChange_);

    std::swap(one.smallGradient_,other.smallGradient_);

    }

LevenbergMarquardtParameters& LevenbergMarquardtParameters::optimize()
{
    std::cout<<"nIter"<<"\t"<<"currSS"<<"\t\t"<<"landa"<<"\t";
    std::cout<<"ParamChange"<<"\t"<<"SSChange"<<"\t"<<"NormGrad"<<"\n";
    initialize();
    while (!meetConvergenceCriteria())
            iterate();

    optimParam_=currParam_;
    return *this;
}


void LevenbergMarquardtParameters::iterate()
{
    computeJacobian();
    computeSearchDirection();
    updateLanda();
    std::cout<<nIter_<<"\t"<<currSS_<<"\t\t"<<landa_<<"\t";
    std::cout<<ParamChange_<<"\t"<<SSChange_<<"\t"<<NormGrad_<<"\n";

    nIter_++;
}


void LevenbergMarquardtParameters::computeJacobian()
{

    for (std::size_t i=0; i<nPar_; i++)
    {
        Parameters x(currParam_);
        x[i]+=dx_;
        std::vector<double> yfit=f_->yfit(x);
        //std::cout<<nFeval_++<<" ";
        for (std::size_t n=0;n<nData_;++n)
        {
            J_[n][i]=(yfit[n]-currYfit_[n])/dx_;
        }
    }
}


void LevenbergMarquardtParameters::computeSearchDirection()
{
    for (std::size_t i=0; i<nPar_; ++i)
        for (std::size_t j=0; j<nPar_; ++j)
        {
            JTJ_[i][j]=0;
            for (std::size_t n=0; n<nData_; ++n)
            {
                JTJ_[i][j]+=J_[n][i]*J_[n][j];
            }
        }



    for (std::size_t i=0; i<nPar_; ++i)
    {
        JTJ_[i][i]*=1+landa_;
    }

    JTJinv_=inv(JTJ_);



    for (std::size_t i=0; i<nPar_; ++i)
    {
        G_[i]=0;
        for (std::size_t n=0; n<nData_;++n)
        {
            G_[i]+=(data_[n]-currYfit_[n])*w_[n]*J_[n][i];
        }

    }
    for (std::size_t i=0; i<nPar_; ++i)
    {
        d_[i]=0;
        for (std::size_t j=0; j<nPar_;++j)
        {
            d_[i]+=JTJinv_[i][j]*G_[j];
        }

    }
    newParam_=currParam_;
    for (std::size_t i=0; i<nPar_; ++i)
        newParam_[i]+=d_[i];

    newYfit_=f_->yfit(newParam_);
    nFeval_++;

    newSS_=0;
    for (std::size_t n=0; n<nData_; n++)
    {
        newSS_+=(newYfit_[n]-data_[n])*(newYfit_[n]-data_[n])*w_[n];
    }
}


void LevenbergMarquardtParameters::updateLanda()
{
    std::size_t ifevalLoop=0;
    if (newSS_>=currSS_)
    {
        while(((newSS_>=currSS_)&&(nFeval_<maxFeval_))||(newSS_!=newSS_))
        {
            if (landa_*v_>=maxLanda_) break;
            landa0_=landa_;
            landa_=landa0_*v_;
            newSS0_=newSS_;
            newParam0_=newParam_;
            computeSearchDirection();
            ifevalLoop++;
            std::cerr<<landa_<<" ";
        }

    }
    else
    {
        landa0_=landa_;
        landa_=landa_/v_;
        newSS0_=newSS_;
        newParam0_=newParam_;
        newYfit0_=newYfit_;
        computeSearchDirection();
        ifevalLoop++;
        if (newSS_>=newSS0_)
        {
            landa_=landa0_;
            newParam_=newParam0_;
            newSS_=newSS0_;
            newYfit_=newYfit0_;
        }
    }
    if (newSS_>=currSS_)
    {
        ParamChange_=0;
        SSChange_=0;

    }
    else
    {
        ParamChange_=0;
        for (std::size_t i=0; i<nPar_; ++i)
            ParamChange_+=(currParam_[i]-newParam_[i])*(currParam_[i]-newParam_[i]);
        ParamChange_=sqrt(ParamChange_);
        SSChange_=currSS_-newSS_;
        currParam_=newParam_;
        currSS_=newSS_;
        currYfit_=newYfit_;
    }
    NormGrad_=0;
    for (std::size_t i=0; i<nPar_; ++i)
        NormGrad_+=G_[i]*G_[i];
    NormGrad_=sqrt(NormGrad_);
 }


bool LevenbergMarquardtParameters::meetConvergenceCriteria()
{
    surpassIter_=bool(nIter_>=maxIter_);
    surpassFeval_=nFeval_>=maxFeval_;
    surpassLanda_=landa_>=maxLanda_;
    smallParamChange_=ParamChange_<minParamChange_;
    smallSSChange_=SSChange_<minSSChange_;
    smallGradient_=NormGrad_<minGradient_;


    return surpassIter_||
            surpassFeval_||
            smallParamChange_||
            smallSSChange_||
            smallGradient_||
            surpassLanda_||
            currSS_!=currSS_;
}

void LevenbergMarquardtParameters::initialize()
{
    nIter_=0;
    nFeval_=0;
    ParamChange_=std::numeric_limits<double>::infinity();
    SSChange_=std::numeric_limits<double>::infinity();
    NormGrad_=std::numeric_limits<double>::infinity();
    currParam_=initialParam_;

    currYfit_=f_->yfit(currParam_);
    nFeval_++;
    currSS_=0;
    for (std::size_t n=0; n<nData_; ++n)
    {
        currSS_+=(currYfit_[n]-data_[n])*(currYfit_[n]-data_[n])*w_[n];
    }


}


Parameters LevenbergMarquardtParameters::OptimParameters()const
{
    return currParam_;
}

std::size_t LevenbergMarquardtParameters::numEval()const
{
    return nFeval_;
}
std::size_t LevenbergMarquardtParameters::numIter()const
{
    return nIter_;
}
double LevenbergMarquardtParameters::SS()const
{
    return currSS_;
}
std::vector<double> LevenbergMarquardtParameters::Gradient()const
{
    return G_;
}
