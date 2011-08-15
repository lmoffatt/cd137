#include <cmath>
#include<limits>
#include "Includes/LevenbergMarquardt.h"
#include "Includes/MatrixInverse.h"

LevenbergMarquardt::LevenbergMarquardt(
    ABC_function* f,
    const std::vector<double>& data,
    const std::vector<double>& initialParam):
    f_(f),
    data_(data),
    initialParam_(initialParam),
    nPar_(initialParam_.size()),
    nData_(data.size()),
    dx_(1e-7),
    maxIter_(100),
    maxFeval_(100),
    minParamChange_(1e-5),
    minSSChange_(1e-4),
    minGradient_(1e-4),
    landa_(1e5),
    nIter_(0),
    nFeval_(0),
    currSS_(std::numeric_limits<double>::infinity()),
    newSS_(std::numeric_limits<double>::infinity()),
    currParam_(std::vector<double>(nPar_)),
    newParam_(currParam_),
    currYfit_(std::vector<double>(nData_)),
    newYfit_(currYfit_),
    J_(std::vector< std::vector< double> >(nPar_,std::vector<double>(nData_))),
    G_(currParam_),
    JTJ_(std::vector< std::vector<double> > (nPar_,std::vector<double>(nPar_))),
    JTJinv_(JTJ_),
    d_(currParam_),
    optimParam_(currParam_),
    surpassIter_(false),
    surpassFeval_(false),
    ParamChange_(std::numeric_limits<double>::infinity()),
    SSChange_(ParamChange_),
    NormGrad_(ParamChange_),
    smallParamChange_(false),
    smallSSChange_(false),
    smallGradient_(false){}



LevenbergMarquardt::LevenbergMarquardt(){}

LevenbergMarquardt::LevenbergMarquardt (const LevenbergMarquardt& other):
    data_(other.data_),
            initialParam_(other.initialParam_),
            nPar_(other.nPar_),
            nData_(other.nData_),
            dx_(other.dx_),
            maxIter_(other.maxIter_),
            maxFeval_(other.maxFeval_),
            minParamChange_(other.minParamChange_),
            minSSChange_(other.minSSChange_),
            minGradient_(other.minGradient_),
            landa_(other.landa_),
            nIter_(other.nIter_),
            nFeval_(other.nFeval_),



            currSS_(other.currSS_),
            newSS_ (other.newSS_),
            currParam_ (other.currParam_),
            newParam_ (other.newParam_),
            currYfit_(other.currYfit_),
            newYfit_(other.newYfit_),

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



LevenbergMarquardt&
LevenbergMarquardt::operator=(const LevenbergMarquardt& other)
{
    if (this!=&other)
    {
        LevenbergMarquardt tmp(other);
        swap(*this,tmp);
    }
    return *this;
}

void swap(LevenbergMarquardt& one, LevenbergMarquardt& other)
{
    std::swap(one.data_,other.data_);
    std::swap(one.initialParam_,other.initialParam_);
    std::swap(one.nPar_,other.nPar_);
    std::swap(one.nData_,other.nData_);
    std::swap(one.dx_,other.dx_);
    std::swap(one.maxIter_,other.maxIter_);
    std::swap(one.maxFeval_,other.maxFeval_);
    std::swap(one.minParamChange_,other.minParamChange_);
    std::swap(one.minSSChange_,other.minSSChange_);
    std::swap(one.minGradient_,other.minGradient_);
    std::swap(one.landa_,other.landa_);
    std::swap(one.nIter_,other.nIter_);
    std::swap(one.nFeval_,other.nFeval_);



    std::swap(one.currSS_,other.currSS_);
    std::swap(one.newSS_,other.newSS_);
    std::swap(one.currParam_,other.currParam_);
    std::swap(one.newParam_,other.newParam_);
    std::swap(one.currYfit_,other.currYfit_);
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

LevenbergMarquardt& LevenbergMarquardt::optimize()
{
    initialize();
    while (!meetConvergenceCriteria())
            iterate();

    optimParam_=currParam_;
    return *this;
}


void LevenbergMarquardt::iterate()
{
    computeJacobian();
    computeSearchDirection();
    updateLanda();
    std::cout<<nIter_<<"\t"<<currSS_<<"\n";
    nIter_++;
}


void LevenbergMarquardt::computeJacobian()
{

    for (std::size_t i=0; i<currParam_.size(); i++)
    {
        std::vector<double> x(currParam_);
        x[i]+=dx_;
	std::vector<double> yfit=f_->yfit(x);
	std::cout<<nFeval_++<<" ";
        for (std::size_t j=0;j<currYfit_.size();++j)
        {
            J_[i][j]=(yfit[j]-currYfit_[j])/dx_;
        }
    }
}


void LevenbergMarquardt::computeSearchDirection()
{
    for (std::size_t i=0; i<nPar_; ++i)
        for (std::size_t j=0; j<nPar_; ++j)
        {
            JTJ_[i][j]=0;
            for (std::size_t n=0; n<nData_; ++n)
            {
                JTJ_[i][j]+=J_[i][n]*J_[j][n];
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
            G_[i]+=J_[i][n]*(data_[n]-currYfit_[i]);
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
    for (std::size_t n=0; n<nData_; ++n)
    {
        newSS_+=(newYfit_[n]-data_[n])*(newYfit_[n]-data_[n]);
    }
}


void LevenbergMarquardt::updateLanda()
{

    while((newSS_>=currSS_)&&(nFeval_<maxFeval_))
    {
        landa_=landa_*10;
        computeSearchDirection();
    }
    ParamChange_=0;
    for (std::size_t i=0; i<nPar_; ++i)
	ParamChange_+=(currParam_[i]-newParam_[i])*(currParam_[i]-newParam_[i]);
    ParamChange_=sqrt(ParamChange_);
    NormGrad_=0;
    for (std::size_t i=0; i<nPar_; ++i)
	NormGrad_+=G_[i]*G_[i];
    NormGrad_=sqrt(NormGrad_);

    SSChange_=currSS_-newSS_;

    currParam_=newParam_;
    currYfit_=newYfit_;
    currSS_=newSS_;
}

bool LevenbergMarquardt::meetConvergenceCriteria()
{
    surpassIter_=bool(nIter_>=maxIter_);
    surpassFeval_=nFeval_>=maxFeval_;

    smallParamChange_=ParamChange_<minParamChange_;
    smallSSChange_=SSChange_<minSSChange_;
    smallGradient_=NormGrad_<minGradient_;

    return surpassIter_||
	    surpassFeval_||
	    smallParamChange_||
	    smallSSChange_||
	    smallGradient_;
}

void LevenbergMarquardt::initialize()
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
	currSS_+=(currYfit_[n]-data_[n])*(currYfit_[n]-data_[n]);
    }


}
