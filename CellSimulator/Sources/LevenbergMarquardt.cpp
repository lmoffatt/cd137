#include <cmath>
#include<limits>
#include "Includes/LevenbergMarquardt.h"
#include "Includes/MatrixInverse.h"

LevenbergMarquardt::LevenbergMarquardt(
    std::vector<double> (*fun) (std::vector<double>),
    const std::vector<double>& data,
    const std::vector<double>& initialParam):
    fun_(fun),
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
    currSS_(std::numeric_limits<double>::quiet_NaN()),
    newSS_(std::numeric_limits<double>::quiet_NaN()),
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
    nIter_++;
}


void LevenbergMarquardt::computeJacobian()
{

    for (std::size_t i=0; i<currParam_.size(); i++)
    {
        std::vector<double> x(currParam_);
        x[i]+=dx_;
        std::vector<double> yfit=(*fun_)(x);
        nFeval_++;
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

    newYfit_=(*fun_)(newParam_);
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

    currYfit_=(*fun_)(currParam_);
    nFeval_++;
    currSS_=0;
    for (std::size_t n=0; n<nData_; ++n)
    {
	currSS_+=(currYfit_[n]-data_[n])*(currYfit_[n]-data_[n]);
    }

}
