#include "Includes/LevenbergMarquardt.h"
#include "Includes/MatrixInverse.h"

LevenbergMarquardt::LevenbergMarquardt(
    std::vector<double> (*fun) (std::vector<double>),
    const std::vector<double>& data,
    const std::vector<double>& initialParam):
    fun_(fun),
    data_(data),
    initialParam_(initialParam){}




LevenbergMarquardt& LevenbergMarquardt::optimize()
{
    initialize();
    while (!meetConvergenceCriteria())
            iterate();

    return *this;
}


void LevenbergMarquardt::iterate()
{
    computeJacobian();
    computeSearchDirection();
    updateLanda();
}


void LevenbergMarquardt::computeJacobian()
{
    currYfit_=(*fun_)(currParam_);
    nFeval_++;

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
        JTJ_[i][i]*=1+landa_[i];
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
        for (std::size_t j=0; n<nPar_;++n)
        {
            d_[i]+=JTJinv_[i][j]*G_[j];
        }

    }
    newParam_(currParam_);
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


}

bool LevenbergMarquardt::meetConvergenceCriteria()
{
    bool smallGradient=true;
    for (std::size_t i=0; i<nPar_; ++i)




}
