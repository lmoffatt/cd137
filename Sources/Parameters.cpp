#include "Parameters.h"
#include <limits>
#include <cmath>
#include <cstdlib>



double Parameters::mean(const std::string& name)const
{
    std::map<std::string,double>::const_iterator it=mean_.find(name);
    if(it!=mean_.end())
        return (*it).second;
    else
        return std::numeric_limits<double>::quiet_NaN();



}
double Parameters::pMean(const std::string& name)const{
    return std::log10(mean(name));
}


double Parameters::pStd(const std::string& name)const
{
    std::map<std::string,double>::const_iterator it=pStd_.find(name);
    if(it!=mean_.end())
        return (*it).second;
    else
        return std::numeric_limits<double>::quiet_NaN();

}


double Parameters::pResidual(const std::string& name, double log10value)const{
    return (log10value-pMean(name))/pStd(name);
}

std::vector<double> Parameters::residuals(const std::vector<std::string> names,const std::vector<double> log10Values)
{
    std::vector<double> result(names.size());
    for (std::size_t i=0;i<names.size();i++)
    {
        result[i]=pResidual(names[i],log10Values[i]);
    }
    return result;

}


bool Parameters::setMeans(const std::vector<std::string> names,const std::vector<double> values)
{
    for (std::size_t i=0;i<names.size();i++)
    {
        if (hasName(names[i]))
            mean_[names[i] ]=values[i];
        else return false;
    }
    return true;

}
bool Parameters::setpMeans(const std::vector<std::string> names,const std::vector<double> log10values)
{
    for (std::size_t i=0;i<names.size();i++)
    {
        if (hasName(names[i]))
            mean_[names[i]]=pow(10,log10values[i]);
        else return false;
    }
    return true;

}


bool Parameters::setpStd(const std::vector<std::string> names,const std::vector<double> values)
{
    for (std::size_t i=0;i<names.size();i++)
    {
        if (hasName(names[i]))
            pStd_[names[i]]=values[i];
        else return false;
    }
    return true;

}


std::vector<double> Parameters::means(const std::vector<std::string> names)
{
    std::vector<double> result(names.size());
    for (std::size_t i=0;i<names.size();i++)
    {
        result[i]=mean(names[i]);
    }
    return result;


}

std::vector<double> Parameters::pMeans(const std::vector<std::string> names)
{
    std::vector<double> result(names.size());
    for (std::size_t i=0;i<names.size();i++)
    {
        result[i]=pMean(names[i]);
    }
    return result;
}


void Parameters::push_back(const std::string& name,double meanValue,double pStdValue)
{
    mean_[name]=meanValue;
    pStd_[name]=pStdValue;
}

bool Parameters::setpMean(const std::string& name, double value)
{
    std::map<std::string,double>::iterator it=mean_.find(name);
    if(it!=mean_.end())
    {
        (*it).second=pow(10,value);
        return true;
    }
    else
        return false;

}

bool Parameters::setpStd(const std::string& name, double value)
{
    std::map<std::string,double>::iterator it=pStd_.find(name);
    if(it!=pStd_.end())
    {
        (*it).second=value;
        return true;
    }
    else
        return false;

}

bool Parameters::hasName(const std::string& name)const{
    return mean_.find(name)!=mean_.end();
}


Parameters Parameters::randomSample()const
{
    Parameters sample;
    for (std::map<std::string,double>::const_iterator it=mean_.begin();
         it!=mean_.end();
         ++it)
    {
        double m=pow(10,randNormal(pMean(it->first),pStd(it->first)));
        sample.push_back(it->first,m,0);

    }
    return sample;




}

double randNormal(double mean,double stddev)
{
    return randNormal()*stddev+mean;
}
double randNormal()
{
    const std::size_t n=20;
    double r=0;
    for (std::size_t i=0;i<n;++i)
        r+=(1.0*rand())/RAND_MAX;

    r=(r-0.5*n)/sqrt(n/12);
    return r;

}



