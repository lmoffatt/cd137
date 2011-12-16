#include "Parameters.h"
#include <limits>
#include <cmath>
#include <cstdlib>



Parameters::Parameters(const Parameters& other):
    name_(other.name_),
    pMean_(other.pMean_),
    pStd_(other.pStd_),
    cov_(other.cov_){}



Parameters& Parameters::operator=(const Parameters& other)
{
    if (this!=&other)
    {
        Parameters tmp(other);
        swap(*this,tmp);
    }
    return *this;
}

void swap(Parameters& one, Parameters& other)
{
    std::swap(one.name_,other.name_);
    std::swap(one.pMean_,other.pMean_);
    std::swap(one.pStd_,other.pStd_);
    std::swap(one.cov_,other.cov_);

}


double Parameters::mean(const std::string& name)const
{
    std::map<std::string,std::size_t>::const_iterator it=name_.find(name);
    if(it!=name_.end())
        return pow(10,pMean_[(*it).second]);
    else
        return std::numeric_limits<double>::quiet_NaN();
}

double Parameters::mean_ratio(const std::string& name)const
{
    std::map<std::string,std::size_t>::const_iterator it=name_.find(name);
    if(it!=name_.end())
        return mean(name)/(mean(name)+1.0);
    else
        return std::numeric_limits<double>::quiet_NaN();
}



double Parameters::pMean(const std::string& name)const{
    return std::log10(mean(name));
}

/// returns the standard deviation
double Parameters::pStd(const std::string& name)const
{
    std::map<std::string,size_t>::const_iterator it=name_.find(name);
    if((it!=name_.end())&& !pStd_.empty())
        return pStd_[(*it).second];
    else
        return std::numeric_limits<double>::quiet_NaN();

}



/// returns the standard deviation in dB (deciBel)
double Parameters::dBStd(const std::string& name)const
{
    std::map<std::string,size_t>::const_iterator it=name_.find(name);
    if((it!=name_.end())&& !pStd_.empty())
        return pStd_[(*it).second]*20;
    else
        return std::numeric_limits<double>::quiet_NaN();

}



double Parameters::pResidual(const std::string& name, double log10value)const{
    return 20.0*(log10value-pMean(name))/pStd(name);
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

std::size_t Parameters::nameIndex(const std::string& name)const
{
    return (*name_.find(name)).second;

}



bool Parameters::setMeans(const std::vector<std::string> names,const std::vector<double> values)
{
    for (std::size_t i=0;i<names.size();i++)
    {
        if (hasName(names[i]))
            pMean_[nameIndex(names[i])]=log10(values[i]);
        else return false;
    }
    return true;

}
bool Parameters::setpMeans(const std::vector<std::string> names,const std::vector<double> log10values)
{
    for (std::size_t i=0;i<names.size();i++)
    {
        if (hasName(names[i]))
            pMean_[nameIndex(names[i])]=log10values[i];
        else return false;
    }
    return true;

}


bool Parameters::setpStd(const std::vector<std::string> names,const std::vector<double> values)
{
    for (std::size_t i=0;i<names.size();i++)
    {
        if (hasName(names[i]))
            pStd_[nameIndex(names[i])]=values[i];
        else return false;
    }
    return true;

}


std::vector<double> Parameters::means(const std::vector<std::string> names)const
{
    std::vector<double> result(names.size());
    for (std::size_t i=0;i<names.size();i++)
    {
        result[i]=mean(names[i]);
    }
    return result;


}

std::vector<double> Parameters::pMeans(const std::vector<std::string> names)const
{
    std::vector<double> result(names.size());
    for (std::size_t i=0;i<names.size();i++)
    {
        result[i]=pMean(names[i]);
    }
    return result;
}


void Parameters::push_back_dB(const std::string& name,double meanValue,double pStdValue)
{
    name_[name]=pMean_.size();
    pMean_.push_back(log10(meanValue));
    pStd_.push_back(pStdValue/10);

}

void Parameters::push_back_1S(const std::string& name,double minValue_p34,double maxValue_p68)
{
    push_back_dB(name,sqrt(minValue_p34*maxValue_p68),log10(maxValue_p68/minValue_p34)*20.0);
}


void Parameters::push_back_2S(const std::string& name,double minValue_p02,double maxValue_p98)
{
    push_back_dB(name,sqrt(minValue_p02*maxValue_p98),log10(maxValue_p98/minValue_p02)*10.0);

}

void Parameters::push_back_3S(const std::string& name,double minValue_p001,double maxValue_p999)
{
    push_back_dB(name,sqrt(minValue_p001*maxValue_p999),log10(maxValue_p999/minValue_p001)*20.0/3.0);

}

std::string Parameters::mode()const
{
    return mode_;
}
Parameters&
Parameters::setMode(const std::string& modeString)
{
    mode_=modeString;
    return *this;
}





bool Parameters::setpMean(const std::string& name, double value)
{
    std::map<std::string,std::size_t>::iterator it=name_.find(name);
    if(it!=name_.end())
    {
        pMean_[(*it).second]=value;
        return true;
    }
    else
        return false;

}

bool Parameters::setpStd(const std::string& name, double value)
{
    std::map<std::string,std::size_t>::iterator it=name_.find(name);
    if(it!=name_.end())
    {
        pStd_[(*it).second]=value;
        return true;
    }
    else
        return false;

}

bool Parameters::hasName(const std::string& name)const{
    return name_.find(name)!=name_.end();
}


std::vector<std::string> Parameters::names()const
{
    std::vector<std::string>  myNames;
    for (std::map<std::string,std::size_t>::const_iterator it=name_.begin();
         it!=name_.end();
         ++it)
    {
        myNames.push_back(it->first);

    }
    return myNames;

}

std::vector<std::string> Parameters::commonNames(const Parameters& other)const
{
    std::vector<std::string>  myCommonNames;
    for (std::map<std::string,std::size_t>::const_iterator it=name_.begin();
         it!=name_.end();
         ++it)
    {
        if (other.hasName(it->first))
        myCommonNames.push_back(it->first);

    }
    return myCommonNames;

}



Parameters Parameters::randomSample()const
{
    Parameters sample;
    for (std::map<std::string,std::size_t>::const_iterator it=name_.begin();
         it!=name_.end();
         ++it)

    {
        double m=pow(10,randNormal(pMean(it->first),pStd(it->first)));
        sample.push_back_dB(it->first,m,0);

    }
    return sample;

}

Parameters Parameters::randomSample(double factor)const{
    Parameters sample;
    for (std::map<std::string,std::size_t>::const_iterator it=name_.begin();
         it!=name_.end();
         ++it)

    {
        double m=pow(10,randNormal(pMean(it->first),factor*pStd(it->first)));
        sample.push_back_dB(it->first,m,0);

    }
    return sample;

}


std::size_t Parameters::size()const
{
    return name_.size();
}


std::vector<double> Parameters::residuals(const Parameters& prior)const
{
    std::vector<std::string> cnames=commonNames(prior);
    std::vector<double> result(cnames.size());
    for (std::size_t i=0;i<cnames.size();i++)
    {
        result[i]=prior.pResidual(cnames[i],pMean(cnames[i]));
    }
    return result;
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

Parameters& Parameters::applyParameters(const Parameters& other)
{
    std::vector<std::string> cnames=commonNames(other);
    setpMeans(cnames,other.pMeans(cnames));
    return *this;
}

Parameters& Parameters::scaleError(double factor)
{
    for (std::size_t i=0; i<pStd_.size(); i++)
    {
        pStd_[i]*=factor;
    }
    return *this;
}



std::vector<double> Parameters::pMeans()const
{
    return pMean_;
}

std::vector<double> Parameters::pStds()const
{
    return pStd_;
}


const double& Parameters::operator[](std::size_t i)const
{
    return pMean_[i];
}
double& Parameters::operator[](std::size_t i){
return pMean_[i];
}
