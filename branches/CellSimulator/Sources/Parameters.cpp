#include "Parameters.h"
#include <limits>
#include <cmath>
#include <cstdlib>
#include <istream>
#include <sstream>
#include <iostream>

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
    std::swap(one.mode_,other.mode_);

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
        return pStd_[(*it).second]*10;
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
    double stds=log10(maxValue_p68/minValue_p34)*10/2.0;
    push_back_dB(name,sqrt(minValue_p34*maxValue_p68),log10(maxValue_p68/minValue_p34)*20.0);
}


void Parameters::push_back_2S(const std::string& name,double minValue_p02,double maxValue_p98)
{
    push_back_dB(name,sqrt(minValue_p02*maxValue_p98),log10(maxValue_p98/minValue_p02)*10.0/4.0);

}

void Parameters::push_back_3S(const std::string& name,double minValue_p001,double maxValue_p999)
{
    push_back_dB(name,sqrt(minValue_p001*maxValue_p999),log10(maxValue_p999/minValue_p001)*10.0/6.0);

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


std::string Parameters::indexToName(std::size_t i)const
{
    for (std::map<std::string,std::size_t>::const_iterator it=name_.begin();
         it!=name_.end();
         ++it)
    {
        if (i==it->second)
            return it->first;
    }
         return "";

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
    Parameters sample(*this);
    for (std::size_t i=0; i<pMean_.size();i++)

    {
        sample.pMean_[i]=randNormal(pMean_[i],pStd_[i]);
        sample.pStd_[i]=0;
    }
    Parameters s(sample);
    return s;

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
Parameters Parameters::randomSample(Parameters prior,double factor)const{
    Parameters sample;
    for (std::map<std::string,std::size_t>::const_iterator it=name_.begin();
         it!=name_.end();
         ++it)

    {
        double m=pow(10,randNormal(pMean(it->first),factor*prior.pStd(it->first)));
        sample.push_back_dB(it->first,m,0);

    }
    return sample;

}




Parameters Parameters::randomSample(Parameters prior,double factor,double probIncludeParameter)const{
    Parameters sample,sampleout;
    for (std::size_t i=0; i<pMean_.size();i++)
    {
        std::map<std::string,std::size_t>::const_iterator it;
    for (it=name_.begin();
         it!=name_.end();
         ++it)

    {
        if (it->second==i)
            break;
    }

        double r=(1.0*rand())/(1.0*RAND_MAX);

        if (r<probIncludeParameter)
        {
            std::string str=it->first;

            double m=pMean(it->first);
            double s=factor*prior.pStd(it->first);
            m=randNormal(m,s);
            m=pow(10,m);

            sample.push_back_dB(str,m,0);
        }
        else
        {
            std::string str=it->first;

            double m=pMean(it->first);
             m=pow(10,m);

            sample.push_back_dB(str,m,0);
        }
    }

     sampleout=sample;
    return sampleout;

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

    //Box-Muller method http://en.wikipedia.org/wiki/Normal_distribution#Generating_values_from_normal_distribution
    double U=(1.0*rand())/RAND_MAX;
    double V=(1.0*rand())/RAND_MAX;
    const double PI=3.1415926;

    double r=sqrt(-2*log(U))*cos(2*PI*V);
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

void Parameters::setCovariance(const std::vector< std::vector <double> >& cov)
{
    if (cov.size()==size())
    {
        cov_=cov;
        for (std::size_t i=0; i<cov.size(); i++)
        {
            pStd_[i]=sqrt(cov_[i][i]);
        }
    }
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


std::ostream& operator<<(std::ostream& s, const Parameters& p){

    s<<"Parameters\n";
    s<<"Begin\n";

    for (std::size_t i=0; i<p.size(); i++)
    {
        s<<p.indexToName(i)<<"\t"<<pow(10,p.pMean_[i])<<"\t"<<p.pStd_[i]*10<<"\t dB\n";
    }
    s<<"End\n";
    return s;
   }

 std::istream& operator>>(std::istream& s, Parameters& p)
{
     std::string line;
     std::getline(s,line);
     if (line.compare("Parameters")==0)
     {
         std::getline(s,line);
         std::getline(s,line);
         while (!line.empty()&& line.compare("End")!=0)
         {
            std::string name;
            double value1,value2;
            std::string errorType;
            std::stringstream ss(line);
            ss>>name;
            ss>>value1;
            ss>>value2;
            ss>>errorType;
            if (errorType.compare("dB")==0)
                p.push_back_dB(name,value1,value2);
            else if (errorType.compare("1S")==0)
                p.push_back_1S(name,value1,value2);
            else if (errorType.compare("2S")==0)
                p.push_back_2S(name,value1,value2);
            else if (errorType.compare("3S")==0)
                p.push_back_3S(name,value1,value2);
            std::getline(s,line);
            if (line.compare("End")==0)
                std::cout<<"here";

         }

     }
      return s;
}


 Parameters::Parameters():
     name_(std::map<std::string, std::size_t> ()),
     pMean_(std::vector<double>()),
     pStd_(std::vector<double> ()),  // not in dB
     cov_(std::vector< std::vector <double> > ()),
       mode_("")

 {}


 bool areTheSame(const Parameters& one, const Parameters& other)
 {
     if (one.size()!=other.size())
         return false;
     for (std::size_t i=0; i<one.size();i++)
         if (one.indexToName(i)!=other.indexToName(i))
             return false;
     return true;
 }


 double dbDistance(const Parameters& one,const Parameters& other)
 {
     double result=0;

     for (std::size_t i=0;i<one.size();i++)
        {
            result+=pow(other.pMean(other.indexToName(i))-one.pMean(one.indexToName(i)),2);
        }
     result=sqrt(result/one.size())*10;
     return result;
 }


 double Parameters::chi2Distance(const Parameters &other)const
 {
     double result=0;

     for (std::size_t i=0;i<size();i++)
        {
            result+=pow((other.pMean(other.indexToName(i))-pMean(indexToName(i))/pStd(indexToName(i))),2);
        }
    // result=sqrt(result/one.size())*10;
     return result;

 }
