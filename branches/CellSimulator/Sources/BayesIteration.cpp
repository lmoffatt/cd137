#include <vector>
#include <fstream>
#include <cmath>
#include "BayesIteration.h"

std::vector<double> ABC_data::getDataWeigth()
{
    std::vector<double> w=getDataStandardError();
    for (std::size_t i=0; i<w.size(); i++)
    {
        w[i]=1.0/w[i]/w[i];
    }
    std::vector<double> wW(w);
    return wW;
}

ABC_data::~ABC_data(){
}

ABC_model::~ABC_model(){}

BayesIteration::~BayesIteration(){
}



Parameters BayesIteration::Posterior()const
{
    return posterior_;
}

Parameters BayesIteration::Prior(std::size_t n)const
{
    return priors_[n];

}


BayesIteration::BayesIteration(const ABC_model* m,
                               Parameters prior,
                               const ABC_data* d,
                               const std::string& filename):
    m_(m),
    priors_(1,prior),
    data_(1,d),
    posterior_(),
    numSeeds_(20),
    filename_(filename)
{
    //getPosterior();
}


BayesIteration& BayesIteration::addNewData(ABC_data* d)
{
    data_.push_back(d);
    priors_.push_back(posterior_);
    getPosterior();
}



BayesIteration::BayesIteration(const BayesIteration& other):
    m_(other.m_),
    data_(other.data_),
    priors_(other.priors_),
    posterior_(other.posterior_)
{}


/*   void swap(BayesIteration& one, BayesIteration& other);

    BayesIteration& operator=(const BayesIteration& other);

    BayesIteration();

    ~BayesIteration(){}
*/
// void reset(const SimParameters& sp,const Treatment& tr);



ABC_model& BayesIteration::setData(const ABC_data& experimentalData)
{}


BayesIteration::BayesIteration(){}




std::vector<double> BayesIteration::yfit (const std::vector<double>& param){}
std::vector<double> BayesIteration::yfit(const Parameters& parameters)const
{

    std::vector<double> simulatedResults=m_->yfit(parameters);
    std::vector<double> param=parameters.pMeans();

    simulatedResults.insert(simulatedResults.end(),param.begin(),param.end());
    return simulatedResults;
}


std::vector<double> BayesIteration::getData()const
{
    std::vector<double>data =data_.back()->getData();
    std::vector<double> param=priors_.back().pMeans();

    data.insert(data.end(),param.begin(),param.end());
    std::vector<double> d(data);
    return d;


}

std::vector<double> BayesIteration::getDataStandardError()const
{

    std::vector<double>se =data_.back()->getDataStandardError();
    std::vector<double> param=priors_.back().pStds();

    se.insert(se.end(),param.begin(),param.end());
    std::vector<double> o(se);
    return o;

}




BayesIteration& BayesIteration::getPosterior()
{
    std::vector<double> data=getData();
    std::vector<double> w=getDataWeigth();
    std::vector<LevenbergMarquardtParameters> LMs;
    std::vector<Parameters> Ps;





    Parameters p=priors_.back();

    p.scaleError(0.2);

    std::size_t factor=100;
    std::map<double,Parameters> seeds=getRandomParameters(numSeeds_*factor);


    LevenbergMarquardtParameters LM(this,
                                    data,
                                    p,
                                    w);
    LM.optimize();
    LMs.push_back(LM);
    Ps.push_back(LM.OptimParameters());
    std::ofstream f;

    f.open(filename_.c_str(),std::ios_base::app);
    f<<"--------------------------------------------------"
       "---------------------------------------------------\n";
    f<<"----------Start from the center-------------------\n";
    f<<"--------------------------------------------------"
       "---------------------------------------------------\n";
    f<<LM;
    put(f,Ps.back());
    f.close();


    std::map<double,Parameters>::iterator it=seeds.begin();

    for (std::size_t i=0; i<numSeeds_; i++)
    {
        Parameters initParam=(*it).second;
        ++it;
        LevenbergMarquardtParameters LM(this,
                                        data,
                                        initParam,
                                        w);
        LM.optimize();
        LMs.push_back(LM);
        Ps.push_back(LM.OptimParameters());
        f.open(filename_.c_str(),std::ios_base::app);
        f<<"--------------------------------------------------"
           "---------------------------------------------------\n";
        f<<"----------Start from the perisphery-------------------\n";
        f<<"--------------------------------------------------"
           "---------------------------------------------------\n";
        f<<LM;
        put(f,Ps.back());
        f.close();

    }


    return *this;
}


std::ostream& BayesIteration::put(std::ostream& s,const Parameters& parameters)const
{
    m_->put(s,parameters);
    return s;
}


 void BayesIteration::setFilename(const std::string filename)
{
    filename_=filename;
}

 double BayesIteration::SumWeighedSquare(const Parameters& p)
 {
     std::vector<double> d=this->getData();
     std::vector<double> w=this->getDataWeigth();

     std::vector<double> y=this->yfit(p);
     double SSW=0;
     for (std::size_t i=0;i<d.size();i++)
     {
         SSW+=pow(y[i]-d[i],2)*w[i];
     };
     return SSW;
 }


 std::map<double,Parameters> BayesIteration::getRandomParameters(std::size_t num)
 {
     std::map<double,Parameters> myMap;
     for (std::size_t i=0; i<num; i++)
     {
         Parameters p=priors_.back().randomSample(0.001);
         double ss=SumWeighedSquare(p);

         std::cout<<ss<<"\n";
     //    std::cout<<p;
         myMap[ss]=p;
     }
     return myMap;
 }
