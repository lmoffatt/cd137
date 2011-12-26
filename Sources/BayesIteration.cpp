#include <vector>
#include <fstream>
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
                               const ABC_data* d):
    m_(m),
    priors_(1,prior),
    data_(1,d),
    posterior_(),
    numSeeds_(20)
{
    getPosterior();
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

    p.scaleError(0.01);

    LevenbergMarquardtParameters LM(this,
                                    data,
                                    p,
                                    w);
    LM.optimize();
    LMs.push_back(LM);
    Ps.push_back(LM.OptimParameters());
    std::string filename="parametersOptimal.txt";
    std::ofstream f;
    f.open(filename.c_str(),std::ios_base::app);
    f<<LM;
    put(f,Ps.back());
    f.close();



    for (std::size_t i=0; i<numSeeds_; i++)
    {
        Parameters initParam=p.randomSample();

        LevenbergMarquardtParameters LM(this,
                                        data,
                                        initParam,
                                        w);
        LM.optimize();
        LMs.push_back(LM);
        Ps.push_back(LM.OptimParameters());
        f.open(filename.c_str(),std::ios_base::app);
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
