#include "BayesIteration.h"

 std::vector<double> ABC_data::getDataWeigth()
{
    std::vector<double> w=getDataStandardError();
    for (std::size_t i=0; i<w.size(); i++)
    {
        w[i]=1.0/w[i]/w[i];
    }
    return w;
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


    BayesIteration::BayesIteration(ABC_model* m,
                   Parameters prior,
                   ABC_data* d):
        m_(m),
        priors_(1,prior),
        data_(1,d),
        posterior_()
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
    std::vector<double> BayesIteration::yfit(const Parameters& parameters)
    {

        std::vector<double> simulatedResults=m_->yfit(parameters);
         std::vector<double> param=parameters.pMeans();

       simulatedResults.insert(simulatedResults.end(),param.begin(),param.end());
        return simulatedResults;
    }


    std::vector<double> BayesIteration::getData()
    {
        std::vector<double>data =data_.back()->getData();
        std::vector<double> param=priors_.back().pMeans();

      data.insert(data.end(),param.begin(),param.end());
       return data;


    }

    std::vector<double> BayesIteration::getDataStandardError()
    {

        std::vector<double>se =data_.back()->getDataStandardError();
        std::vector<double> param=priors_.back().pStds();

       se.insert(se.end(),param.begin(),param.end());
       return se;

    }




    void BayesIteration::getPosterior()
    {
        std::vector<double> data=getData();
        std::vector<double> w=getDataWeigth();
        std::vector<LevenbergMarquardtParameters> LMs;


        for (std::size_t i=0; i<numSeeds_; i++)
        {
            Parameters initParam=priors_.back().randomSample();

            LevenbergMarquardtParameters LM(this,
                                            data,
                                            initParam,
                                            w);
            LM.OptimParameters();
            LMs.push_back(LM);
        }








        }


