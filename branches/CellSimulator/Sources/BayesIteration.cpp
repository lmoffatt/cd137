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

BayesIteration& BayesIteration::getPosterior(const Parameters& startingPoint)
{
    std::vector<double> data=getData();
    std::vector<double> w=getDataWeigth();

    Parameters p=priors_.back();

   std::size_t numIterations=200;

    LevenbergMarquardtParameters LM(this,
                                    data,
                                    startingPoint,
                                    w,
                                    numIterations);
    LM.optimize();
     std::ofstream f;

    f.open(filename_.c_str(),std::ios_base::app);

    f<<"--------------------------------------------------"
       "---------------------------------------------------\n";
    f<<"----------The prior is-----------------------------\n";
    f<<"--------------------------------------------------"
       "---------------------------------------------------\n";

    put(f,p);



    f<<"--------------------------------------------------"
       "---------------------------------------------------\n";
    f<<"----------Start from the following point------------------\n";
    f<<"--------------------------------------------------"
       "---------------------------------------------------\n";
    put(f,startingPoint);


    f<<"--------------------------------------------------"
       "---------------------------------------------------\n";
    f<<"----------Result of Levenberg Marquardt------------------\n";
    f<<"--------------------------------------------------"
       "---------------------------------------------------\n";
    f<<LM;
    put(f,LM.OptimParameters());


    f.close();
}

BayesIteration& BayesIteration::getPosterior(const Parameters& startingPoint,double factor, std::size_t numSeeds,double probParChange)
{
    std::vector<double> data=getData();
    std::vector<double> w=getDataWeigth();

    Parameters p=priors_.back();

    std::size_t numIterations=200;

    std::ofstream f;

    f.open(filename_.c_str(),std::ios_base::app);

    f<<"--------------------------------------------------"
       "---------------------------------------------------\n";
    f<<"----------The prior is-----------------------------\n";
    f<<"--------------------------------------------------"
       "---------------------------------------------------\n";

    put(f,p);

    for (std::size_t i=0;i<numSeeds;i++)
    {

        Parameters seed;
        seed=startingPoint.randomSample(p,factor,probParChange);
        LevenbergMarquardtParameters LM(this,
                                        data,
                                        seed,
                                        w,
                                        numIterations);
        LM.optimize();



        f<<"--------------------------------------------------"
           "---------------------------------------------------\n";
        f<<"----------Start from the following point------------------\n";
        f<<"--------------------------------------------------"
           "---------------------------------------------------\n";
        put(f,seed);


        f<<"--------------------------------------------------"
           "---------------------------------------------------\n";
        f<<"----------Result of Levenberg Marquardt------------------\n";
        f<<"--------------------------------------------------"
           "---------------------------------------------------\n";
        f<<LM;
        put(f,LM.OptimParameters());
    }

    f.close();
}



BayesIteration& BayesIteration::getPosterior()
{
    std::vector<double> data=getData();
    std::vector<double> w=getDataWeigth();
    std::vector<LevenbergMarquardtParameters> LMs;
    std::vector<Parameters> Ps;


    Parameters p=priors_.back();


    std::size_t factor=2;
    std::size_t numIterations=200;
    std::size_t numSeeds=2;
    std::map<double,Parameters> seeds=getRandomParameters(numSeeds*factor,0.2);

    std::map<double,Parameters> friuts;

    LevenbergMarquardtParameters LM(this,
                                    data,
                                    p,
                                    w,
                                    numIterations);
    LM.optimize();
    LMs.push_back(LM);
    Ps.push_back(LM.OptimParameters());

    friuts[LM.SS()]=LM.OptimParameters();
    std::ofstream f;

    f.open(filename_.c_str(),std::ios_base::app);

    f<<"--------------------------------------------------"
       "---------------------------------------------------\n";
    f<<"----------The prior is-----------------------------\n";
    f<<"--------------------------------------------------"
       "---------------------------------------------------\n";

    put(f,p);



    f<<"--------------------------------------------------"
       "---------------------------------------------------\n";
    f<<"----------Start from the center-------------------\n";
    f<<"--------------------------------------------------"
       "---------------------------------------------------\n";
    f<<LM;
    put(f,Ps.back());
    f.close();


    std::map<double,Parameters>::iterator it=seeds.begin();

    for (std::size_t i=0; i<numSeeds; i++)
    {
        Parameters initParam=(*it).second;
        double ss=(*it).first;
        ++it;
        LevenbergMarquardtParameters LM(this,
                                        data,
                                        initParam,
                                        w,
                                        numIterations);
        LM.optimize();
        friuts[LM.SS()]=LM.OptimParameters();

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
    double errorFactor=0.2;
    double errorShrinkFactor=3;

    std::size_t numCycle=3;
    for (size_t iCycle=0;
         iCycle<numCycle;
         iCycle++)
    {
        seeds=friuts;
        it=seeds.begin();
        size_t j=0;
        size_t jfactor=5;
        Parameters initParam=(*it).second;
        errorFactor/=errorShrinkFactor;



        for (std::size_t i=0; i<numSeeds; i++)
        {
            ++j;
            if(j==jfactor)
            {
                ++it;
                initParam=(*it).second;
                j=0;
            }
            initParam=initParam.randomSample(errorFactor);

            LevenbergMarquardtParameters LM(this,
                                            data,
                                            initParam,
                                            w,
                                            numIterations);
            LM.optimize();
            friuts[LM.SS()]=LM.OptimParameters();

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



 std::map<double,Parameters> BayesIteration::getRandomParameters(const Parameters& per,
                                                                 std::size_t num,
                                                                 double factor)
 {
     std::map<double,Parameters> myMap;
     for (std::size_t i=0; i<num; i++)
     {
         Parameters p=per.randomSample(factor);
         double ss=SumWeighedSquare(p);

         std::cout<<ss<<"\n";
     //    std::cout<<p;
         if (!(ss!=ss))
             myMap[ss]=p;
     }


     for (std::map<double,Parameters>::iterator it=myMap.begin();it!=myMap.end();++it)
     {
         double ss=(*it).first;
         std::cout<<(*it).first<<"\n";

     }


     return myMap;
 }

 std::map<double,Parameters> BayesIteration::getRandomParameters(std::size_t num,
                                                                 double factor)
 {
     std::map<double,Parameters> myMap;
     for (std::size_t i=0; i<num; i++)
     {
         Parameters p=priors_.back().randomSample(factor);
         double ss=SumWeighedSquare(p);

         std::cout<<ss<<"\n";
     //    std::cout<<p;
         if (!(ss!=ss))
             myMap[ss]=p;
     }


     for (std::map<double,Parameters>::iterator it=myMap.begin();it!=myMap.end();++it)
     {
         double ss=(*it).first;
         std::cout<<(*it).first<<"\n";

     }


     return myMap;
 }

