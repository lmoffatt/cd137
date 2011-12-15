#ifndef BAYESITERATION_H
#define BAYESITERATION_H
#include <vector>
#include "LevenbergMarquardtParameters.h"


class ABC_data
{
public:
    virtual std::vector<double> getData()const=0;
    virtual std::vector<double> getDataStandardError()const=0;
    virtual std::vector<double> getDataWeigth();
    virtual ~ABC_data();
};


class ABC_model:public ABC_function
{
public:
    virtual ABC_model& setData(const ABC_data& experimentalData)=0;
    virtual ~ABC_model();

};




class BayesIteration:public ABC_model, public ABC_data
{
public:

    Parameters Posterior()const;

    Parameters Prior(std::size_t n=0)const;


    BayesIteration(ABC_model* f,
                   Parameters prior,
                   ABC_data* d);


    BayesIteration& addNewData(ABC_data* d);

    virtual ~BayesIteration();

    virtual std::vector<double> yfit (const std::vector<double>& param);
    virtual std::vector<double> yfit(const Parameters& parameters);

    virtual std::vector<double> getData();
    virtual std::vector<double> getDataStandardError();

    virtual ABC_model& setData(const ABC_data& experimentalData);

    BayesIteration(const BayesIteration& other);
    /*

    friend void swap(BayesIteration& one, BayesIteration& other);

    BayesIteration& operator=(const BayesIteration& other);
*/
    BayesIteration();



private:
    void getPosterior();

    ABC_model* m_;

    std::vector<ABC_data*> data_;

    std::vector<Parameters> priors_;

    Parameters posterior_;
    std::size_t numSeeds_;

    std::vector<LevenbergMarquardtParameters> LM_;


    };

#endif // BAYESITERATION_H
