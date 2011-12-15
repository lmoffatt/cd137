#ifndef LEVENBERGMARQUARDT_H
#define LEVENBERGMARQUARDT_H
#include <vector>
#include "Parameters.h"

class ABC_function
{
public:
    virtual std::vector<double> yfit(const std::vector<double>& parameters)=0;
    virtual std::vector<double> yfit(const Parameters& parameters)=0;


};




class LevenbergMarquardt
{
public:

    std::vector<double> OptimParameters()const;
    std::vector< std::vector<double> > OptimParCov()const;

    std::size_t numEval()const;
    std::size_t numIter()const;
    double SS()const;
    std::vector<double> Gradient()const;




    LevenbergMarquardt& optimize();


    LevenbergMarquardt(ABC_function* f,
                       const std::vector<double>& data,
                       const std::vector<double>& initialParam);


    LevenbergMarquardt(const LevenbergMarquardt& other);

    friend void swap(LevenbergMarquardt& one, LevenbergMarquardt& other);

    LevenbergMarquardt& operator=(const LevenbergMarquardt& other);

    LevenbergMarquardt();

    ~LevenbergMarquardt(){}

    std::string report();

   // void reset(const SimParameters& sp,const Treatment& tr);

private:
    ABC_function* f_;
    std::vector<double> data_;

    std::vector<double> initialParam_;

    std::size_t nPar_;
    std::size_t nData_;



    // parameters of the optimization
    /// delta x used for Jacobian approximation
    double dx_;
    std::size_t maxIter_;
    std::size_t maxFeval_;

    double minParamChange_;
    double minSSChange_;
    double minGradient_;

    double maxLanda_;

    // variables that change on each iteration

    double landa_;
    double landa0_;
    double v_;

    std::size_t nIter_;
    std::size_t nFeval_;



    double currSS_;
    double newSS_;
    double newSS0_;
    std::vector<double> currParam_;
    std::vector<double> newParam_;
    std::vector<double> newParam0_;

    std::vector<double> currYfit_;
    std::vector<double> newYfit_;
    std::vector<double> newYfit0_;

    std::vector< std::vector< double> > J_;
    std::vector<double> G_;
    std::vector< std::vector<double> > JTJ_;
    std::vector< std::vector<double> > JTJinv_;

    std::vector<double> d_;


    std::vector<double> optimParam_;

    bool surpassIter_;
    bool surpassFeval_;
    bool surpassLanda_;

    double ParamChange_;
    double SSChange_;
    double NormGrad_;

    bool smallParamChange_;

    bool smallSSChange_;

    bool smallGradient_;


    bool meetConvergenceCriteria();

    void initialize();
    void iterate();

    void computeJacobian();
    void computeSearchDirection();
    void updateLanda();


    double SumSquare();




};


#endif // LEVENBERGMARQUARDT_H
