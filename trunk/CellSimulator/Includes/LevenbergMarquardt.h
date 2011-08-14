#ifndef LEVENBERGMARQUARDT_H
#define LEVENBERGMARQUARDT_H
#include <vector>


class LevenbergMarquardt
{
public:

    LevenbergMarquardt& optimize();


    LevenbergMarquardt(std::vector<double> (*fun) (std::vector<double>),
                       const std::vector<double>& data,
                       const std::vector<double>& initialParam);

private:
    std::vector<double> (*fun_) (std::vector<double> );
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

    // variables that change on each iteration

    double landa_;


    std::size_t nIter_;
    std::size_t nFeval_;



    double currSS_;
    double newSS_;
    std::vector<double> currParam_;
    std::vector<double> newParam_;

    std::vector<double> currYfit_;
    std::vector<double> newYfit_;

    std::vector< std::vector< double> > J_;
    std::vector<double> G_;
    std::vector< std::vector<double> > JTJ_;
    std::vector< std::vector<double> > JTJinv_;

    std::vector<double> d_;


    std::vector<double> optimParam_;

    bool surpassIter_;
    bool surpassFeval_;

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
