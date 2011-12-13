#ifndef PARAMETERS_H
#define PARAMETERS_H

#include<map>
#include <vector>

#include <string>
class Parameters {
public:
    double mean(const std::string& name)const;
    double pMean(const std::string& name)const;


    double pStd(const std::string& name)const;


    double residual(const std::string& name, double value)const;
    double pResidual(const std::string& name, double log10value)const;

    std::vector<double> residuals(const std::vector<std::string> names,const std::vector<double> log10Values);


    bool setMeans(const std::vector<std::string> names,const std::vector<double> values);
    bool setpMeans(const std::vector<std::string> names,const std::vector<double> log10values);
    bool setpStd(const std::vector<std::string> names,const std::vector<double> values);


    std::vector<double> means(const std::vector<std::string> names);

    std::vector<double> pMeans(const std::vector<std::string> names);

    void push_back(const std::string& name,double meanValue,double pStdValue);

    bool setpMean(const std::string& name, double value);

    bool setpStd(const std::string& name, double value);

    bool hasName(const std::string& name)const;


    Parameters randomSample()const;



private:
    std::map<std::string,double> mean_;
    std::map<std::string,double> pStd_;


};

double randNormal(double mean,double stddev);
double randNormal();


#endif // PARAMETERS_H
