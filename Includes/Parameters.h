#ifndef PARAMETERS_H
#define PARAMETERS_H

#include<map>
#include <vector>

#include <string>
class Parameters
{
public:
    double mean(const std::string& name)const;

    /// pMean("someParameter")=log10(mean("someParameter")
    double pMean(const std::string& name)const;


    double pStd(const std::string& name)const;

    double dBStd(const std::string& name)const;


    double residual(const std::string& name, double value)const;
    double pResidual(const std::string& name, double log10value)const;

    std::vector<double> residuals(const std::vector<std::string> names,const std::vector<double> log10Values);

    std::vector<double> residuals(const Parameters& sample)const;



    bool setMeans(const std::vector<std::string> names,const std::vector<double> values);
    bool setpMeans(const std::vector<std::string> names,const std::vector<double> log10values);
    bool setpStd(const std::vector<std::string> names,const std::vector<double> values);


    std::vector<double> means(const std::vector<std::string> names)const;

    std::vector<double> pMeans()const;

    std::vector<double> pStds()const;



    std::vector<double> pMeans(const std::vector<std::string> names)const;



    void push_back(const std::string& name,double meanValue,double pStdValue);

    void push_back_1S(const std::string& name,double minValue_p34,double maxValue_p68);

    void push_back_2S(const std::string& name,double minValue_p02,double maxValue_p98);

    void push_back_3S(const std::string& name,double minValue_p001,double maxValue_p999);


    bool setpMean(const std::string& name, double value);

    bool setpStd(const std::string& name, double value);

    bool hasName(const std::string& name)const;

    std::vector<std::string> names()const;

    std::size_t nameIndex(const std::string& name)const;

    std::vector<std::string> commonNames(const Parameters& other)const;

    std::size_t size()const;


    ///returns the log of the mean
    const double& operator[](std::size_t)const;
    double& operator[](std::size_t);


    //applies the commonNames mean values to *this
    Parameters& applyParameters(const Parameters& other);


    Parameters randomSample()const;

    Parameters(const Parameters& other);
    Parameters(){}
    ~Parameters(){}

    Parameters& operator=(const Parameters& other);

    void friend swap(Parameters& one, Parameters& other);



private:
    std::map<std::string, std::size_t> name_;
    std::vector<double> pMean_;
    std::vector<double> pStd_;
    std::vector< std::vector <double> > cov_;


};

double randNormal(double mean,double stddev);
double randNormal();


#endif // PARAMETERS_H
