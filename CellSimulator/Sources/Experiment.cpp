#include "Includes/Experiment.h"
#include <fstream>
#include <iostream>
void Experiment::push_back(const Treatment& treatment,
               const Results& result){
    treatments_.push_back(treatment);
    results_.push_back(result);
}

const Treatment& Experiment::Treatment_i(std::size_t i)const
{
    return treatments_[i];
}

const Results& Experiment::Result_i(std::size_t i)const{
    return results_[i];
}


std::size_t Experiment::size()const
{
    return results_.size();
}

std::vector<double> SumSquare_i(const Experiment& one, const Experiment& two)
{
    std::vector<double> ss=SumSquare_i(one.Result_i(0),two.Result_i(0));
    for (std::size_t i=1;i<one.size();i++)
    {
        std::vector<double> ss2=SumSquare_i(one.Result_i(i),two.Result_i(i));
        ss.insert(ss.end(),ss2.begin(),ss2.end());
    }

    return ss;
}


double SumSquare(const Experiment& one, const Experiment& two)
{
    std::vector<double> ss=SumSquare_i(one,two);
     double s=0;
    for (std::size_t i=0;i<ss.size();i++)
        s+=ss[i];
    return s;
}

void SumSquareTXT(const Experiment& one, const Experiment& two)
{
    std::string file="SumSquare.txt";
    std::ofstream f;
    f.open(file.c_str());
    f<<"La suma de cuadrados obtenida para sus par�metros es ";
    f<<SumSquare(one, two);

}
