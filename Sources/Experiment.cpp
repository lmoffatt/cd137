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
    f<<"La suma de cuadrados obtenida para sus parámetros es ";
    f<<SumSquare(one, two);

}

std::vector<double> Experiment::getData()const
{
    std::vector<double> data;
    for (std::size_t i=0; i<results_.size(); i++)
    {
        std::vector<double> datai=Result_i(i).getData();
        for (std::size_t j=0;j<datai.size();++j)
            data.push_back(datai[j]);
     }
    return data;
}

std::vector<double> Experiment::getDataStandardError()const
{
    std::vector<double> data;
    for (std::size_t i=0; i<results_.size(); i++)
    {
        std::vector<double> datai=Result_i(i).getDataStandardError();
        for (std::size_t j=0;j<datai.size();++j)
            data.push_back(datai[j]);
     }
    return data;
}



Experiment::Experiment(const Experiment& other):
    treatments_(other.treatments_),
    results_(other.results_)
{}
Experiment::~Experiment(){}
Experiment::Experiment(){}
Experiment& Experiment::operator=(const Experiment& other)
{
    if (this!=&other)
    {
	Experiment tmp(other);
	swap(*this, tmp);
    }
    return *this;
}
 void swap(Experiment& one, Experiment& two)
{
    std::swap(one.treatments_,two.treatments_);
    std::swap(one.results_,two.results_);

}

 std::ostream& operator<<(std::ostream& s, const Experiment& e)
 {
     s<<"\n Experiment\n";
     for (std::size_t i=0; i<e.size(); ++i)
     {
	 s<<"\n Treatment number "<<i<<"\n";
	 s<<e.Treatment_i(i);

	 s<<"\n Results \n";
	 s<<e.Result_i(i);
     }

     return s;

 }
