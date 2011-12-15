#ifndef EXPERIMENT_H
#define EXPERIMENT_H
#include <cstddef>
#include "Includes/Results.h"
#include "Includes/Treatment.h"


class Experiment: public ABC_data
{
public:
    const Treatment& Treatment_i(std::size_t i)const;

    const Results& Result_i(std::size_t i)const;

    std::size_t size()const;

    void push_back(const Treatment& treatment,
                   const Results& result);

    std::vector<double> getData()const;

    std::vector<double> getDataStandardError()const;

    Experiment(const Experiment& other);
    ~Experiment();
    Experiment();
    Experiment& operator=(const Experiment& other);
    friend void swap(Experiment& one, Experiment& two);

    friend std::ostream& operator<<(std::ostream& s, const Experiment& e);


private:
    std::vector<Treatment> treatments_;
    std::vector<Results> results_;


};

double SumSquare(const Experiment& one, const Experiment& two);
void SumSquareTXT(const Experiment& one, const Experiment& two);
std::vector<double>SumSquare_i(const Experiment& one, const Experiment& two);


#endif // EXPERIMENT_H
