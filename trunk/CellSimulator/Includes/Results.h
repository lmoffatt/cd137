#ifndef RESULTS_H
#define RESULTS_H
#include <vector>
#include <string>
#include <iostream>

#include "Measurement.h"

class Results
{
public:
    const std::vector<Measurement>& TNF()const;
    const std::vector<Measurement>& IFN() const;
    const std::vector<Measurement>& APC_expression() const;
    const std::vector<Measurement>& NK_expression() const;
    const std::vector<Measurement>& LT_expression () const;

    Results();
    Results(std::string experimentName);
    Results(std::vector<Measurement> myTNF,
            std::vector<Measurement> myIFN,
            std::vector<Measurement> myAPCexpression,
            std::vector<Measurement> myNKexpression,
            std::vector<Measurement> myLTexpression,
            double duration);



    double Duration()const;
 //   std::string filename;





private:
    std::vector<Measurement> TNF_;
    std::vector<Measurement> IFN_;
    std::vector<Measurement> APC_expression_;
    std::vector<Measurement> NK_expression_;
    std::vector<Measurement> LT_expression_;
    double duration_;


};

double SumSquare(const Results& one, const Results& two);
std::vector<double> SumSquare_i(const Results& one, const Results& two);
void SumSquareTXT(const Results& one, const Results& two);


std::ostream& operator<<(std::ostream& s, const Results& res);

#endif // RESULTS_H
