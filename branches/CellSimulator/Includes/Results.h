#ifndef RESULTS_H
#define RESULTS_H
#include <vector>
#include <string>
#include <iostream>

#include "Measurement.h"

#include "Includes/BayesIteration.h"

class Results:public ABC_data
{
public:
    const std::vector<Measurement>& TNF()const;
    const std::vector<Measurement>& IFN() const;
    const std::vector<Measurement>& APC_expression() const;
    const std::vector<Measurement>& NK_expression() const;
    const std::vector<Measurement>& LT_expression () const;
    const std::vector<Measurement>& APC_IFNg()const;
    const std::vector<Measurement>& APC_TNFa()const;
    const std::vector<Measurement>& NK_IFNg()const;
    const std::vector<Measurement>& NK_TNFa()const;
    const std::vector<Measurement>& LT_IFNg()const;
    const std::vector<Measurement>& LT_TNFa()const;
    const std::vector<Measurement>& LT_Apoptosis()const;
    const std::vector<Measurement>& Proliferation()const;
    const std::vector<Measurement>& num_cells() const;

    Results();
    Results(std::string experimentName);
    Results(const std::vector<Measurement>& myTNF,
	    const std::vector<Measurement>& myIFN,
	    const std::vector<Measurement>& myAPCexpression,
	    const std::vector<Measurement>& myNKexpression,
	    const std::vector<Measurement>& myLTexpression,
            const std::vector<Measurement>& myAPC_IFNg,
            const std::vector<Measurement>& myAPC_TNFa,
            const std::vector<Measurement>& myNK_IFNg,
            const std::vector<Measurement>& myNK_TNFa,
            const std::vector<Measurement>& myLT_IFNg,
            const std::vector<Measurement>& myLT_TNFa,
            const std::vector<Measurement>& myLT_Apoptosis,
            const std::vector<Measurement>& myProliferation,
            const std::vector<Measurement>& mynum_cells,
            double duration);

    virtual std::vector<double> getData()const;
    virtual std::vector<double> getDataStandardError()const;




    double Duration()const;
 //   std::string filename;


    Results(const Results& other);
    Results& operator=(const Results& other);
    ~Results();
    friend void swap(Results& one, Results& other);



private:
    std::vector<Measurement> TNF_;
    std::vector<Measurement> IFN_;
    std::vector<Measurement> APC_expression_;
    std::vector<Measurement> NK_expression_;
    std::vector<Measurement> LT_expression_;
    std::vector<Measurement> APC_IFNg_;
    std::vector<Measurement> APC_TNFa_;
    std::vector<Measurement> NK_IFNg_;
    std::vector<Measurement> NK_TNFa_;
    std::vector<Measurement> LT_IFNg_;
    std::vector<Measurement> LT_TNFa_;
    std::vector<Measurement> LT_Apoptosis_;
    std::vector<Measurement> Proliferation_;
    std::vector<Measurement> num_cells_;
    double duration_;


};

double SumSquare(const Results& one, const Results& two);
std::vector<double> SumSquare_i(const Results& one, const Results& two);
void SumSquareTXT(const Results& one, const Results& two);


std::ostream& operator<<(std::ostream& s, const Results& res);

#endif // RESULTS_H
