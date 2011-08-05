#include <cmath>
#include "Results.h"

const std::vector<Measurement>& Results::TNF()const
{
    return TNF_;
}

Results::Results(std::string experimentName):
    TNF_(),duration_()

{

    if (experimentName=="Fig_1_2")
    {
        duration_=120.0;
        TNF_.push_back(Measurement(16.0,50.0));
        TNF_.push_back(Measurement(48.0,40.0));
        TNF_.push_back(Measurement(120.0,30.0));
    }
}

Results::Results(std::vector<Measurement> myTNF):
    TNF_(myTNF),
    duration_(0)
{
    if (duration_<TNF_[TNF_.size()-1].Time())
        duration_=TNF_[TNF_.size()-1].Time();

}

const std::vector<Measurement>& Results::IFN()const
{
    return IFN_;
}

Results::Results(std::string experimentName):
    IFN_(),duration_()

{

    if (experimentName=="Fig_1_2")
    {
        duration_=120.0;
        IFN_.push_back(Measurement(16.0,12.0));
        IFN_.push_back(Measurement(48.0,30.0));
        IFN_.push_back(Measurement(120.0,50.0));
    }
}

Results::Results(std::vector<Measurement> myIFN):
    IFN_(myIFN),
    duration_(0)
{
    if (duration_<IFN_[IFN_.size()-1].Time())
        duration_=IFN_[IFN_.size()-1].Time();

}

std::vector<double> SumSquare_TNF(const Results one, const Results two)
{
    std::vector<double> ss_TNF=0;
    std::size_t n=0;
    std::size_t i=0;
    std::size_t j=0;

    for (;;)
    {
        if(two.TNF()[j].Time()==one.TNF()[i].Time())
        {
            ss_TNF.push_back(std::pow(two.TNF()[j].Measure()-one.TNF()[i].Measure(),2));
            n++;
        }
        if (one.TNF()[i].Time()<=two.TNF()[j].Time())
        {
            i++;
            if (i>=one.TNF().size())
                break;
        }
        else
        {
            j++;
            if (j>=two.TNF().size())
                break;
        }

    }

    if (n>0)
        return ss_TNF;
    else
        return -1;
}

std::vector<double> SumSquare_IFN(const Results one, const Results two)
{
    std::vector<double> ss=0;
    std::size_t n=0;
    std::size_t i=0;
    std::size_t j=0;

    for (;;)
    {
        if(two.IFN()[j].Time()==one.IFN()[i].Time())
        {
            ss.push_back(std::pow(two.IFN()[j].Measure()-one.IFN()[i].Measure(),2));
            n++;
        }
        if (one.IFN()[i].Time()<=two.IFN()[j].Time())
        {
            i++;
            if (i>=one.IFN().size())
                break;
        }
        else
        {
            j++;
            if (j>=two.IFN().size())
                break;
        }

    }

    if (n>0)
        return ss_IFN;
    else
        return -1;
};

std::vector<double> SumSquare_APCexpression(const Results one, const Results two)
{
    std::vector<double> ss=0;
    std::size_t n=0;
    std::size_t i=0;
    std::size_t j=0;

    for (;;)
    {
        if(two.APC_expression()[j]==one.APC_expression()[i].Time())
        {
            ss.push_back(std::pow(two.APC_expression()[j].Measure()-one.APC_expression()[i].Measure(),2));
            n++;
        }
        if (one.APC_expression()[i].Time()<=two.APC_expression()[j].Time())
        {
            i++;
            if (i>=one.APC_expression().size())
                break;
        }
        else
        {
            j++;
            if (j>=two.APC_expression().size())
                break;
        }

    }

    if (n>0)
        return ss_APC_expression;
    else
        return -1;
};

std::vector<double> SumSquare_NKexpression(const Results one, const Results two)
{
    std::vector<double> ss=0;
    std::size_t n=0;
    std::size_t i=0;
    std::size_t j=0;

    for (;;)
    {
        if(two.NK_expression()[j]==one.NK_expression()[i].Time())
        {
            ss.push_back(std::pow(two.NK_expression()[j].Measure()-one.NK_expression()[i].Measure(),2));
            n++;
        }
        if (one.NK_expression()[i].Time()<=two.NK_expression()[j].Time())
        {
            i++;
            if (i>=one.NK_expression().size())
                break;
        }
        else
        {
            j++;
            if (j>=two.NK_expression().size())
                break;
        }

    }

    if (n>0)
        return ss_NK_expression();
    else
        return -1;
};

std::vector<double> SumSquare_LTexpression(const Results one, const Results two)
{
    std::vector<double> ss=0;
    std::size_t n=0;
    std::size_t i=0;
    std::size_t j=0;

    for (;;)
    {
        if(two.LT_expression()[j]==one.LT_expression()[i].Time())
        {
            ss.push_back(std::pow(two.LT_expression()[j].Measure()-one.LT_expression()[i].Measure(),2));
            n++;
        }
        if (one.LT_expression()[i].Time()<=two.LT_expression()[j].Time())
        {
            i++;
            if (i>=one.LT_expression().size())
                break;
        }
        else
        {
            j++;
            if (j>=two.LT_expression().size())
                break;
        }

    }

    if (n>0)
        return ss_LT_expression();
    else
        return -1;
};



std::vector<double> SumSquare_i(const Results one, const Results two)
{
    std::vector<double> ss=SumSquare_TNF(one,two);
    std::vector<double> ss_other=SumSquare_IFN(one,two);
    ss.insert(ss.end(),ss_other.begin(),ss_other.end());

    ss_other=SumSquare_APCexpression(IFN(one,two);
    ss.insert(ss.end(),ss_other.begin(),ss_other.end());

    std::vector<double> ss_APC=SumSquare_APCexpression(one,two);
    std::vector<double> ss_NK=SumSquare_NKexpression(one,two);
    std::vector<double> ss_LT=SumSquare_LTexpression(one,two);
    double s=0;
    return ss;
}

double SumSquare(const Results one, const Results two)
{
    std::vector<double> ss=SumSquare_i(one,two);
     double s=0;
    for (std::size_t i=0;i>ss.size();i++)
        s+=ss[i];
    return s;
}
