#include "Results.h"
#include <cmath>
#include <iostream>
#include <string>
#include <fstream>

const std::vector<Measurement>& Results::TNF()const
{
    return TNF_;
}

const std::vector<Measurement>& Results::IFN()const
{
    return IFN_;
}

const std::vector<Measurement>& Results::APC_expression()const
{
    return APC_expression_;
}

const std::vector<Measurement>& Results::NK_expression()const
{
    return NK_expression_;
}

const std::vector<Measurement>& Results::LT_expression()const
{
    return LT_expression_;
}


double Results::Duration() const
{
    return duration_;
}

Results::Results(std::string experimentName):
    TNF_(),IFN_(), APC_expression_(), NK_expression_ (), LT_expression_(), duration_()

{

    if (experimentName=="media")
    {
        duration_=120.0;
        TNF_.push_back(Measurement(16.0,2.1));
        TNF_.push_back(Measurement(48.0,1.0));
        TNF_.push_back(Measurement(119.0,0.5));
        IFN_.push_back(Measurement(16.0,0.02));
        IFN_.push_back(Measurement(48.0,0.03));
        IFN_.push_back(Measurement(119.0,0.01));
        APC_expression_.push_back(Measurement(0.0,3.2));
        APC_expression_.push_back(Measurement(16.0,3.4));
        APC_expression_.push_back(Measurement(119.0,2.8));
        NK_expression_.push_back (Measurement (24.0, 2.5));
        LT_expression_.push_back (Measurement (0.0,0.2));
        LT_expression_.push_back (Measurement (16.0,0.9));
        LT_expression_.push_back(Measurement (24.0, 0.8));
        LT_expression_.push_back(Measurement (119.0,2.1));
    }
    if (experimentName=="mtb")
    {
        duration_=120.0;
        TNF_.push_back(Measurement(16.0,50.5));
        TNF_.push_back(Measurement(48.0,40.3));
        TNF_.push_back(Measurement(119.0,30.2));
        IFN_.push_back(Measurement(16.0,8.1));
        IFN_.push_back(Measurement(48.0,12.0));
        IFN_.push_back(Measurement(119.0,28.3));
        APC_expression_.push_back(Measurement(0.0,3.2));
        APC_expression_.push_back(Measurement(16.0,16.4));
        APC_expression_.push_back(Measurement(119.0,2.3));
        NK_expression_.push_back (Measurement (24.0, 11.3));
        LT_expression_.push_back (Measurement (0.0,0.7));
        LT_expression_.push_back (Measurement (16.0,1.1));
        LT_expression_.push_back(Measurement (24.0, 1.4));
        LT_expression_.push_back(Measurement (119.0,9.1));
    }
    if (experimentName=="block")
    {
        duration_=120.0;
        TNF_.push_back(Measurement(16.0,63.7));
        TNF_.push_back(Measurement(48.0,58.6));
        TNF_.push_back(Measurement(119.0,52.4));
        IFN_.push_back(Measurement(16.0,10.3));
        IFN_.push_back(Measurement(48.0,8.1));
        IFN_.push_back(Measurement(119.0,12.4));
        /*  APC_expression_.push_back(Measurement(0.0,3.2));
        APC_expression_.push_back(Measurement(16.0,11.4));
        APC_expression_.push_back(Measurement(120.0,2.3));
        NK_expression_.push_back (Measurement (24.0, 11.3));
        LT_expression_.push_back (Measurement (0.0,1.2));
        LT_expression_.push_back (Measurement (16.0,1.1));
        LT_expression_.push_back(Measurement (24.0, 1.4));
        LT_expression_.push_back(Measurement (120.0,9.1));*/
    }
}




Results::Results(const std::vector<Measurement>& myTNF,
		 const std::vector<Measurement>& myIFN,
		 const std::vector<Measurement>& myAPCexpression,
		 const std::vector<Measurement>& myNKexpression,
		 const std::vector<Measurement>& myLTexpression,
		 double duration):
    TNF_(myTNF),
    IFN_(myIFN),
    APC_expression_(myAPCexpression),
    NK_expression_(myNKexpression),
    LT_expression_(myLTexpression),
    duration_(duration)
{if ((!TNF_.empty())&&duration_<TNF_[TNF_.size()-1].Time())
        duration_=TNF_[TNF_.size()-1].Time();


    if ((!IFN_.empty())&&duration_<IFN_[IFN_.size()-1].Time())
        duration_=IFN_[IFN_.size()-1].Time();

    if ((!APC_expression_.empty())&&duration_<APC_expression_[APC_expression_.size()-1].Time())
        duration_=APC_expression_[APC_expression_.size()-1].Time();

    if ((!NK_expression_.empty())&&duration_<NK_expression_[NK_expression_.size()-1].Time())
        duration_=NK_expression_[NK_expression_.size()-1].Time();

    if ((!LT_expression_.empty())&&duration_<LT_expression_[LT_expression_.size()-1].Time())
        duration_=LT_expression_[LT_expression_.size()-1].Time();

}

Results::Results(){}









std::vector<double> SumSquare_TNF(const Results one, const Results two)
{
    std::vector<double> ss_TNF;
    std::size_t n=0;
    std::size_t i=0;
    std::size_t j=0;

    for (;;)
    {
        if (one.TNF().empty())
            break;
        if (two.TNF().empty())
            break;

        if(two.TNF()[j].Time()==one.TNF()[i].Time())
        {
            ss_TNF.push_back(pow(two.TNF()[j].Measure()-one.TNF()[i].Measure(),2));
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

    return ss_TNF;
}

std::vector<double> SumSquare_IFN(const Results one, const Results two)
{
    std::vector<double> ss;
    std::size_t n=0;
    std::size_t i=0;
    std::size_t j=0;

    for (;;)
    {
        if (one.IFN().empty())
            break;
        if (two.IFN().empty())
            break;
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

    return ss;

}

std::vector<double> SumSquare_APCexpression(const Results one, const Results two)
{
    std::vector<double> ss;
    std::size_t n=0;
    std::size_t i=0;
    std::size_t j=0;

    for (;;)
    {
        if (one.APC_expression().empty())
            break;
        if (two.APC_expression().empty())
            break;

        if(two.APC_expression()[j].Time()==one.APC_expression()[i].Time())
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

    return ss;

}

std::vector<double> SumSquare_NKexpression(const Results one, const Results two)
{
    std::vector<double> ss;
    std::size_t n=0;
    std::size_t i=0;
    std::size_t j=0;

    for (;;)
    {
        if (one.NK_expression().empty())
            break;
        if (two.NK_expression().empty())
            break;

        if(two.NK_expression()[j].Time()==one.NK_expression()[i].Time())
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
    return ss;

}

std::vector<double> SumSquare_LTexpression(const Results one, const Results two)
{
    std::vector<double> ss;
    std::size_t n=0;
    std::size_t i=0;
    std::size_t j=0;

    for (;;)
    {
        if (one.LT_expression().empty())
            break;
        if (two.LT_expression().empty())
            break;

        if(two.LT_expression()[j].Time()==one.LT_expression()[i].Time())
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

    return ss;

}



std::vector<double> SumSquare_i(const Results& one, const Results& two)
{
    std::vector<double> ss=SumSquare_TNF(one,two);
    std::vector<double> ss_other=SumSquare_IFN(one,two);
    ss.insert(ss.end(),ss_other.begin(),ss_other.end());

    ss_other=SumSquare_APCexpression(one,two);
    ss.insert(ss.end(),ss_other.begin(),ss_other.end());

    ss_other=SumSquare_NKexpression(one,two);
    ss.insert(ss.end(),ss_other.begin(),ss_other.end());

    ss_other=SumSquare_LTexpression(one,two);
    ss.insert(ss.end(),ss_other.begin(),ss_other.end());

    return ss;
}

double SumSquare(const Results& one, const Results& two)
{
    std::vector<double> ss=SumSquare_i(one,two);
    double s=0;
    for (std::size_t i=0;i>ss.size();i++)
        s+=ss[i];
    return s;
}

void SumSquareTXT(const Results& one, const Results& two)
{
    std::string file="SumSquare.txt";
    std::ofstream fi;
    fi.open(file.c_str());
    fi<<"La suma de cuadrados obtenida para sus parámetros es ";
    fi<<SumSquare(one, two);

}


std::ostream& operator<<(std::ostream& s, const Results& res)
{
    s<<"Results \n\n";

    s<<"TNF\n";
    for (std::size_t i=0; i<res.TNF().size(); i++)
        s<<res.TNF()[i].Time()<<"\t"<<res.TNF()[i].Measure()<<"\n";

    s<<"IFN\n";
    for (std::size_t i=0; i<res.IFN().size(); i++)
        s<<res.IFN()[i].Time()<<"\t"<<res.IFN()[i].Measure()<<"\n";

    s<<"APC_expression \n";
    for (std::size_t i=0; i<res.APC_expression().size(); i++)
        s<<res.APC_expression()[i].Time()<<"\t"<<
           res.APC_expression()[i].Measure()<<"\n";

    s<<"NK_expression \n";
    for (std::size_t i=0; i<res.NK_expression().size(); i++)
        s<<res.NK_expression()[i].Time()<<"\t"<<
           res.NK_expression()[i].Measure()<<"\n";

    s<<"LT_expression \n";
    for (std::size_t i=0; i<res.LT_expression().size(); i++)
        s<<res.LT_expression()[i].Time()<<"\t"<<
           res.LT_expression()[i].Measure()<<"\n";

    return s;

}

std::vector<double> Results::getData()const
{
    std::vector<double> data;
    for (std::size_t i=0; i<TNF_.size(); ++i)
	data.push_back(TNF_[i].Measure());

    for (std::size_t i=0; i<IFN_.size(); ++i)
	data.push_back(IFN_[i].Measure());

    for (std::size_t i=0; i<APC_expression_.size(); ++i)
	data.push_back(APC_expression_[i].Measure());

    for (std::size_t i=0; i<NK_expression_.size(); ++i)
	data.push_back(NK_expression_[i].Measure());

    for (std::size_t i=0; i<LT_expression_.size(); ++i)
	data.push_back(LT_expression_[i].Measure());

     return data;
 }



Results::Results(const Results& other):
    TNF_(other.TNF_),
    IFN_(other.IFN_),
    APC_expression_(other.APC_expression_),
    NK_expression_(other.NK_expression_),
    LT_expression_(other.LT_expression_),
    duration_(other.duration_)
{}
Results& Results::operator=(const Results& other)
{
    if (this!=&other)
    {
	Results tmp(other);
	swap(*this,tmp);

    }
    return *this;
}

Results::~Results(){}
void swap(Results& one, Results& other)
{
    std::swap(one.TNF_,other.TNF_);
    std::swap(one.IFN_,other.IFN_);
    std::swap(one.APC_expression_,other.APC_expression_);
    std::swap(one.NK_expression_,other.NK_expression_);
    std::swap(one.LT_expression_,other.LT_expression_);
    std::swap(one.duration_,other.duration_);
}
