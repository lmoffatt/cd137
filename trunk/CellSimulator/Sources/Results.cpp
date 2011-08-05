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

    if (experimentName=="Fig_1_2")
    {
        duration_=120.0;
        TNF_.push_back(Measurement(16.0,50.0));
        TNF_.push_back(Measurement(48.0,40.0));
        TNF_.push_back(Measurement(120.0,30.0));
        IFN_.push_back(Measurement(16.0,8.0));
        IFN_.push_back(Measurement(48.0,20));
        IFN_.push_back(Measurement(120.0,50));
        APC_expression_.push_back(Measurement(0.0,3.2));
        APC_expression_.push_back(Measurement(16.0,11.4));
        APC_expression_.push_back(Measurement(120.0,2.3));
        NK_expression_.push_back (Measurement (24.0, 11.3));
        LT_expression_.push_back (Measurement (0.0,1.2));
        LT_expression_.push_back (Measurement (16.0,1.1));
        LT_expression_.push_back(Measurement (24.0, 1.4));
        LT_expression_.push_back(Measurement (120.0,9.1));
    }
}

Results::Results(std::vector<Measurement> myTNF,
                 std::vector<Measurement> myIFN,
                 std::vector<Measurement> myAPCexpression_,
                 std::vector<Measurement> myNKexpression_,
                 std::vector<Measurement> myLTexpression_,
                 double duration_):
                 TNF_(myTNF),
                 IFN_(myIFN),
                 APC_expression_(myAPCexpression_),
                 NK_expression_(myNKexpression_),
                 LT_expression_(myLTexpression_),
    duration_(0)
{if (duration_<TNF_[TNF_.size()-1].Time())
        duration_=TNF_[TNF_.size()-1].Time();


 if (duration_<IFN_[IFN_.size()-1].Time())
            duration_=IFN_[IFN_.size()-1].Time();

 if (duration_<APC_expression_[APC_expression_.size()-1].Time())
 duration_=APC_expression_[APC_expression_.size()-1].Time();

 if (duration_<NK_expression_[NK_expression_.size()-1].Time())
 duration_=NK_expression_[NK_expression_.size()-1].Time();

 if (duration_<LT_expression_[LT_expression_.size()-1].Time())
 duration_=LT_expression_[LT_expression_.size()-1].Time();

 };

 Results::Results(){}









 std::vector<double> SumSquare_TNF(const Results one, const Results two)
 {
     std::vector<double> ss_TNF;
     std::size_t n=0;
     std::size_t i=0;
     std::size_t j=0;

     for (;;)
     {
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



 std::vector<double> SumSquare_i(const Results one, const Results two)
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

 double SumSquare(const Results one, const Results two)
 {
     std::vector<double> ss=SumSquare_i(one,two);
      double s=0;
     for (std::size_t i=0;i>ss.size();i++)
         s+=ss[i];
     return s;
 }

 void SumSquareTXT(const Results one, const Results two)
 {
     std::string file="SumSquare.txt";
     std::ofstream fi;
     fi.open(file.c_str());
     fi<<"La suma de cuadrados obtenida para sus parámetros es ";
     fi<<SumSquare(one, two);

 }

