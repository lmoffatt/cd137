#include "Results.h"
#include <cmath>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>

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

const std::vector<Measurement>& Results::APC_IFNg()const
{
    return APC_IFNg_;

}

const std::vector<Measurement>& Results::APC_TNFa()const
{
    return APC_TNFa_;

}

const std::vector<Measurement>& Results::NK_IFNg()const
{
    return NK_IFNg_;

}

const std::vector<Measurement>& Results::NK_TNFa()const
{
    return NK_TNFa_;

}

const std::vector<Measurement>& Results::LT_IFNg()const
{
    return LT_IFNg_;

}

const std::vector<Measurement>& Results::LT_TNFa()const
{
    return LT_TNFa_;

}

const std::vector<Measurement>& Results::LT_Apoptosis()const
{
    return LT_Apoptosis_;

}

const std::vector<Measurement>& Results::Proliferation()const
{
    return Proliferation_;

}

const std::vector<Measurement>& Results::num_cells() const
{
    return num_cells_;
}

double Results::Duration() const
{
    return duration_;
}

Results::Results(std::string experimentName):
    TNF_(),IFN_(), APC_expression_(), NK_expression_ (), LT_expression_(), APC_IFNg_(), APC_TNFa_(),
    NK_IFNg_(), NK_TNFa_(), LT_IFNg_(), LT_TNFa_(), LT_Apoptosis_(),Proliferation_(), num_cells_(),duration_()

{

    if (experimentName=="media")
    {
        duration_=120.0;
/*1*/        TNF_.push_back(Measurement(16.0,2.402,0.722));//
/*2*/        TNF_.push_back(Measurement(48.0,0.738,0.446));//
/*3*/        TNF_.push_back(Measurement(119.0,0.176,0.113));//
/*4*/        IFN_.push_back(Measurement(16.0,0.034,0.034));//
/*5*/        IFN_.push_back(Measurement(48.0,0.219,0.149));//
/*6*/        IFN_.push_back(Measurement(119.0,0.146,0.094));//
/*7*/        APC_expression_.push_back(Measurement(0.0,3.0,0.6));//
/*8*/        APC_expression_.push_back(Measurement(16.0,3.5,0.34));//
/*9*/        APC_expression_.push_back(Measurement(119.0,4.0,1.5));//
             NK_expression_.push_back(Measurement(0.0,1.54,0.67));
/*10*/       NK_expression_.push_back (Measurement (24.0,1.68,0.70));//

/*11*/        LT_expression_.push_back (Measurement (0.0,2.0,1.1));//
/*12*/        LT_expression_.push_back (Measurement (16.0,1.6,0.7));//
/*13*/        LT_expression_.push_back(Measurement (24.0, 3.1, 1.8));//
/*14*/        LT_expression_.push_back(Measurement (119.0,6.3, 2.6));//
/*15*/        APC_IFNg_.push_back(Measurement (16.0,2.83,1.16));//
/*16*/        APC_IFNg_.push_back (Measurement (119.0,2.73,2.05));//
/*17*/        APC_TNFa_.push_back(Measurement (16.0,5.93,2.73));//
/*18*/        APC_TNFa_.push_back (Measurement (119.0,4.12,2.91));//
/*19*/        NK_IFNg_.push_back(Measurement (24.0,3.55,1.43));//
/*20*/        NK_TNFa_.push_back (Measurement (24.0,3.06,1.44));//
/*21*/        LT_IFNg_.push_back (Measurement (119.0,3.4,1.1));//
/*22*/        LT_TNFa_.push_back(Measurement (119.0,1.54,1.33));//
/*23*/        LT_Apoptosis_.push_back(Measurement(119.0,16.55,4.10));//
/*24*/        Proliferation_.push_back(Measurement(119,3516.0,935.0));//
/*25*/        num_cells_.push_back(Measurement(24.4,1.0e6,5.0e5));//
/*26*/        num_cells_.push_back(Measurement(118.0,1.0e6,5.0e5));//

    }
    if (experimentName=="mtb")
    {
        duration_=120.0;
          /*27*/            TNF_.push_back(Measurement(16.0,50.29,8.82));//
          /*28*/            TNF_.push_back(Measurement(48.0,42.00,7.737));//
           /*29*/           TNF_.push_back(Measurement(119.0,28.309,5.613));//
                           IFN_.push_back(Measurement(16.0,4.794,0.89));//
   /*30*/                  IFN_.push_back(Measurement(48.0,12.29,1.873));//
     /*31*/                 IFN_.push_back(Measurement(119.0,28.209,5.584));//
   /*32*/                APC_expression_.push_back(Measurement(0.0,3.0,0.6));//  //:means true error
        /*33*/        APC_expression_.push_back(Measurement(16.0,14.5,2.0));
        /*34*/        APC_expression_.push_back(Measurement(119.0,6.9,3.6));
        /*35*/        NK_expression_.push_back (Measurement (24.0, 11.32,1.36));
                       LT_expression_.push_back(Measurement (0.0, 2.1, 1.9));
      /*36*/        LT_expression_.push_back(Measurement (24.0, 4.1, 2.3));//
         /*37*/             LT_expression_.push_back(Measurement (119.0,29.7, 7.77));
        /*39*/              APC_IFNg_.push_back(Measurement (16.0,7.66,2.20));//
       /*40*/             APC_TNFa_.push_back (Measurement (119.0,7.05,4.03));
        /*41*/             NK_IFNg_.push_back(Measurement (24.0,27.18,2.32));
        /*42*/             NK_TNFa_.push_back (Measurement (24.0,5.81,0.97));
                         LT_TNFa_.push_back(Measurement (119.0,4.27,0.59));
                         LT_IFNg_.push_back(Measurement (119.0,9.4,1.5));
                         LT_Apoptosis_.push_back(Measurement(119.0,27.61,2.57));
        /*48*/           Proliferation_.push_back(Measurement(119.0,14173.0,1240.0));
        /*49*/       num_cells_.push_back(Measurement(24.0,1.0e6,0.5e6));//
        /*50*/       num_cells_.push_back(Measurement(119.0,1.0e6,0.5e6));//


    }
    if (experimentName=="block")
    {
        duration_=120.0;

              TNF_.push_back(Measurement(16.0,61.503,8.527));
            TNF_.push_back(Measurement(48.0,54.45,8.102));
/*53*/        TNF_.push_back(Measurement(119.0,38.86,6.632));
//
/*54*/        IFN_.push_back(Measurement(16.0,8.2280,1.251));//
//
              IFN_.push_back(Measurement(48.0,7.911,1.208));
//
/*56*/        IFN_.push_back(Measurement(119.0,13.091,2.24));//
/*57*/        APC_IFNg_.push_back(Measurement (16.0,14.11,4.02));//
/*58*/
/*59*/        APC_TNFa_.push_back(Measurement (16.0,41.7,6.99));//
//
/*60*/        APC_TNFa_.push_back (Measurement (119.0,12.03,4.13));//
/*61*/
/*62*/        NK_TNFa_.push_back (Measurement (24.0,8.66,1.48));//
              NK_IFNg_.push_back(Measurement (24.0,36.11,3.2));


              LT_IFNg_.push_back (Measurement (119.0,3.1,1.21));
/*64*/        LT_TNFa_.push_back(Measurement (119.0,1.93,0.893));//
//
              LT_Apoptosis_.push_back(Measurement(119.0,43.13,3.86));//
/*66*/                   Proliferation_.push_back(Measurement(119.0,5740.0,4326.0));
/*67*/       num_cells_.push_back(Measurement(24.0,1.0e6,5.0e5));//
/*68*/       num_cells_.push_back(Measurement(119.0,1.e6,5.0e5));//

         }
}
///Empieza

///Termina



Results::Results(const std::vector<Measurement>& myTNF,
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
		 double duration):
    TNF_(myTNF),
    IFN_(myIFN),
    APC_expression_(myAPCexpression),
    NK_expression_(myNKexpression),
    LT_expression_(myLTexpression),
    APC_IFNg_(myAPC_IFNg),
    APC_TNFa_(myAPC_TNFa),
    NK_IFNg_(myNK_IFNg),
    NK_TNFa_(myNK_TNFa),
    LT_IFNg_(myLT_IFNg),
    LT_TNFa_(myLT_TNFa),
    LT_Apoptosis_(myLT_Apoptosis),
    Proliferation_(myProliferation),
    num_cells_(mynum_cells),
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


    if ((!APC_IFNg_.empty())&&duration_<APC_IFNg_[APC_IFNg_.size()-1].Time())
        duration_=APC_IFNg_[APC_IFNg_.size()-1].Time();

    if ((!APC_TNFa_.empty())&&duration_<APC_TNFa_[APC_TNFa_.size()-1].Time())
        duration_=APC_TNFa_[APC_TNFa_.size()-1].Time();

    if ((!NK_IFNg_.empty())&&duration_<NK_IFNg_[NK_IFNg_.size()-1].Time())
        duration_=NK_IFNg_[NK_IFNg_.size()-1].Time();

    if ((!NK_TNFa_.empty())&&duration_<NK_TNFa_[NK_TNFa_.size()-1].Time())
        duration_=NK_TNFa_[NK_TNFa_.size()-1].Time();

    if ((!LT_IFNg_.empty())&&duration_<LT_IFNg_[LT_IFNg_.size()-1].Time())
        duration_=NK_IFNg_[NK_IFNg_.size()-1].Time();

    if ((!LT_TNFa_.empty())&&duration_<LT_TNFa_[LT_TNFa_.size()-1].Time())
        duration_=LT_TNFa_[LT_TNFa_.size()-1].Time();

    if ((!LT_Apoptosis_.empty())&&duration_<LT_Apoptosis_[LT_Apoptosis_.size()-1].Time())
        duration_=LT_Apoptosis_[LT_Apoptosis_.size()-1].Time();

    if ((!Proliferation_.empty())&&duration_<Proliferation_[Proliferation_.size()-1].Time())
        duration_=Proliferation_[Proliferation_.size()-1].Time();

    if ((!num_cells_.empty())&&duration_<num_cells_[num_cells_.size()-1].Time())
        duration_=num_cells_[num_cells_.size()-1].Time();

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

std::vector<double> SumSquare_APC_IFN(const Results one, const Results two)
{
    std::vector<double> ss;
    std::size_t n=0;
    std::size_t i=0;
    std::size_t j=0;

    for (;;)
    {
        if (one.APC_IFNg().empty())
            break;
        if (two.APC_IFNg().empty())
            break;
        if(two.APC_IFNg()[j].Time()==one.APC_IFNg()[i].Time())
        {
            ss.push_back(std::pow(two.APC_IFNg()[j].Measure()-one.APC_IFNg()[i].Measure(),2));
            n++;
        }
        if (one.APC_IFNg()[i].Time()<=two.APC_IFNg()[j].Time())
        {
            i++;
            if (i>=one.APC_IFNg().size())
                break;
        }
        else
        {
            j++;
            if (j>=two.APC_IFNg().size())
                break;
        }

    }

    return ss;

}

std::vector<double> SumSquare_APC_TNF(const Results one, const Results two)
{
    std::vector<double> ss;
    std::size_t n=0;
    std::size_t i=0;
    std::size_t j=0;

    for (;;)
    {
        if (one.APC_TNFa().empty())
            break;
        if (two.APC_TNFa().empty())
            break;
        if(two.APC_TNFa()[j].Time()==one.APC_TNFa()[i].Time())
        {
            ss.push_back(std::pow(two.APC_TNFa()[j].Measure()-one.APC_TNFa()[i].Measure(),2));
            n++;
        }
        if (one.APC_TNFa()[i].Time()<=two.APC_TNFa()[j].Time())
        {
            i++;
            if (i>=one.APC_TNFa().size())
                break;
        }
        else
        {
            j++;
            if (j>=two.APC_TNFa().size())
                break;
        }

    }

    return ss;

}

std::vector<double> SumSquare_NK_IFN(const Results one, const Results two)
{
    std::vector<double> ss;
    std::size_t n=0;
    std::size_t i=0;
    std::size_t j=0;

    for (;;)
    {
        if (one.NK_IFNg().empty())
            break;
        if (two.NK_IFNg().empty())
            break;
        if(two.NK_IFNg()[j].Time()==one.NK_IFNg()[i].Time())
        {
            ss.push_back(std::pow(two.NK_IFNg()[j].Measure()-one.NK_IFNg()[i].Measure(),2));
            n++;
        }
        if (one.NK_IFNg()[i].Time()<=two.NK_IFNg()[j].Time())
        {
            i++;
            if (i>=one.NK_IFNg().size())
                break;
        }
        else
        {
            j++;
            if (j>=two.NK_IFNg().size())
                break;
        }

    }

    return ss;

}

std::vector<double> SumSquare_NK_TNF(const Results one, const Results two)
{
    std::vector<double> ss;
    std::size_t n=0;
    std::size_t i=0;
    std::size_t j=0;

    for (;;)
    {
        if (one.NK_TNFa().empty())
            break;
        if (two.NK_TNFa().empty())
            break;
        if(two.NK_TNFa()[j].Time()==one.NK_TNFa()[i].Time())
        {
            ss.push_back(std::pow(two.NK_TNFa()[j].Measure()-one.NK_TNFa()[i].Measure(),2));
            n++;
        }
        if (one.NK_TNFa()[i].Time()<=two.NK_TNFa()[j].Time())
        {
            i++;
            if (i>=one.NK_TNFa().size())
                break;
        }
        else
        {
            j++;
            if (j>=two.NK_TNFa().size())
                break;
        }

    }

    return ss;

}


std::vector<double> SumSquare_LT_IFN(const Results one, const Results two)
{
    std::vector<double> ss;
    std::size_t n=0;
    std::size_t i=0;
    std::size_t j=0;

    for (;;)
    {
        if (one.LT_IFNg().empty())
            break;
        if (two.LT_IFNg().empty())
            break;
        if(two.LT_IFNg()[j].Time()==one.LT_IFNg()[i].Time())
        {
            ss.push_back(std::pow(two.LT_IFNg()[j].Measure()-one.LT_IFNg()[i].Measure(),2));
            n++;
        }
        if (one.LT_IFNg()[i].Time()<=two.LT_IFNg()[j].Time())
        {
            i++;
            if (i>=one.LT_IFNg().size())
                break;
        }
        else
        {
            j++;
            if (j>=two.LT_IFNg().size())
                break;
        }

    }

    return ss;

}

std::vector<double> SumSquare_LT_TNF(const Results one, const Results two)
{
    std::vector<double> ss;
    std::size_t n=0;
    std::size_t i=0;
    std::size_t j=0;

    for (;;)
    {
        if (one.LT_TNFa().empty())
            break;
        if (two.LT_TNFa().empty())
            break;
        if(two.LT_TNFa()[j].Time()==one.LT_TNFa()[i].Time())
        {
            ss.push_back(std::pow(two.LT_TNFa()[j].Measure()-one.LT_TNFa()[i].Measure(),2));
            n++;
        }
        if (one.LT_TNFa()[i].Time()<=two.LT_TNFa()[j].Time())
        {
            i++;
            if (i>=one.LT_TNFa().size())
                break;
        }
        else
        {
            j++;
            if (j>=two.LT_TNFa().size())
                break;
        }

    }

    return ss;

}



std::vector<double> SumSquare_LT_Apoptosis(const Results one, const Results two)
{
    std::vector<double> ss;
    std::size_t n=0;
    std::size_t i=0;
    std::size_t j=0;

    for (;;)
    {
        if (one.LT_Apoptosis().empty())
            break;
        if (two.LT_Apoptosis().empty())
            break;
        if(two.LT_Apoptosis()[j].Time()==one.LT_Apoptosis()[i].Time())
        {
            ss.push_back(std::pow(two.LT_Apoptosis()[j].Measure()-one.LT_Apoptosis()[i].Measure(),2));
            n++;
        }
        if (one.LT_Apoptosis()[i].Time()<=two.LT_Apoptosis()[j].Time())
        {
            i++;
            if (i>=one.LT_Apoptosis().size())
                break;
        }
        else
        {
            j++;
            if (j>=two.LT_Apoptosis().size())
                break;
        }

    }

    return ss;

}


std::vector<double> SumSquare_Proliferation(const Results one, const Results two)
{
    std::vector<double> ss;
    std::size_t n=0;
    std::size_t i=0;
    std::size_t j=0;

    for (;;)
    {
        if (one.Proliferation().empty())
            break;
        if (two.Proliferation().empty())
            break;
        if(two.Proliferation()[j].Time()==one.Proliferation()[i].Time())
        {
            ss.push_back(std::pow(two.Proliferation()[j].Measure()-one.Proliferation()[i].Measure(),2));
            n++;
        }
        if (one.Proliferation()[i].Time()<=two.Proliferation()[j].Time())
        {
            i++;
            if (i>=one.Proliferation().size())
                break;
        }
        else
        {
            j++;
            if (j>=two.Proliferation().size())
                break;
        }

    }

    return ss;

}




std::vector<double> SumSquare_num_cells(const Results one, const Results two)
{
    std::vector<double> ss;
    std::size_t n=0;
    std::size_t i=0;
    std::size_t j=0;

    for (;;)
    {
        if (one.num_cells().empty())
            break;
        if (two.num_cells().empty())
            break;
        if(two.num_cells()[j].Time()==one.num_cells()[i].Time())
        {
            ss.push_back(std::pow(two.num_cells()[j].Measure()-one.num_cells()[i].Measure(),2));
            n++;
        }
        if (one.num_cells()[i].Time()<=two.num_cells()[j].Time())
        {
            i++;
            if (i>=one.num_cells().size())
                break;
        }
        else
        {
            j++;
            if (j>=two.num_cells().size())
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

    ss_other=SumSquare_APC_IFN(one,two);
    ss.insert(ss.end(),ss_other.begin(),ss_other.end());

    ss_other=SumSquare_APC_TNF(one,two);
    ss.insert(ss.end(),ss_other.begin(),ss_other.end());

    ss_other=SumSquare_NK_IFN(one,two);
    ss.insert(ss.end(),ss_other.begin(),ss_other.end());

    ss_other=SumSquare_NK_TNF(one,two);
    ss.insert(ss.end(),ss_other.begin(),ss_other.end());

    ss_other=SumSquare_LT_IFN(one,two);
    ss.insert(ss.end(),ss_other.begin(),ss_other.end());

    ss_other=SumSquare_LT_TNF(one,two);
    ss.insert(ss.end(),ss_other.begin(),ss_other.end());

    ss_other=SumSquare_LT_Apoptosis(one,two);
    ss.insert(ss.end(),ss_other.begin(),ss_other.end());

    ss_other=SumSquare_Proliferation(one,two);
    ss.insert(ss.end(),ss_other.begin(),ss_other.end());

    ss_other=SumSquare_num_cells(one,two);
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
        s<<res.TNF()[i].Time()<<"\t"<<
           res.TNF()[i].Measure()<<"\t"<<
           res.TNF()[i].StdError()<<"\n";

    s<<"IFN\n";
    for (std::size_t i=0; i<res.IFN().size(); i++)
        s<<res.IFN()[i].Time()<<"\t"<<
        res.IFN()[i].Measure()<<"\t"<<
           res.IFN()[i].StdError()<<"\n";

    s<<"APC_expression \n";
    for (std::size_t i=0; i<res.APC_expression().size(); i++)
        s<<res.APC_expression()[i].Time()<<"\t"<<
           res.APC_expression()[i].Measure()<<"\t"<<
           res.APC_expression()[i].StdError()<<"\n";

    s<<"NK_expression \n";
    for (std::size_t i=0; i<res.NK_expression().size(); i++)
        s<<res.NK_expression()[i].Time()<<"\t"<<
           res.NK_expression()[i].Measure()<<"\t"<<
           res.NK_expression()[i].StdError()<<"\n";

    s<<"LT_expression \n";
    for (std::size_t i=0; i<res.LT_expression().size(); i++)
        s<<res.LT_expression()[i].Time()<<"\t"<<
           res.LT_expression()[i].Measure()<<"\t"<<
           res.LT_expression()[i].StdError()<<"\n";
    s<<"APC_IFNg \n";
    for (std::size_t i=0; i<res.APC_IFNg().size(); i++)
        s<<res.APC_IFNg()[i].Time()<<"\t"<<
           res.APC_IFNg()[i].Measure()<<"\t"<<
           res.APC_IFNg()[i].StdError()<<"\n";

    s<<"APC_TNFa \n";
    for (std::size_t i=0; i<res.APC_TNFa().size(); i++)
        s<<res.APC_TNFa()[i].Time()<<"\t"<<
           res.APC_TNFa()[i].Measure()<<"\t"<<
    res.APC_TNFa()[i].StdError()<<"\n";

    s<<"NK_IFNg \n";
    for (std::size_t i=0; i<res.NK_IFNg().size(); i++)
        s<<res.NK_IFNg()[i].Time()<<"\t"<<
           res.NK_IFNg()[i].Measure()<<"\t"<<
    res.NK_IFNg()[i].StdError()<<"\n";

    s<<"NK_TNFa \n";
    for (std::size_t i=0; i<res.NK_TNFa().size(); i++)
        s<<res.NK_TNFa()[i].Time()<<"\t"<<
           res.NK_TNFa()[i].Measure()<<"\t"<<
           res.NK_TNFa()[i].StdError()<<"\n";

    s<<"LT_IFNg \n";
    for (std::size_t i=0; i<res.LT_IFNg().size(); i++)
        s<<res.LT_IFNg()[i].Time()<<"\t"<<
           res.LT_IFNg()[i].Measure()<<"\t"<<
    res.LT_IFNg()[i].StdError()<<"\n";

    s<<"LT_TNFa \n";
    for (std::size_t i=0; i<res.LT_TNFa().size(); i++)
        s<<res.LT_TNFa()[i].Time()<<"\t"<<
           res.LT_TNFa()[i].Measure()<<"\t"<<
    res.LT_TNFa()[i].StdError()<<"\n";

    s<<"LT_Apoptosis \n";
    for (std::size_t i=0; i<res.LT_Apoptosis().size(); i++)
        s<<res.LT_Apoptosis()[i].Time()<<"\t"<<
           res.LT_Apoptosis()[i].Measure()<<"\t"<<
    res.LT_Apoptosis()[i].StdError()<<"\n";

    s<<"Proliferation \n";
    for (std::size_t i=0; i<res.Proliferation().size(); i++)
        s<<res.Proliferation()[i].Time()<<"\t"<<
           res.Proliferation()[i].Measure()<<"\t"<<
    res.Proliferation()[i].StdError()<<"\n";

    s<<"num_cells \n";
    for (std::size_t i=0; i<res.num_cells().size(); i++)
        s<<res.num_cells()[i].Time()<<"\t"<<
           res.num_cells()[i].Measure()<<"\t"<<
    res.num_cells()[i].StdError()<<"\n";



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

    for (std::size_t i=0; i<APC_IFNg_.size(); ++i)
        data.push_back(APC_IFNg_[i].Measure());

    for (std::size_t i=0; i<APC_TNFa_.size(); ++i)
        data.push_back(APC_TNFa_[i].Measure());


    for (std::size_t i=0; i<NK_IFNg_.size(); ++i)
        data.push_back(NK_IFNg_[i].Measure());

    for (std::size_t i=0; i<NK_TNFa_.size(); ++i)
        data.push_back(NK_TNFa_[i].Measure());


    for (std::size_t i=0; i<LT_IFNg_.size(); ++i)
        data.push_back(LT_IFNg_[i].Measure());

    for (std::size_t i=0; i<LT_TNFa_.size(); ++i)
        data.push_back(LT_TNFa_[i].Measure());

    for (std::size_t i=0; i<LT_Apoptosis_.size(); ++i)
        data.push_back(LT_Apoptosis_[i].Measure());

    for (std::size_t i=0; i<Proliferation_.size(); ++i)
        data.push_back(Proliferation_[i].Measure());

    for (std::size_t i=0; i<num_cells_.size(); ++i)
        data.push_back(num_cells_[i].Measure());

    std::vector<double> o(data);
     return o;
 }



std::vector<double> Results::getDataStandardError()const
{
    std::vector<double> data;
            data.clear();
    for (std::size_t i=0; i<TNF_.size(); ++i)
        data.push_back(TNF_[i].StdError());

    for (std::size_t i=0; i<IFN_.size(); ++i)
        data.push_back(IFN_[i].StdError());

    for (std::size_t i=0; i<APC_expression_.size(); ++i)
        data.push_back(APC_expression_[i].StdError());

    for (std::size_t i=0; i<NK_expression_.size(); ++i)
        data.push_back(NK_expression_[i].StdError());

    for (std::size_t i=0; i<LT_expression_.size(); ++i)
        data.push_back(LT_expression_[i].StdError());

    for (std::size_t i=0; i<APC_IFNg_.size(); ++i)
        data.push_back(APC_IFNg_[i].StdError());

    for (std::size_t i=0; i<APC_TNFa_.size(); ++i)
        data.push_back(APC_TNFa_[i].StdError());


    for (std::size_t i=0; i<NK_IFNg_.size(); ++i)
        data.push_back(NK_IFNg_[i].StdError());

    for (std::size_t i=0; i<NK_TNFa_.size(); ++i)
        data.push_back(NK_TNFa_[i].StdError());


    for (std::size_t i=0; i<LT_IFNg_.size(); ++i)
        data.push_back(LT_IFNg_[i].StdError());

    for (std::size_t i=0; i<LT_TNFa_.size(); ++i)
        data.push_back(LT_TNFa_[i].StdError());

    for (std::size_t i=0; i<LT_Apoptosis_.size(); ++i)
        data.push_back(LT_Apoptosis_[i].StdError());

    for (std::size_t i=0; i<Proliferation_.size(); ++i)
        data.push_back(Proliferation_[i].StdError());

    for (std::size_t i=0; i<num_cells_.size(); ++i)
        data.push_back(num_cells_[i].StdError());



     return data;
 }


Results::Results(const Results& other):
    TNF_(other.TNF_),
    IFN_(other.IFN_),

    APC_expression_(other.APC_expression_),
    NK_expression_(other.NK_expression_),
    LT_expression_(other.LT_expression_),

    APC_IFNg_(other.APC_IFNg_),
    APC_TNFa_(other.APC_TNFa_),

    NK_IFNg_(other.NK_IFNg_),
    NK_TNFa_(other.NK_TNFa_),

    LT_IFNg_(other.LT_IFNg_),
    LT_TNFa_(other.LT_TNFa_),
    LT_Apoptosis_(other.LT_Apoptosis_),

    Proliferation_(other.Proliferation_),
    num_cells_(other.num_cells_),
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

    std::swap(one.APC_IFNg_,other.APC_IFNg_);
    std::swap(one.APC_TNFa_,other.APC_TNFa_);

    std::swap(one.NK_IFNg_,other.NK_IFNg_);
    std::swap(one.NK_TNFa_,other.NK_TNFa_);

    std::swap(one.LT_IFNg_,other.LT_IFNg_);
    std::swap(one.LT_TNFa_,other.LT_TNFa_);
    std::swap(one.LT_Apoptosis_,other.LT_Apoptosis_);

    std::swap(one.Proliferation_,other.Proliferation_);
    std::swap(one.num_cells_,other.num_cells_);
    std::swap(one.duration_,other.duration_);
}
