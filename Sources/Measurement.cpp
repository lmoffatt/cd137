#include "Measurement.h"
#include <algorithm>

double Measurement::Time()const
{
    return time_;
}

double Measurement::Measure()const
{
    return measure_;
}

double Measurement::StdError()const
{
    return stdError_;
}


Measurement::~Measurement(){}



Measurement::Measurement(double Time,double Measure,double Std
                         ):
    time_(Time),
    measure_(Measure),
    stdError_(Std){}

Measurement::Measurement(double Time):
    time_(Time),
    measure_(0.0/0.0)
{}


Measurement& Measurement::setMeasurement(double newMeasure)
{
    measure_=newMeasure;
    return *this;
}

Measurement::Measurement(const Measurement& other):
    time_(other.time_),
    measure_(other.measure_),
    stdError_(other.stdError_){}

Measurement& Measurement::operator=(const Measurement& other)
{
    if (this!=&other)
    {
        Measurement tmp(other);
        swap(*this,tmp);
    }
    return *this;
}

Measurement::Measurement():
    time_(),
    measure_(),
    stdError_(){}

void  swap(Measurement& one, Measurement& two)
{
   std::swap(one.time_,two.time_);
   std::swap(one.measure_,two.measure_);
   std::swap(one.stdError_,two.stdError_);
}

