#include "Measurement.h"
#

double Measurement::Time()const
{
    return time_;
}

double Measurement::Measure()const
{
    return measure_;
}

Measurement::Measurement(double Time,double Measure):
    time_(Time),
    measure_(Measure){}

Measurement::Measurement(double Time):
    time_(Time),
    measure_(0.0/0.0)
{}


void Measurement::setMeasurement(double newMeasure)
{
    measure_=newMeasure;
}

Measurement::Measurement(const Measurement& other):
    time_(other.time_),
    measure_(other.measure_){}

/*Measurement& operator=(const Measurement& other)
{
    if (this!=&other)
    {
        Measurement tmp(other);
        swap(*this,tmp);
    }
    return *this;
}*/

//Measurement();

/*void friend swap(Measurement& one, Measurement& two)
{
   std::swap(one.time_,two.time_);
   std::swap(one.measure_,two.measure_);
}
*/
