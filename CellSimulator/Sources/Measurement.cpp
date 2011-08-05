#include "Measurement.h"

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

