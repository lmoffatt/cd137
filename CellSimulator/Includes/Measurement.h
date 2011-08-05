#ifndef MEASUREMENT_H
#define MEASUREMENT_H

class Measurement
{
public:
    double Time()const;
    double Measure()const;

    Measurement(double Time,double Measure);
private:
    double time_;
    double measure_;
};

#endif // MEASUREMENT_H
