#ifndef MEASUREMENT_H
#define MEASUREMENT_H


class Measurement
{
public:

    double Time()const;
    double Measure()const;
    double StdError()const;

    Measurement& setMeasurement(double newMeasure);

    Measurement(double Time,double Measure,double Std);
    Measurement(double Time);
    //copia
    Measurement(const Measurement& other);

    Measurement& operator=(const Measurement& other);

    Measurement();

    ~Measurement();

    void friend swap(Measurement& one, Measurement& two);
private:
    double time_;
    double measure_;
    double stdError_;
};



#endif // MEASUREMENT_H

