#ifndef MEASUREMENT_H
#define MEASUREMENT_H

class Measurement
{
public:
    double Time()const;
    double Measure()const;

    void setMeasurement(double newMeasure);

    Measurement(double Time,double Measure);
    Measurement(double Time);
    //copia
    Measurement(const Measurement& other);

    Measurement& operator=(const Measurement& other);

    Measurement();

    virtual ~Measurement(){}

    void friend swap(Measurement& one, Measurement& two);
private:
    double time_;
    double measure_;
};

#endif // MEASUREMENT_H
