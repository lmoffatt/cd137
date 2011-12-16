#ifndef TREATMENT_H
#define TREATMENT_H


#include <vector>
#include <string>
#include <iostream>

#include "Measurement.h"

struct Treatment
{


    double Ag;
    double Ab;

    double sim_duration_d;

    double init_cells;

    double time_step_d;

    double t_apop_meas_d;

};

inline std::ostream& operator<<(std::ostream& s, const Treatment& tr)
{
    s<<"\n Ag \t"<<tr.Ag;
    s<<"\n Ab \t"<<tr.Ab;
    s<<"\n sim_duration_d \t"<<tr.sim_duration_d;
    s<<"\n init_cells \t"<<tr.init_cells;
    s<<"\n time_step_d \t"<<tr.time_step_d;
return s;
}

#endif // TREATMENT_H
