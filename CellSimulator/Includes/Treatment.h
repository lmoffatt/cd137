#ifndef TREATMENT_H
#define TREATMENT_H


#include <vector>
#include <string>

#include "Measurement.h"

struct Treatment
{


    double Ag;
    double Ab;

    double sim_duration_d;

    double init_cells;

    double time_step_d;

};

#endif // TREATMENT_H
