#ifndef PARAMS_H
#define PARAMS_H

#include <vector>
#include <string>

struct Params {
    double recovery;  // rate; I1->R & IR->R
    double beta;      // infection rate per capita; S->I1, P->IR
    double birth;     // birth rate; equivalent to death in const pop model
    double death;     // rate; all classes -> S (via birth)
    double kappa;     // scaling parameter for pre-exposed class
    double rho;       // rate; waning R->P
    std::vector<double> Population;
    // TODO: sort of weird that population appears as part of Parameters
    // or more particularly that it's a vector instead of single value
};

Params parseParams(std::string jsonfile);

#endif