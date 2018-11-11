#ifndef GILLESPIE_H
#define GILLESPIE_H

#include <gsl/gsl_rng.h>
#include <vector>

struct GEvent {
  double deltaT;
  int which;
};

GEvent gillespie_ran_event(gsl_rng* rng, std::vector<double> event_rates);

#endif