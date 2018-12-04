#include "Gillespie.h"

template int discrete_ran_weighted(gsl_rng* rng, std::vector<int> weights);
template int discrete_ran_weighted(gsl_rng* rng, std::vector<unsigned int> weights);
template int discrete_ran_weighted(gsl_rng* rng, std::vector<double> weights);

GEvent gillespie_ran_event(gsl_rng* rng, std::vector<double> event_rates) {
  double overall_rate = std::accumulate(event_rates.begin(), event_rates.end(), 0);
  
  // time is an exponential draw
  double t = gsl_ran_exponential(rng, 1.0/overall_rate); // the 1.0/sm converts gsl parameterization (mean waiting time) to rate-based 
  int which = discrete_ran_weighted(rng, event_rates);

  return { t, which };
}

bool gillespie_ran_bool(gsl_rng* rng, double ptrue) {
  return gsl_ran_flat(rng, 0.0, 1) < ptrue;
}