#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <vector>
#include <algorithm>
#include <numeric>

#include "Gillespie.h"

int discrete_ran_weighted(gsl_rng* rng, std::vector<double> weights) {
  std::partial_sum(weights.begin(), weights.end(), weights.begin());
  std::vector<double>::iterator up =
    std::upper_bound(
      weights.begin(), weights.end(),
      gsl_ran_flat(rng, 0.0, weights.back())
  );
  return (up - weights.begin());
}

GEvent gillespie_ran_event(gsl_rng* rng, std::vector<double> event_rates) {
  double overall_rate = std::accumulate(event_rates.begin(), event_rates.end(), 0);
  
  // time is an exponential draw
  double t = gsl_ran_exponential(rng, 1.0/overall_rate); // the 1.0/sm converts gsl parameterization (mean waiting time) to rate-based 
  int which = discrete_ran_weighted(rng, event_rates);

  return { t, which };
}