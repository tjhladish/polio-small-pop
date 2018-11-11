#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <vector>
#include <algorithm>
#include <numeric>

#include "Gillespie.h"

GEvent gillespie_ran_event(gsl_rng* rng, std::vector<double> event_rates) {
  double overall_rate = std::accumulate(event_rates.begin(), event_rates.end(), 0);
  
  // time is an exponential draw
  double t = gsl_ran_exponential(rng, 1.0/overall_rate); // the 1.0/sm converts gsl parameterization (mean waiting time) to rate-based 
  
  // which event is a proportional draw, so need cdf by event
  // this creates the cumulative sum over event rates
  std::partial_sum(event_rates.begin(), event_rates.end(), event_rates.begin());
  std::vector<double>::iterator up =
    std::upper_bound(
      event_rates.begin(), event_rates.end(),
      gsl_ran_flat(rng, 0.0, overall_rate)
  );
  int which = (up - event_rates.begin());

  return { t, which };
}