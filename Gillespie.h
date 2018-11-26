#ifndef GILLESPIE_H
#define GILLESPIE_H

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <vector>
#include <numeric>
#include <algorithm>

struct GEvent {
  double deltaT;
  int which;
};

template <typename NUMERIC>
int discrete_ran_weighted(gsl_rng* rng, NUMERIC * weights, size_t size) {
  std::vector<NUMERIC> tmp(size);
  std::partial_sum(weights, weights + size, tmp.begin());
  auto up = std::upper_bound(
    tmp.begin(), tmp.end(),
    gsl_ran_flat(rng, 0.0, tmp.back())
  );
  return (up - tmp.begin());
}

template <typename NUMERIC>
int discrete_ran_weighted(gsl_rng* rng, std::vector<NUMERIC> weights) {
  return discrete_ran_weighted(rng, &weights[0], weights.size());
}

// TODO also a begin / end interface?

GEvent gillespie_ran_event(gsl_rng* rng, std::vector<double> event_rates);
bool   gillespie_ran_bool(gsl_rng* rng, double ptrue);

#endif