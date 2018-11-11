#include <vector>
#include <iostream>

#include "Gillespie.h"

using namespace std;

gsl_rng* initialize_rng(int seed = 0) {
  gsl_rng_env_setup(); // sets default rng from environment
  const gsl_rng_type* rng_type;
  rng_type = gsl_rng_default;
  gsl_rng* rng;
  rng = gsl_rng_alloc(rng_type);
  gsl_rng_set(rng, 0);
  return rng;
}

int main(){
    auto rng = initialize_rng();

    vector<double> weights = { 1.0, 2.0, 3.0, 4.0 };
    GEvent ge;
    for (int i=0; i < 1000; i++) {
      ge = gillespie_ran_event(rng, weights);
      std::cout << ge.deltaT << ", " << ge.which << '\n';
    }

    gsl_rng_free(rng);

    return 0;
}
