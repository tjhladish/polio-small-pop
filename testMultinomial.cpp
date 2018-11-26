#include <vector>
#include <string>
#include <gsl/gsl_rng.h>

#include "Params.h"
#include "States.h"

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
    string jsfile = "testeq.json";
    Params refparams = parseParams(jsfile);
    //int numVillages = refparams.Population.size();
    auto rng = initialize_rng();

    auto proportions = equilibrium_fraction(refparams);
    vector<vector<unsigned int>> compartments(refparams.Population.size());
    for (int i=0; i < compartments.size(); i++) {
      compartments[i] = multinomial_compartments(rng, proportions, refparams.Population[i]);
      printDiscretePop(compartments[i]);
    }

    gsl_rng_free(rng);

    return 0;
}
