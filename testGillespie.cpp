#include <vector>
#include <iostream>

#include "Gillespie.h"
#include "Utils.h"

using namespace std;

int main(){
    auto rng = initialize_rng();

    vector<double> weights = { 1.0, 2.0, 3.0, 4.0 };
    GEvent ge;
    for (int i=0; i < 1000; i++) {
      ge = gillespie_ran_event(rng, weights);
      std::cout << ge.deltaT << ", " << ge.which << '\n';
    }

    vector<unsigned int> iweights = {4, 3, 2, 1};
    int which;
    for (int i=0; i < 1000; i++) {
      which = discrete_ran_weighted(rng, &iweights[1], iweights.size()-1)+1;
      std::cout << which << '\n';
    }

    gsl_rng_free(rng);

    return 0;
}
