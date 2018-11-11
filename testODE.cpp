#include <vector>
#include <string>

#include "Params.h"
#include "States.h"

using namespace std;

vector<double> initialize_compartment(int villageId, Params ref) {
    //initial population from equilibrium values
    Params params = ref;
    params.Population = {ref.Population[villageId]};

    vector<double> compartments = equilibrium_fraction(ref);

    printResults(compartments);

    return compartments;
}

vector<double> findEquilibrium(Params params) {
    auto res = equilibrium_fraction(params, false);
    printResults(res);
    return res;
}

int main(){
    string jsfile = "testeq.json";
    Params refparams = parseParams(jsfile);
    int numVillages = refparams.Population.size();

    //initialize size of vector of vector of compartments
    vector<vector<double>> compartments(numVillages);

    //find expected compartment size for each village
    for(int i = 0; i < numVillages; i++){
        compartments[i] = initialize_compartment(i, refparams);
    }

    findEquilibrium(refparams);

    return 0;
}
