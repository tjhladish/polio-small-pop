#include <iostream>
#include <math.h>
#include <random>
#include <tuple>
#include <array>
#include <vector>
#include <chrono>
#include <fstream>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>
#include <assert.h>
#include <sstream>
#include <map>
#include <sys/time.h>
#include <algorithm>

#include "Params.h"
#include "States.h"

using namespace std;

double func_s(double I1, Params* p) {
  return 1.0 - ((p->recovery + p->death)/(p->death))*I1;
}

double func_ir(double I1, double S, Params* p) {
  return (((p->recovery+p->death)/(p->beta * S)) - 1.0)*(I1 / p->kappa);
}

double func_r(double I1, double IR, Params* p) {
  return (p->recovery)/(p->rho + p->death)*(I1 + IR/p->kappa);
}

double func_p(double I1, double IR, Params* p) {
  return ((p->recovery)/(p->kappa) + p->death) * IR/(p->kappa * p->beta * (I1 + p->kappa * IR));
}

double func_i1(double I1, void * p){
    Params * ps = (Params *)p;

    double S = func_s(I1, ps);
    double IR = func_ir(I1, S, ps);
    double R = func_r(I1, IR, ps);
    double P = func_p(I1, IR, ps);

    return 1 - S - I1 - IR - R - P;
}

vector<double> initialize_compartment(int villageId, Params ref) {
    //initial population from equilibrium values
    Params params = ref;
    params.Population = {ref.Population[villageId]};

    vector<double> compartments = equilibrium_fraction(ref);

    printResults(compartments);

    return compartments;
}

//vector<double>
void findEquilibrium(Params params) {

    gsl_function F;
    F.function = &func_i1;
    F.params = &params;
    
    gsl_root_fsolver * workspace_F = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    printf("F solver: %s\n", gsl_root_fsolver_name(workspace_F));

    double i1lo = params.death*(1.0/(params.death+params.recovery) - 1.0/params.beta);
    double i1hi = (params.death / (params.recovery + params.death)) - 1e-8;

    /* set solver */
    gsl_root_fsolver_set(workspace_F, &F, i1lo, i1hi);

    int status;
    int iter = 0;
    int max_iter = 100; // should be an arg to solving function
    double I1;
    
    /* main loop */
    do {
      iter++;
      status = gsl_root_fsolver_iterate(workspace_F);
      I1 = gsl_root_fsolver_root (workspace_F);
      i1lo = gsl_root_fsolver_x_lower (workspace_F);
      i1hi = gsl_root_fsolver_x_upper (workspace_F);
      status = gsl_root_test_interval (i1lo, i1hi, 0, 0.001);
    } while (status == GSL_CONTINUE && iter < max_iter);

    // assert(status == GSL_ENOPROG);
    
    // double I1 = gsl_root_fsolver_root(workspace_F);
    double S = func_s(I1, &params);
    double IR = func_ir(I1, S, &params);
    double R = func_r(I1, IR, &params);
    double P = func_p(I1, IR, &params);

    vector<double> res = { S, I1, R, P, IR };

    printResults(res);

    gsl_root_fsolver_free(workspace_F);

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
