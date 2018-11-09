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

// const double KAPPA                 = 0.4179; //waning depth parameter
// const double RHO                   = 0.2; //waning speed parameter

//other parameters
// const vector<double> village_pop   = {1000};
// const int numVillages              = village_pop.size(); //total number of villages under consideration
// const int numDaysToRecover         = 28;
// const double RECOVERY              = 365/numDaysToRecover;    //recovery rate (/year)
// const double BETA                  = 135;   //contact rate (individuals/year)
// const double lifespan              = 50;
// const double BIRTH                 = 1/lifespan; //birth rate (per year)
// const double DEATH                 = 1/lifespan; //death rate (per year)

int func_m(const gsl_vector * x, void * p, gsl_vector * f){
    Params * params = (Params *)p;
    const double recovery = (params->recovery);
    const double beta = (params->beta);
    const double birth = (params->birth);
    const double death = (params->death);
    const double kappa = (params->kappa);
    const double rho = (params->rho);
    const double Population = (params->Population[0]);
    const double S  = gsl_vector_get(x,S_STATE);
    const double I1 = gsl_vector_get(x,I1_STATE);
    const double R  = gsl_vector_get(x,R_STATE);
    const double P  = gsl_vector_get(x,P_STATE);
    const double Ir = gsl_vector_get(x,IR_STATE);

    const double birthdt = birth*Population;
    const double foi = ((beta*(I1+kappa*Ir))/Population);
    const double infection1 = S*foi;
    const double infection2 = P*foi*kappa;
    const double recovery1 = recovery*I1;
    const double recovery2 = (recovery/kappa)*Ir;
    const double waning = rho*R;

    gsl_vector_set (f, S_STATE,  birthdt - infection1 - death*S);
    gsl_vector_set (f, I1_STATE, infection1 - recovery1 - death*I1);
    gsl_vector_set (f, R_STATE,  recovery1 + recovery2 - waning - death*R);
    gsl_vector_set (f, P_STATE,  waning - infection2 - death*P);
    gsl_vector_set (f, IR_STATE, infection2 - recovery2 - death*Ir);

    return GSL_SUCCESS;
}

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

    int i, times, status;
    gsl_multiroot_function F;
    gsl_multiroot_fsolver *workspace_F;
    gsl_vector *x;
    int num_dimensions = 5;

    x = gsl_vector_alloc(num_dimensions);

    workspace_F = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrids,num_dimensions);
    printf("F solver: %s\n", gsl_multiroot_fsolver_name(workspace_F));
    F.f=&func_m;
    F.n=num_dimensions;
    F.params = &params;
    int MAXTIMES = 100;
    /* set initial value */
    for(i = 0; i < num_dimensions; i++){
        gsl_vector_set(x, i, params.Population[0]/num_dimensions);
    }

    /* set solver */
    gsl_multiroot_fsolver_set(workspace_F,&F, x);

    /* main loop */
    for(times = 0; times < MAXTIMES; times++) {
        status = gsl_multiroot_fsolver_iterate(workspace_F);

        /*fprintf(stderr, "%d times: ", times);
         for(i = 0; i < num_dimensions; i++) {
         fprintf(stderr, "%10.3e ", gsl_vector_get(workspace_F->x, i));
         }
         fprintf(stderr, "\n");
         */
        if((status == GSL_EBADFUNC) || (status == GSL_ENOPROG))
        {
            fprintf(stderr, "Status: %s\n", gsl_strerror(status));
            break;
        }
    }

    assert(status == GSL_ENOPROG);
    /* print answer */
    for(i = 0; i < num_dimensions; i++) {
        fprintf(stderr, "%s %25.17e\n", statestr[i].c_str(), gsl_vector_get(workspace_F->x, i));
    }

    assert(num_dimensions==5);
    vector<double> compartments = {
      gsl_vector_get(workspace_F->x,S_STATE),
      gsl_vector_get(workspace_F->x,I1_STATE),
      gsl_vector_get(workspace_F->x,R_STATE),
      gsl_vector_get(workspace_F->x,P_STATE),
      gsl_vector_get(workspace_F->x,IR_STATE)
    };

    gsl_multiroot_fsolver_free(workspace_F);

    gsl_vector_free(x);

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

    /* print answer */
    fprintf(stderr, "I1: %25.17e\n", I1);
    fprintf(stderr, "S: %25.17e\n", S);
    fprintf(stderr, "IR: %25.17e\n", IR);
    fprintf(stderr, "R: %25.17e\n", R);
    fprintf(stderr, "P: %25.17e\n", P);
    
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
