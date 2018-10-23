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
#include <gsl/gsl_math.h>
#include <assert.h>
#include <sstream>
#include <map>
#include <sys/time.h>
#include <algorithm>

using namespace std;

enum StateType {
  S_STATE,
  I1_STATE,
  R_STATE,
  P_STATE,
  IR_STATE,
  NUM_OF_STATE_TYPES
};

struct Params{
    double recovery;
    double beta;
    double birth;
    double death;
    double kappa;
    double rho;
    vector<double> Population;
};

const double KAPPA                 = 0.4179; //waning depth parameter
const double RHO                   = 0.2; //waning speed parameter

//other parameters
const vector<double> village_pop   = {1000};
const int numVillages              = village_pop.size(); //total number of villages under consideration
const int numDaysToRecover         = 28;
const double RECOVERY              = 365/numDaysToRecover;    //recovery rate (/year)
const double BETA                  = 135;   //contact rate (individuals/year)
const double lifespan              = 50;
const double BIRTH                 = 1/lifespan; //birth rate (per year)
const double DEATH                 = 1/lifespan; //death rate (per year)

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

vector<double> initialize_compartment(int villageId) {
    //initial population from equilibrium values
    Params params ={};
    params.recovery = RECOVERY;
    params.beta = BETA;
    params.birth = BIRTH;
    params.death = DEATH;
    params.kappa = KAPPA;
    params.rho = RHO;
    params.Population = {village_pop[villageId]};

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
        gsl_vector_set(x,i,100);
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
        fprintf(stderr, "%3d %25.17e\n", i, gsl_vector_get(workspace_F->x, i));
    }

    assert(num_dimensions==5);
    vector<double> compartments = {gsl_vector_get(workspace_F->x,S_STATE),
        gsl_vector_get(workspace_F->x,I1_STATE),
        gsl_vector_get(workspace_F->x,R_STATE),
        gsl_vector_get(workspace_F->x,P_STATE),
        gsl_vector_get(workspace_F->x,IR_STATE)};

    gsl_multiroot_fsolver_free(workspace_F);

    gsl_vector_free(x);

    return compartments;
}

int main(){

    //initialize size of vector of vector of compartments
    vector<vector<double>> compartments(numVillages);

    //find expected compartment size for each village
    for(int i = 0; i < numVillages; i++){
        compartments[i] = initialize_compartment(i);
    }

    return 0;
}
