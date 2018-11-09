#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>

#include "States.h"

std::vector<double> equilibrium_fraction(Params p) {
  
}

// MULTIROOT CODE

std::vector<double> multiroot_solver(Params p) {
  
}

int multiroot_evolve(const gsl_vector * x, void * p, gsl_vector * f){
    Params * params = (Params *)p;
    const double recovery = (params->recovery);
    const double beta     = (params->beta);
    const double birth    = (params->birth);
    const double death    = (params->death);
    const double kappa    = (params->kappa);
    const double rho      = (params->rho);
    const double Population = 1.0; // (params->Population[0]); // assume const population
    const double S  = gsl_vector_get(x,S_STATE);
    const double I1 = gsl_vector_get(x,I1_STATE);
    const double R  = gsl_vector_get(x,R_STATE);
    const double P  = gsl_vector_get(x,P_STATE);
    const double Ir = gsl_vector_get(x,IR_STATE);
    // assert sum(S,I1,R,P,Ir)

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

// UNIROOT CODE

std::vector<double> uniroot_solver(Params p) {
  
}